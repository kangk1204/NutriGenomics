import csv
import gzip
import io
import json
import sqlite3
from pathlib import Path
import pytest
from nutriomics_evidence.ctd import ingest_ctd,read_ctd,mesh_id,gene_id
from nutriomics_evidence.graph import build_graph,connect,query,validate
from nutriomics_evidence.audit import audit_ctd

FIELDS={
 'CTD_chem_gene_ixns.csv.gz':'ChemicalName,ChemicalID,CasRN,GeneSymbol,GeneID,GeneForms,Organism,OrganismID,Interaction,InteractionActions,PubMedIDs',
 'CTD_curated_chemicals_diseases.csv.gz':'ChemicalName,ChemicalID,CasRN,DiseaseName,DiseaseID,DirectEvidence,PubMedIDs',
 'CTD_curated_genes_diseases.csv.gz':'GeneSymbol,GeneID,DiseaseName,DiseaseID,DirectEvidence,OmimIDs,PubMedIDs',
 'CTD_chemicals_diseases.csv.gz':'ChemicalName,ChemicalID,CasRN,DiseaseName,DiseaseID,DirectEvidence,InferenceGeneSymbol,InferenceScore,OmimIDs,PubMedIDs',
 'CTD_chemicals.csv.gz':'ChemicalName,ChemicalID,CasRN,PubChemCID,PubChemSID,DTXSID,InChIKey,Definition,ParentIDs,TreeNumbers,ParentTreeNumbers,MESHSynonyms,CTDCuratedSynonyms'}

def write_ctd(path,fields,rows,date='Tue Sep 29 13:22:17 EDT 2026'):
    content=io.StringIO(newline='');content.write('# CTD\n# Report created: '+date+'\n# Fields:\n# '+fields+'\n#\n')
    csv.writer(content,lineterminator='\n').writerows(rows)
    path.write_bytes(gzip.compress(content.getvalue().encode()))

@pytest.fixture
def fixture_dir(tmp_path):
    directory=tmp_path/'source';directory.mkdir()
    write_ctd(directory/'CTD_curated_chemicals_diseases.csv.gz',FIELDS['CTD_curated_chemicals_diseases.csv.gz'],[
        ['X','C123','','Hypertension','MESH:D006973','therapeutic','123'],
        ['X','C123','','Hypertension','MESH:D006973','marker/mechanism','456'],
        ['Y','C456','','Hypertension, Pulmonary','MESH:D006976','therapeutic','321']])
    write_ctd(directory/'CTD_curated_genes_diseases.csv.gz',FIELDS['CTD_curated_genes_diseases.csv.gz'],[
        ['ACE','1636','Hypertension','MESH:D006973','marker/mechanism|therapeutic','','123|456']])
    interaction='X increases ACE expression, and a quoted\n# phrase is preserved'
    write_ctd(directory/'CTD_chem_gene_ixns.csv.gz',FIELDS['CTD_chem_gene_ixns.csv.gz'],[
        ['X','C123','','ACE','1636','protein','Homo sapiens','9606',interaction,'increases^expression','123'],
        ['X','C123','','ACE','1636','protein','Homo sapiens','9606',interaction,'increases^expression','456'],
        ['X','C123','','ACE','1636','protein','Mus musculus','10090',interaction,'increases^expression','789']])
    write_ctd(directory/'CTD_chemicals_diseases.csv.gz',FIELDS['CTD_chemicals_diseases.csv.gz'],[
        ['X','C123','','Hypertension','MESH:D006973','therapeutic','','','','123'],
        ['X','C123','','Hypertension','MESH:D006973','','ACE','12.5','','456']])
    return directory

def test_stream_quoting_metadata_and_crc(fixture_dir):
    data=list(read_ctd(fixture_dir/'CTD_chem_gene_ixns.csv.gz'))
    assert data[0]['report_created'].endswith('2026')
    assert len(data)==4 and '\n# phrase' in data[1]['Interaction']


def test_complete_audit_reports_corruption_without_materializing(fixture_dir,tmp_path):
    output=tmp_path/'audit.json'
    result=audit_ctd(fixture_dir,output)
    assert result['complete'] and result['validated_files']==4
    assert not list(tmp_path.glob('*.parquet'))
    damaged=fixture_dir/'CTD_chem_gene_ixns.csv.gz';payload=bytearray(damaged.read_bytes());payload[-8]^=1;damaged.write_bytes(payload)
    result=audit_ctd(fixture_dir,output)
    assert not result['complete'] and result['validated_files']==3
    assert result['errors'][0]['file']==damaged.name

def test_mesh_and_gene_normalization():
    assert mesh_id('C123')==mesh_id('MESH:C123')=='MESH:C123'
    assert mesh_id('D006973')=='MESH:D006973'
    assert gene_id('1636')=='NCBIGene:1636'
    with pytest.raises(ValueError):mesh_id('Hypertension')
    with pytest.raises(ValueError):gene_id('ACE')


def test_dictionary_hierarchy_roots_not_molecular_entities(fixture_dir,tmp_path):
    filename='CTD_chemicals.csv.gz'
    write_ctd(fixture_dir/filename,FIELDS[filename],[
        ['Descriptor hierarchy','MESH:D','','','','','','','','','','',''],
        ['X','MESH:C123','','123','','','AAAAAAAAAAAAAA-BBBBBBBBBB-C','','','','','','']])
    out=tmp_path/'out';ingest_ctd(fixture_dir,out,[*list(FIELDS)[:3],filename])
    database=tmp_path/'graph.sqlite';build_graph(out,database)
    db=connect(database,readonly=True)
    assert not db.execute("SELECT 1 FROM node WHERE id='MESH:D'").fetchone()
    assert db.execute("SELECT identifier FROM identifier_mapping WHERE subject='MESH:C123' AND namespace='PubChem'").fetchone()[0]=='123'
    db.close()

def test_damaged_gzip_never_publishes(fixture_dir,tmp_path):
    path=fixture_dir/'CTD_chem_gene_ixns.csv.gz';data=bytearray(path.read_bytes());data[-8]^=1;path.write_bytes(data)
    with pytest.raises((gzip.BadGzipFile,EOFError)):
        ingest_ctd(fixture_dir,tmp_path/'out',[path.name])
    assert not (tmp_path/'out'/'CTD_chem_gene_ixns.parquet').exists()
    assert not (tmp_path/'out'/'manifest.json').exists()

def test_missing_schema_column_rejected(tmp_path):
    path=tmp_path/'CTD_chem_gene_ixns.csv.gz';write_ctd(path,'ChemicalName,ChemicalID',[['X','C123']])
    with pytest.raises(ValueError,match='schema'):list(read_ctd(path))

def test_stale_dictionary_date_is_not_download_date(tmp_path):
    path=tmp_path/'CTD_types.csv.gz';write_ctd(path,'TypeName,Code',[['binding','b']],date='Mon Feb 12 14:09:00 EST 2024')
    result=ingest_ctd(tmp_path,tmp_path/'out',[path.name])
    assert result['files'][path.name]['version']=='2024-02-12'

def test_selected_graph_separation_dedup_taxon(fixture_dir,tmp_path):
    out=tmp_path/'parquet';ingest_ctd(fixture_dir,out,list(FIELDS)[:4])
    report=build_graph(out,tmp_path/'graph.sqlite',include_inferred=True)
    assert report['errors']==[]
    db=connect(tmp_path/'graph.sqlite',readonly=True)
    assert not db.execute("SELECT 1 FROM node WHERE id='MESH:D006976'").fetchone()
    influence=db.execute("SELECT * FROM edge WHERE relation='chemical_gene_influence'").fetchall()
    assert len(influence)==1 and influence[0]['taxon']=='9606'
    assert db.execute('SELECT count(*) FROM citation WHERE edge_id=?',(influence[0]['id'],)).fetchone()[0]==2
    therapeutic=db.execute("SELECT id FROM edge WHERE relation='chemical_disease_therapeutic'").fetchone()[0]
    assert db.execute('SELECT count(*) FROM edge_source WHERE edge_id=?',(therapeutic,)).fetchone()[0]==2
    inferred=query(db,tier='inferred');assert len(inferred)==1 and inferred[0]['context']['score_is_not_probability']
    assert not db.execute("SELECT 1 FROM edge WHERE relation IN ('dti_positive','direct_binding')").fetchone()
    db.close()

def test_parquet_tampering_rejected(fixture_dir,tmp_path):
    out=tmp_path/'out';ingest_ctd(fixture_dir,out,list(FIELDS)[:3]);p=out/'CTD_curated_genes_diseases.parquet';p.write_bytes(p.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='hash mismatch'):build_graph(out,tmp_path/'graph.sqlite')

def test_truncated_rows_rejected(tmp_path):
    path=tmp_path/'CTD_types.csv.gz';write_ctd(path,'TypeName,Code',[['binding']])
    with pytest.raises(ValueError,match='expected 2'):list(read_ctd(path))

def test_crdownload_rejected(tmp_path):
    with pytest.raises(ValueError):next(read_ctd(tmp_path/'516865.crdownload'))

def test_replay_is_idempotent_and_new_snapshot_requires_new_dir(fixture_dir,tmp_path):
    out=tmp_path/'out';name='CTD_curated_chemicals_diseases.csv.gz'
    first=ingest_ctd(fixture_dir,out,[name]);second=ingest_ctd(fixture_dir,out,[name]);assert first==second
    path=fixture_dir/name;path.write_bytes(gzip.compress(gzip.decompress(path.read_bytes()).replace(b'123',b'999')))
    with pytest.raises(ValueError,match='new source release'):ingest_ctd(fixture_dir,out,[name])
