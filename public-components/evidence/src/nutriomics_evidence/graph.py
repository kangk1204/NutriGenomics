"""Selected SQLite graph over versioned Parquet, with separate evidence tiers."""
from __future__ import annotations
import hashlib
import json
import os
import sqlite3
from pathlib import Path
from .ctd import gene_id, mesh_id
from .sources import atomic_json, now, sha256

DISEASES = ('MESH:D006973', 'MESH:D000075222')
SCHEMA = '''
PRAGMA foreign_keys=ON;
CREATE TABLE IF NOT EXISTS source(id TEXT PRIMARY KEY, metadata_json TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS node(id TEXT PRIMARY KEY, kind TEXT NOT NULL, name TEXT NOT NULL, metadata_json TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS edge(id TEXT PRIMARY KEY, subject TEXT NOT NULL REFERENCES node(id), object TEXT NOT NULL REFERENCES node(id), relation TEXT NOT NULL, tier TEXT NOT NULL CHECK(tier IN ('curated','inferred','observed_composition','regulatory_identity','predicted','extracted')), taxon TEXT, context_json TEXT NOT NULL);
CREATE TABLE IF NOT EXISTS citation(edge_id TEXT NOT NULL REFERENCES edge(id), pmid TEXT NOT NULL, PRIMARY KEY(edge_id,pmid));
CREATE TABLE IF NOT EXISTS edge_source(edge_id TEXT NOT NULL REFERENCES edge(id),source_id TEXT NOT NULL REFERENCES source(id),PRIMARY KEY(edge_id,source_id));
CREATE TABLE IF NOT EXISTS identifier_mapping(subject TEXT NOT NULL REFERENCES node(id),namespace TEXT NOT NULL,identifier TEXT NOT NULL,relation TEXT NOT NULL CHECK(relation IN ('source_reported','exact_structure','related','candidate_name','regulatory_identity')),source_id TEXT NOT NULL REFERENCES source(id),basis TEXT NOT NULL,PRIMARY KEY(subject,namespace,identifier,source_id));
CREATE TABLE IF NOT EXISTS measurement(id TEXT PRIMARY KEY, food TEXT NOT NULL REFERENCES node(id),nutrient TEXT NOT NULL REFERENCES node(id),amount REAL,unit TEXT,basis TEXT NOT NULL,status TEXT NOT NULL,source_id TEXT NOT NULL REFERENCES source(id),raw_json TEXT NOT NULL);
CREATE INDEX IF NOT EXISTS edge_subject ON edge(subject);
CREATE INDEX IF NOT EXISTS edge_object ON edge(object);
CREATE INDEX IF NOT EXISTS edge_relation ON edge(relation,tier);
'''


def connect(path: Path, readonly: bool = False):
    if readonly:
        db = sqlite3.connect(path.resolve().as_uri() + '?mode=ro', uri=True)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        db = sqlite3.connect(path)
        db.executescript(SCHEMA)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA foreign_keys=ON')
    return db


def canonical(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def node(db, identifier: str, kind: str, name: str, metadata: dict | None = None) -> None:
    existing = db.execute('SELECT kind,metadata_json FROM node WHERE id=?', (identifier,)).fetchone()
    if existing and existing['kind'] != kind:
        raise ValueError('node kind conflict: ' + identifier)
    merged = json.loads(existing['metadata_json']) if existing else {}
    merged.update(metadata or {})
    db.execute('INSERT INTO node VALUES (?,?,?,?) ON CONFLICT(id) DO UPDATE SET name=excluded.name,metadata_json=excluded.metadata_json', (identifier, kind, name, canonical(merged)))


def source(db, record: dict) -> None:
    existing = db.execute('SELECT metadata_json FROM source WHERE id=?', (record['id'],)).fetchone()
    if existing and json.loads(existing[0])['sha256'] != record['sha256']:
        raise ValueError('source identity conflict')
    db.execute('INSERT OR IGNORE INTO source VALUES (?,?)', (record['id'], canonical(record)))


def edge(db, subject: str, target: str, relation: str, tier: str, taxon: str | None, context: dict, source_id: str, pmids: str = '') -> str:
    if relation in {'direct_binding', 'dti_positive'} and context.get('source') == 'CTD':
        raise ValueError('CTD influence edges cannot become direct-binding labels')
    context_text = canonical(context)
    key = canonical([subject, target, relation, tier, taxon, context])
    identifier = hashlib.sha256(key.encode()).hexdigest()
    db.execute('INSERT OR IGNORE INTO edge VALUES (?,?,?,?,?,?,?)', (identifier, subject, target, relation, tier, taxon, context_text))
    db.execute('INSERT OR IGNORE INTO edge_source VALUES (?,?)', (identifier, source_id))
    for pmid in filter(None, pmids.split('|')):
        if not pmid.isdigit():
            raise ValueError('invalid PMID: ' + pmid)
        db.execute('INSERT OR IGNORE INTO citation VALUES (?,?)', (identifier, pmid))
    return identifier


def parquet_rows(path: Path, sql_where: str = '', parameters=()):
    import duckdb
    db = duckdb.connect()
    try:
        cursor = db.execute('SELECT * FROM read_parquet(?) ' + sql_where, [str(path), *parameters])
        names = [r[0] for r in cursor.description]
        while batch := cursor.fetchmany(10000):
            for row in batch:
                yield dict(zip(names, row, strict=True))
    finally:
        db.close()


def build_graph(input_dir: Path, database: Path, overwrite: bool = False, include_inferred: bool = False) -> dict:
    manifest = json.loads((input_dir / 'manifest.json').read_text(encoding='utf-8'))
    files = manifest['files']
    required = ('CTD_curated_chemicals_diseases.csv.gz', 'CTD_curated_genes_diseases.csv.gz', 'CTD_chem_gene_ixns.csv.gz')
    if any(name not in files for name in required):
        raise ValueError('all three validated core CTD artifacts are required')
    if database.exists() and not overwrite:
        raise FileExistsError('graph exists; use --overwrite for an explicit rebuild')
    # Verify artifacts before using them, including relocated directories.
    for record in files.values():
        if sha256(input_dir / record['parquet']) != record['parquet_sha256']:
            raise ValueError('Parquet hash mismatch: ' + record['parquet'])
    temporary = database.with_suffix(database.suffix + '.building')
    temporary.unlink(missing_ok=True)
    db = connect(temporary)
    chemicals, genes = set(), set()
    try:
        with db:
            for record in files.values():
                source(db, record)
            for filename, kind in ((required[0], 'chemical'), (required[1], 'gene')):
                record = files[filename]
                for row in parquet_rows(input_dir / record['parquet'], 'WHERE DiseaseID IN (?,?)', DISEASES):
                    disease = mesh_id(row['DiseaseID'])
                    subject = mesh_id(row['ChemicalID']) if kind == 'chemical' else gene_id(row['GeneID'])
                    (chemicals if kind == 'chemical' else genes).add(subject)
                    node(db, subject, kind, row['ChemicalName'] if kind == 'chemical' else row['GeneSymbol'])
                    node(db, disease, 'disease', row['DiseaseName'])
                    for evidence in row['DirectEvidence'].split('|'):
                        if evidence not in {'marker/mechanism', 'therapeutic'}:
                            raise ValueError('unexpected curated disease evidence')
                        edge(db, subject, disease, kind + '_disease_' + evidence.replace('/', '_'), 'curated', None,
                             {'source': 'CTD', 'direct_evidence': evidence, 'species': 'not_encoded_in_disease_export'}, record['id'], row['PubMedIDs'])
            record = files[required[2]]
            # Both ends must be disease-relevant; this is candidate graph selection, not a binding label.
            chem_tokens = [x.removeprefix('MESH:') for x in sorted(chemicals)]
            gene_tokens = [x.removeprefix('NCBIGene:') for x in sorted(genes)]
            if chem_tokens and gene_tokens:
                where = 'WHERE OrganismID=? AND ChemicalID IN (' + ','.join('?' for _ in chem_tokens) + ') AND GeneID IN (' + ','.join('?' for _ in gene_tokens) + ')'
                for row in parquet_rows(input_dir / record['parquet'], where, ['9606', *chem_tokens, *gene_tokens]):
                    subject, target = mesh_id(row['ChemicalID']), gene_id(row['GeneID'])
                    edge(db, subject, target, 'chemical_gene_influence', 'curated', '9606',
                         {'source': 'CTD', 'interaction': row['Interaction'], 'actions': sorted(set(row['InteractionActions'].split('|'))), 'gene_forms': row['GeneForms'], 'binding_label': False,
                          'human_assay_does_not_imply_clinical_trial': True}, record['id'], row['PubMedIDs'])
            for filename, kind, id_key in (('CTD_chemicals.csv.gz', 'chemical', 'ChemicalID'), ('CTD_genes.csv.gz', 'gene', 'GeneID')):
                if filename not in files:
                    continue
                record = files[filename]
                for row in parquet_rows(input_dir / record['parquet']):
                    # CTD dictionaries also contain hierarchy roots such as MESH:D.
                    # They remain in the bulk snapshot but are not molecular entities.
                    raw_id=row[id_key]
                    candidate=(raw_id if raw_id.startswith('MESH:') else 'MESH:'+raw_id) if kind=='chemical' else 'NCBIGene:'+raw_id
                    if candidate not in (chemicals if kind == 'chemical' else genes):
                        continue
                    identifier = mesh_id(raw_id) if kind == 'chemical' else gene_id(raw_id)
                    name = row['ChemicalName'] if kind == 'chemical' else row['GeneSymbol']
                    node(db, identifier, kind, name, {'dictionary': row, 'structure_verified': False})
                    keys = [('PubChemCID', 'PubChem'), ('InChIKey', 'InChIKey'), ('DTXSID', 'DTXSID')] if kind == 'chemical' else [('UniProtIDs', 'UniProt')]
                    for key, namespace in keys:
                        for value in filter(None, row.get(key, '').split('|')):
                            db.execute('INSERT OR IGNORE INTO identifier_mapping VALUES (?,?,?,?,?,?)', (identifier, namespace, value, 'source_reported', record['id'], 'CTD source cross-reference; chemical form and sequence identity require review'))
            filename = 'CTD_genes_pathways.csv.gz'
            if filename in files:
                record = files[filename]
                for row in parquet_rows(input_dir / record['parquet']):
                    identifier = gene_id(row['GeneID'])
                    if identifier in genes:
                        node(db, row['PathwayID'], 'pathway', row['PathwayName'])
                        edge(db, identifier, row['PathwayID'], 'annotated_to_pathway', 'curated', None, {'source': 'CTD', 'clinical_effect': False}, record['id'])
            if include_inferred:
                filename = 'CTD_chemicals_diseases.csv.gz'
                if filename not in files:
                    raise ValueError('--include-inferred needs the full validated chemical-disease Parquet')
                record = files[filename]
                for row in parquet_rows(input_dir / record['parquet'], 'WHERE DiseaseID IN (?,?)', DISEASES):
                    subject, disease = mesh_id(row['ChemicalID']), mesh_id(row['DiseaseID'])
                    node(db, subject, 'chemical', row['ChemicalName']); node(db, disease, 'disease', row['DiseaseName'])
                    for evidence in filter(None, row['DirectEvidence'].split('|')):
                        if evidence not in {'marker/mechanism', 'therapeutic'}:
                            raise ValueError('invalid DirectEvidence')
                        edge(db, subject, disease, 'chemical_disease_' + evidence.replace('/', '_'), 'curated', None,
                             {'source': 'CTD', 'direct_evidence': evidence, 'species': 'not_encoded_in_disease_export'}, record['id'], row['PubMedIDs'])
                    if row['InferenceGeneSymbol']:
                        score = float(row['InferenceScore']) if row['InferenceScore'] else None
                        edge(db, subject, disease, 'chemical_disease_inference', 'inferred', None, {'source': 'CTD', 'inference_genes': row['InferenceGeneSymbol'].split('|'), 'inference_score': score, 'score_is_not_probability': True}, record['id'], row['PubMedIDs'])
        result = validate(db)
        if result['errors']:
            raise ValueError(result['errors'])
        result.update({'created_at': now(), 'disease_scope': list(DISEASES), 'human_interaction_taxon': '9606', 'direct_binding_labels_from_CTD': 0, 'food_to_health_causal_paths': 0})
        db.close()
        os.replace(temporary, database)
        atomic_json(database.with_suffix('.coverage.json'), result)
        return result
    except BaseException:
        db.close(); temporary.unlink(missing_ok=True)
        raise


def structure_coverage(db) -> dict:
    """Own molecular structure availability differs from cross-database identity."""
    rows=db.execute("SELECT metadata_json FROM node WHERE kind='chemical'").fetchall()
    verified=sum(bool(json.loads(r[0]).get('full_inchikey')) for r in rows)
    return {'chemical_nodes':len(rows),'chemical_nodes_with_full_structure':verified,
            'chemical_nodes_without_full_structure':len(rows)-verified,
            'cross_database_exact_structure_mappings':db.execute("SELECT count(*) FROM identifier_mapping WHERE relation='exact_structure'").fetchone()[0]}

def validate(db) -> dict:
    errors = []
    integrity = db.execute('PRAGMA integrity_check').fetchone()[0]
    if integrity != 'ok': errors.append(integrity)
    if db.execute('PRAGMA foreign_key_check').fetchall(): errors.append('foreign key violations')
    if db.execute("SELECT count(*) FROM edge WHERE relation='chemical_gene_influence' AND (taxon!='9606' OR tier!='curated')").fetchone()[0]: errors.append('invalid human influence tier/taxon')
    if db.execute("SELECT count(*) FROM edge WHERE relation LIKE '%inference' AND tier!='inferred'").fetchone()[0]: errors.append('inferred promoted to curated')
    if db.execute("SELECT count(*) FROM measurement WHERE (status IN ('missing','below_detection','unparsed','invalid_source_value') AND amount IS NOT NULL) OR (status='measured_zero' AND amount!=0) OR amount<0").fetchone()[0]: errors.append('invalid measured/missing values')
    return {'integrity': integrity, 'errors': errors, 'tables': {table: db.execute(f'SELECT count(*) FROM {table}').fetchone()[0] for table in ('source','node','edge','citation','edge_source','identifier_mapping','measurement')},
            'relations': {f'{r[0]}:{r[1]}':r[2] for r in db.execute('SELECT relation,tier,count(*) FROM edge GROUP BY relation,tier')},
            'unresolved_structure_nodes': db.execute("SELECT count(*) FROM node WHERE kind='chemical' AND id NOT IN (SELECT subject FROM identifier_mapping WHERE relation='exact_structure')").fetchone()[0],
            'unresolved_structure_nodes_definition':'Legacy count lacking an exact_structure mapping; does not mean all nodes lack their own structure',
            'structure_coverage':structure_coverage(db)}


def query(db, subject: str | None = None, disease: str | None = None, tier: str | None = None, limit: int = 20, offset: int = 0) -> list:
    if not 1 <= limit <= 100000 or offset < 0: raise ValueError('invalid pagination')
    clauses, args = [], []
    for column, value in (('subject',subject),('object',disease),('tier',tier)):
        if value is not None: clauses.append(column+'=?'); args.append(value)
    sql = 'SELECT * FROM edge' + (' WHERE '+' AND '.join(clauses) if clauses else '') + ' ORDER BY id LIMIT ? OFFSET ?'
    result=[]
    for row in db.execute(sql,[*args,limit,offset]):
        item=dict(row);item['context']=json.loads(item.pop('context_json'))
        item['pmids']=[r[0] for r in db.execute('SELECT pmid FROM citation WHERE edge_id=? ORDER BY pmid',(row['id'],))]
        item['source_ids']=[r[0] for r in db.execute('SELECT source_id FROM edge_source WHERE edge_id=? ORDER BY source_id',(row['id'],))]
        result.append(item)
    return result
