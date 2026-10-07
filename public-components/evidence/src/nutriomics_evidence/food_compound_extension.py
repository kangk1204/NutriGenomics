"""Strict native food ID/structure/occurrence bridge for an isolated source."""
from __future__ import annotations
import bz2
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import duckdb
import pandas as pd

FULLKEY=re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
AGGREGATE=re.compile(r'^(total\b|protein\b|carbohydrates\b|dietary fiber\b)|\bmixture\b',re.I)

def verify_structure(row):
    from rdkit import Chem
    from rdkit.Chem import inchi
    source_id,public_id,name,raw_smiles,source_inchi,full_key=row
    valid_key=bool(full_key and FULLKEY.fullmatch(full_key))
    computed_inchi=inchi.InchiToInchiKey(source_inchi) if source_inchi and source_inchi.startswith('InChI=') else ''
    molecule=Chem.MolFromSmiles(raw_smiles) if raw_smiles else None
    computed_smiles=Chem.MolToInchiKey(molecule) if molecule is not None else ''
    canonical=Chem.MolToSmiles(molecule,isomericSmiles=True) if molecule is not None else ''
    single=molecule is not None and len(Chem.GetMolFrags(molecule))==1
    stereo=[] if molecule is None else list(Chem.FindPotentialStereo(molecule))
    unresolved=any(item.specified==Chem.StereoSpecified.Unspecified for item in stereo)
    aggregate=bool(name and AGGREGATE.search(name))
    verified=valid_key and computed_inchi==full_key and computed_smiles==full_key and single and not aggregate
    return {'source_compound_id':source_id,'public_compound_id':public_id,'name':name,'source_smiles':raw_smiles,'source_inchi':source_inchi,'full_inchikey':full_key,
        'canonical_isomeric_smiles':canonical,'computed_inchi_key':computed_inchi,'computed_smiles_key':computed_smiles,'full_key_format_valid':valid_key,
        'single_component':single,'aggregate_concept':aggregate,'potential_stereo_elements':len(stereo),'standard_unspecified_stereo':unresolved,
        'standard_structure_verified':verified,'native_analytical_stereo_or_isomer_resolved':False,'source_id':'foodb_sep2022_author_snapshot'}

def curated_parts(path):
    text=path.read_text()
    ids=re.search(r'food_id\s*=\s*c\((.*?)\),\s*desired_part',text,re.S)
    parts=re.search(r'desired_part\s*=\s*c\((.*?)\)\s*\)',text,re.S)
    if not ids or not parts:raise ValueError('published edible-part mapping schema missing')
    identifiers=re.findall(r'\b\d+\b',ids[1]);names=re.findall(r'[\"\']([^\"\']+)[\"\']',parts[1])
    if len(identifiers)!=len(names):raise ValueError('unpaired published edible parts')
    return pd.DataFrame(sorted(set(zip(identifiers,names))),columns=['food_id','orig_food_part'])

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def literal(value):return "'"+str(value).replace("'","''")+"'"

def verified_inputs(root):
    m=json.loads((root/'acquisition.json').read_text())
    for item in m['sources']:
        if digest(root/'raw'/item['name'])!=item['sha256']:raise ValueError('source hash mismatch:'+item['name'])
    return m

def relation(path):return f"read_csv({literal(path)},header=true,all_varchar=true,strict_mode=true,max_line_size=33554432)"

def extract(root):
    verified_inputs(root);destination=root/'extracted';destination.mkdir(exist_ok=True)
    for name in ['Compound.csv','Content.csv']:
        target=destination/name
        if not target.exists():
            import shutil
            with bz2.open(root/'raw'/(name+'.bz2'),'rb') as source,target.open('wb') as f:shutil.copyfileobj(source,f,8*1024*1024)
    return destination

def inspect(root):
    extracted=extract(root);con=duckdb.connect();con.execute('SET threads=2');con.execute('SET enable_progress_bar=false')
    info={}
    for name,path in [('compound',extracted/'Compound.csv'),('content',extracted/'Content.csv'),('food',root/'raw/Food.csv')]:
        con.execute(f'CREATE TABLE {name} AS SELECT * FROM {relation(path)}')
        columns=[row[0] for row in con.execute(f'DESCRIBE {name}').fetchall()]
        info[name]={'rows':con.execute(f'SELECT count(*) FROM {name}').fetchone()[0],'columns':columns,'sample':con.execute(f'SELECT * FROM {name} LIMIT 2').fetchdf().to_dict(orient='records')}
    for table,column in [('content','source_type'),('content','citation'),('content','orig_unit'),('content','orig_food_part'),('food','food_group'),('food','food_type')]:
        if column in info[table]['columns']:info.setdefault('value_counts',{})[table+'.'+column]=con.execute(f'SELECT "{column}",count(*) FROM {table} GROUP BY "{column}" ORDER BY count(*) DESC LIMIT 35').fetchall()
    (root/'source_inspection.json').write_text(json.dumps(info,indent=2),encoding='utf-8');con.close()
    print(json.dumps({k:v for k,v in info.items() if k=='value_counts'}),flush=True)
    return info

def stage(root):
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.warning');RDLogger.DisableLog('rdApp.error')
    manifest=verified_inputs(root);extracted=extract(root);derived=root/'derived';derived.mkdir(exist_ok=True)
    database=derived/'food_compound_extension.duckdb';con=duckdb.connect(str(database));con.execute('SET threads=2');con.execute("SET memory_limit='12GB'");con.execute('SET enable_progress_bar=false')
    for name,path in [('compound',extracted/'Compound.csv'),('content',extracted/'Content.csv'),('food',root/'raw/Food.csv')]:
        con.execute(f'CREATE OR REPLACE TABLE native_{name} AS SELECT * FROM {relation(path)}')
    # The study's published data dictionary explicitly documents shifted labels.
    # No name matching or guessed CAS lookup substitutes for exact structures.
    structure_path=derived/'food_compound_extension_structure_records.parquet'
    source_signature=derived/'structure_input_signature.json'
    import inspect as code_inspect
    signature={x['name']:x['sha256'] for x in manifest['sources'] if x['name']=='Compound.csv.bz2'}
    signature['verifier_sha256']=hashlib.sha256(code_inspect.getsource(verify_structure).encode()).hexdigest()
    cached_signature=json.loads(source_signature.read_text()) if source_signature.exists() else {}
    if structure_path.exists() and cached_signature.get('inputs')==signature and cached_signature.get('output_sha256')==digest(structure_path):
        structure_frame=pd.read_parquet(structure_path)
    else:
        frame=con.execute('SELECT id,public_id,name,cas_number,moldb_inchikey,moldb_smiles FROM native_compound').fetchdf()
        normalized=[verify_structure(tuple(None if pd.isna(x) else x for x in row)) for row in frame.itertuples(index=False,name=None)]
        structure_frame=pd.DataFrame(normalized);structure_frame.to_parquet(structure_path,index=False)
        source_signature.write_text(json.dumps({'inputs':signature,'output_sha256':digest(structure_path)},sort_keys=True),encoding='utf-8')
    print(json.dumps({'stage':'exact_structure_verification','native':len(structure_frame),'valid_standard_rows':int(structure_frame.standard_structure_verified.sum())}),flush=True)
    def view(name,path):con.execute(f'CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet({literal(path)})')
    def materialize(name,query):
        path=derived/(name+'.parquet');con.execute(f'COPY ({query}) TO {literal(path)} (FORMAT PARQUET,COMPRESSION ZSTD)');view(name,path)
    view('food_compound_extension_structure_records',derived/'food_compound_extension_structure_records.parquet')
    parts=curated_parts(root/'raw/edible_part_rules.Rmd');parts.to_parquet(derived/'food_compound_extension_curated_parts.parquet',index=False)
    view('food_compound_extension_curated_parts',derived/'food_compound_extension_curated_parts.parquet')
    materialize('food_compound_extension_native_food','SELECT * FROM native_food')
    materialize('food_compound_extension_native_content','SELECT * FROM native_content')
    positive="greatest("+','.join(f"CASE WHEN isfinite(try_cast({field} AS DOUBLE)) THEN try_cast({field} AS DOUBLE) END" for field in ['orig_content','orig_min','orig_max','standard_content'])+")"
    con.execute(f"""CREATE OR REPLACE TEMP VIEW occurrence_candidates AS SELECT c.*,s.* EXCLUDE(source_compound_id,source_id),
      f.public_id public_food_id,f.name food_name,f.food_group,f.food_type,f.export_to_foodb,
      CASE WHEN p.food_id IS NOT NULL THEN true ELSE false END published_edible_part_verified,
      {positive} positive_value,
      CASE WHEN orig_unit IN ('mg/100g','mg/100 g','mg/100 g fresh weight','mg/100 g freshweight') AND {positive}>100000 THEN true ELSE false END impossible_mass,
      CASE WHEN regexp_full_match(coalesce(c.citation,''),'[0-9]{{6,10}}') THEN c.citation ELSE NULL END source_pmid,
      'foodb_sep2022_author_snapshot' source_dataset_id
      FROM native_content c LEFT JOIN food_compound_extension_structure_records s ON (CASE WHEN c.source_type='Compound' THEN c.source_id END)=s.source_compound_id
      LEFT JOIN native_food f ON c.food_id=f.id LEFT JOIN food_compound_extension_curated_parts p USING(food_id,orig_food_part)""")
    con.execute("""CREATE OR REPLACE TEMP VIEW classified_occurrences AS SELECT *,
      CASE WHEN source_type!='Compound' THEN 'nutrient_namespace_not_single_compound'
        WHEN upper(coalesce(citation_type,''))='PREDICTED' OR upper(coalesce(citation,'')) IN ('HMDB','PATHBANK') OR regexp_matches(lower(coalesce(orig_method,'')),'predict|infer|in.silico') THEN 'predicted_not_observed_food'
        WHEN coalesce(standard_structure_verified,false)=false THEN 'unverified_or_mixture_or_aggregate_structure'
        WHEN public_food_id IS NULL OR food_name IS NULL OR food_group IS NULL OR food_type='Unknown' THEN 'unresolved_food_record'
        WHEN impossible_mass THEN 'impossible_mass_per100g'
        WHEN citation='DUKE' AND NOT published_edible_part_verified THEN 'duke_part_not_in_published_food_specific_curation'
        WHEN citation='DUKE' AND published_edible_part_verified AND positive_value>0 THEN 'eligible_positive_curated_edible_part'
        WHEN citation='DUKE' AND published_edible_part_verified AND positive_value IS NULL THEN 'eligible_qualitative_curated_edible_part'
        WHEN citation IN ('USDA','DTU','PHENOL EXPLORER') AND positive_value>0 AND orig_food_common_name IS NOT NULL THEN 'eligible_positive_food_composition'
        WHEN citation_type IN ('ARTICLE','TEXTBOOK') AND citation NOT IN ('MANUAL','UNKNOWN') AND length(trim(coalesce(citation,'')))>0 AND positive_value>0 THEN 'eligible_positive_published_food_record'
        WHEN citation_type IN ('ARTICLE','TEXTBOOK') AND citation NOT IN ('MANUAL','UNKNOWN') AND length(trim(coalesce(citation,'')))>0 AND positive_value IS NULL
          AND orig_food_common_name IS NOT NULL AND NOT regexp_matches(lower(citation),'leaves|\\bleaf\\b|\\bbark\\b|tissue.culture') THEN 'eligible_qualitative_published_food_record'
        WHEN citation='MANUAL' AND citation_type='EXPERIMENTAL' AND positive_value>0 AND length(trim(coalesce(orig_citation,'')))>0
          AND NOT regexp_matches(lower(orig_citation),'in preparation|unpublished|submitted') THEN 'eligible_positive_referenced_experiment'
        WHEN positive_value<=0 THEN 'no_positive_occurrence_value'
        ELSE 'missing_specific_observed_food_evidence' END eligibility
      FROM occurrence_candidates""")
    materialize('food_compound_extension_eligible_occurrences',"SELECT * FROM classified_occurrences WHERE starts_with(eligibility,'eligible_')")
    materialize('food_compound_extension_rejection_ledger','SELECT id source_content_id,source_type,source_id source_compound_id,food_id,citation,citation_type,eligibility FROM classified_occurrences WHERE NOT starts_with(eligibility,\'eligible_\')')
    materialize('food_compound_extension_standards',"""SELECT * EXCLUDE(rn) FROM
       (SELECT s.*,row_number() OVER(PARTITION BY s.full_inchikey ORDER BY s.source_compound_id) rn
        FROM food_compound_extension_structure_records s WHERE standard_structure_verified AND full_inchikey IN
        (SELECT DISTINCT full_inchikey FROM food_compound_extension_eligible_occurrences)) WHERE rn=1""")
    materialize('food_compound_extension_source_ledger','SELECT DISTINCT citation,citation_type,orig_citation,orig_method,orig_unit,preparation_type,source_pmid FROM food_compound_extension_eligible_occurrences')
    queries={
      'duplicate_full_structure':'SELECT count(*)-count(DISTINCT full_inchikey) FROM food_compound_extension_standards',
      'bad_full_key_or_smiles_or_mixture':'SELECT count(*) FROM food_compound_extension_standards WHERE NOT standard_structure_verified OR NOT single_component OR computed_inchi_key!=full_inchikey OR computed_smiles_key!=full_inchikey',
      'predicted_occurrences_counted':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE citation_type='PREDICTED' OR citation IN ('PATHBANK','HMDB')",
      'unpaired_duke_edible_part':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE citation='DUKE' AND NOT published_edible_part_verified",
      'impossible_mass_counted':'SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE impossible_mass',
      'eligible_orphan_food':'SELECT count(*) FROM food_compound_extension_eligible_occurrences o ANTI JOIN food_compound_extension_native_food f ON o.food_id=f.id',
      'eligible_orphan_standard':'SELECT count(*) FROM food_compound_extension_eligible_occurrences o ANTI JOIN food_compound_extension_standards s USING(full_inchikey)',
      'native_exposure_falsely_resolved':'SELECT count(*) FROM food_compound_extension_standards WHERE native_analytical_stereo_or_isomer_resolved'}
    checks={name:con.execute(q).fetchone()[0] for name,q in queries.items()}
    native={name:con.execute(f'SELECT count(*) FROM native_{name}').fetchone()[0] for name in ['compound','content','food']}
    counts={name:con.execute('SELECT count(*) FROM '+name).fetchone()[0] for name in ['food_compound_extension_eligible_occurrences','food_compound_extension_standards','food_compound_extension_rejection_ledger','food_compound_extension_source_ledger']}
    compound_count=counts['food_compound_extension_standards']
    source_counts=con.execute('SELECT citation,citation_type,count(*),count(DISTINCT full_inchikey) FROM food_compound_extension_eligible_occurrences GROUP BY ALL ORDER BY count(*) DESC').fetchall()
    tier_counts=con.execute('SELECT eligibility,count(*),count(DISTINCT full_inchikey) FROM food_compound_extension_eligible_occurrences GROUP BY eligibility ORDER BY eligibility').fetchall()
    rejects=dict(con.execute('SELECT eligibility,count(*) FROM food_compound_extension_rejection_ledger GROUP BY eligibility').fetchall())
    positive_standards=con.execute('SELECT count(DISTINCT full_inchikey) FROM food_compound_extension_eligible_occurrences WHERE positive_value>0').fetchone()[0]
    unspecified=con.execute('SELECT count(*) FROM food_compound_extension_standards WHERE standard_unspecified_stereo').fetchone()[0]
    for name in ['native_compound','native_content','native_food']:con.execute('DROP TABLE '+name)
    con.execute('CHECKPOINT');con.close()
    receipt={'schema_version':1,'frozen_utc':datetime.now(timezone.utc).isoformat(),'source_manifest_sha256':digest(root/'acquisition.json'),'native_source_counts':native,
       'eligible_counts':counts,'verified_food_single_compound_standard_count':compound_count,'positive_numeric_standard_count':positive_standards,
       'qualitative_only_standard_count':compound_count-positive_standards,'standard_unspecified_stereo_count':unspecified,
       'meets_2000_occurrence_verified_standard_target':compound_count>=2000,'individual_native_analytical_forms_resolved':False,
       'quantity_measured_for_every_standard':False,'independently_validated_new_food_occurrence_experiments':False,
       'occurrence_definition':'Exact native food ID and Compound namespace ID joined to source standard structure; observed/nonpredicted food composition or specific published record; DUKE only exact foodID/part pairs from published dietary study curation; qualitative database/primary paper occurrence distinct from positive numeric records',
       'structure_schema_binding':{'source_SMILES':'cas_number','source_InChI':'moldb_inchikey','source_fullInChIKey':'moldb_smiles','authority':'Published author FooDB dictionary documents shifted labels; all keys recomputed independently from InChI and SMILES; single connected molecule; stereo retained'},
       'native_form_limits':['Keys verify standardized chemical representations, not resolved native food stereoisomers','Milk PMID30994344 groups296 measured metabolites/species corresponding to1447 possible unique structures; source values can be molecular-species assignments','Original units/preparation/food parts/source citations retained; no dose or common participant linkage inferred'],
       'source_group_counts':source_counts,'evidence_tier_counts':tier_counts,'rejection_counts':rejects,'checks':checks,'audit_queries':queries,'all_checks_zero':all(value==0 for value in checks.values()),
       'frozen_previous_database_modified':False,'raw_git_commit_allowed':False,'code_sha256':{path.name:digest(path) for path in sorted((root/'code').glob('*.py'))},
       'outputs':[{'path':str(path),'sha256':digest(path),'bytes':path.stat().st_size} for path in sorted(derived.glob('*')) if path.suffix in ['.parquet','.duckdb']]}
    (derived/'receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['outputs','code_sha256','audit_queries','native_form_limits','structure_schema_binding']}),flush=True)
    if not receipt['all_checks_zero']:raise ValueError('extension structural audit failed')
    return receipt


def finalize(root):
    """Audit every staged row read-only, then write a compact separate receipt."""
    manifest=verified_inputs(root);derived=root/'derived';stage_receipt=json.loads((derived/'receipt.json').read_text())
    database=derived/'food_compound_extension.duckdb';con=duckdb.connect(str(database),read_only=True);con.execute('SET threads=2')
    q={
      'duplicate_native_content_id':'SELECT count(*)-count(DISTINCT id) FROM food_compound_extension_native_content',
      'duplicate_native_food_id':'SELECT count(*)-count(DISTINCT id) FROM food_compound_extension_native_food',
      'duplicate_eligible_content_id':'SELECT count(*)-count(DISTINCT id) FROM food_compound_extension_eligible_occurrences',
      'source_row_partition_difference':'SELECT (SELECT count(*) FROM food_compound_extension_native_content)-(SELECT count(*) FROM food_compound_extension_eligible_occurrences)-(SELECT count(*) FROM food_compound_extension_rejection_ledger)',
      'noncompound_namespace_counted':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE source_type IS DISTINCT FROM 'Compound'",
      'unresolved_food_counted':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE public_food_id IS NULL OR food_name IS NULL OR food_group IS NULL OR food_type='Unknown'",
      'positive_tier_without_positive_value':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE starts_with(eligibility,'eligible_positive_') AND (positive_value IS NULL OR NOT isfinite(positive_value) OR positive_value<=0)",
      'qualitative_tier_with_numeric_value':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE starts_with(eligibility,'eligible_qualitative_') AND positive_value IS NOT NULL",
      'taxonomy_or_uncited_groups_counted':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE citation IN ('KNAPSACK','PHYTOHUB','DFC','HMDB','PATHBANK','UNKNOWN')",
      'noneligible_tier_counted':"SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE NOT starts_with(eligibility,'eligible_')"
    }
    q.update(stage_receipt['audit_queries'])
    checks={name:con.execute(sql).fetchone()[0] for name,sql in q.items()}
    count_sql='SELECT count(DISTINCT full_inchikey) FROM food_compound_extension_standards WHERE standard_structure_verified AND single_component AND computed_inchi_key=full_inchikey AND computed_smiles_key=full_inchikey'
    unique=con.execute(count_sql).fetchone()[0]
    measured=con.execute("SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE starts_with(eligibility,'eligible_positive_')").fetchone()[0]
    qualitative=con.execute("SELECT count(*) FROM food_compound_extension_eligible_occurrences WHERE starts_with(eligibility,'eligible_qualitative_')").fetchone()[0]
    foods=con.execute('SELECT count(DISTINCT food_id) FROM food_compound_extension_eligible_occurrences').fetchone()[0]
    milk=con.execute("SELECT count(*),count(DISTINCT full_inchikey) FROM food_compound_extension_eligible_occurrences WHERE source_pmid='30994344'").fetchone()
    parts=con.execute('SELECT count(*) FROM food_compound_extension_curated_parts').fetchone()[0]
    source_units=con.execute('SELECT orig_unit,count(*) FROM food_compound_extension_eligible_occurrences GROUP BY orig_unit ORDER BY count(*) DESC').fetchall()
    nonmilk=con.execute("SELECT count(DISTINCT full_inchikey) FROM food_compound_extension_eligible_occurrences WHERE source_pmid IS DISTINCT FROM '30994344'").fetchone()[0]
    con.close()
    previous=root.parent/'derived/expansion.duckdb';expected='91ac8845e1addd2b94191af48116ecd1cf6762c7343b6070a85229bd6207d033'
    previous_hash=digest(previous) if previous.exists() else None
    audit={'all_row_checks':checks,'audit_queries':q,'all_row_checks_passed':all(x==0 for x in checks.values()),'check_count':len(checks),'read_only':True}
    (derived/'allrow_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    outputs=[{'path':str(p),'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(derived.glob('*')) if p.suffix in ['.parquet','.duckdb']]
    receipt={'schema_version':1,'protocol_id':'food_compound_extension_20261003_v1','frozen_utc':datetime.now(timezone.utc).isoformat(),
       'metrics':{'eligible_unique_food_compounds':unique,'measured_food_occurrence_rows':measured,'qualitative_food_occurrence_rows':qualitative,'eligible_food_occurrence_rows':measured+qualitative,'valid_foods':foods,
          'positive_numeric_unique_food_compounds':stage_receipt['positive_numeric_standard_count'],'qualitative_only_unique_food_compounds':stage_receipt['qualitative_only_standard_count'],
          'native_source_compounds':stage_receipt['native_source_counts']['compound'],'native_source_content_rows':stage_receipt['native_source_counts']['content'],'native_source_foods':stage_receipt['native_source_counts']['food'],
          'rejected_source_content_rows':stage_receipt['eligible_counts']['food_compound_extension_rejection_ledger'],'published_food_specific_part_pairs':parts,'standard_unspecified_stereo':stage_receipt['standard_unspecified_stereo_count'],
          'milk_species_assigned_occurrence_rows':milk[0],'milk_species_assigned_unique_standard_keys':milk[1],'unique_keys_with_nonmilk_evidence':nonmilk},
       'count_unit':'Distinct exact full InChIKey standardized single chemical structures with specific nonpredicted food occurrence; native analytical stereoisomer resolution is not claimed',
       'count_sql':count_sql,'meets_2000_verified_food_compound_target':unique>=2000 and audit['all_row_checks_passed'],
       'physical_database_path':str(database),'physical_database_sha256':digest(database),
       'audit_passed':audit['all_row_checks_passed'],'all_row_checks':checks,'all_row_audit_path':str(derived/'allrow_audit.json'),'all_row_audit_sha256':digest(derived/'allrow_audit.json'),
       'source_manifest_path':str(root/'acquisition.json'),'source_manifest_sha256':digest(root/'acquisition.json'),'stage_receipt_path':str(derived/'receipt.json'),'stage_receipt_sha256':digest(derived/'receipt.json'),
       'sources':manifest['sources'],'source_repository':'https://github.com/SWi1/FooDB_polyphenol_analysis','source_repository_commit':'3d2bcf9523911fabc8a07c731226ac87a2af262e','source_version':'September2022 historical author-provided FooDB snapshot; not claimed as current2026 FooDB',
       'source_paper_doi':'10.1016/j.tjnut.2024.08.010','source_author_verification_url':'https://www.ars.usda.gov/research/publications/publication/?seqNo115=410641',
       'licenses':{'data':'CC-BY-NC-4.0','data_url':'https://foodb.ca/compliance','terms_url':'https://creativecommons.org/licenses/by-nc/4.0/','code_in_author_repository':'MIT',
          'code_license_file_sha256':next(x['sha256'] for x in manifest['sources'] if x['name']=='source_code_LICENSE'),'project_raw_git_commit_allowed':False,'noncommercial_research_reuse_allowed':True,'commercial_permission_required':True},
       'eligibility_definition':stage_receipt['occurrence_definition'],'source_group_counts':stage_receipt['source_group_counts'],'evidence_tier_counts':stage_receipt['evidence_tier_counts'],'rejections':stage_receipt['rejection_counts'],'native_units':source_units,
       'mapping_limitations':stage_receipt['native_form_limits']+['Full-key agreement verifies the source standard representation; it does not prove every native food exposure is stereochemically or positionally resolved','Positive numeric source rows can refer to molecular-species assignments and censored ranges; this receipt counts source composition records rather than new independent measurements','Qualitative occurrence is retained only with specific published food records or exact author-curated edible-part pairs; amounts are not inferred','No automatic extension of the immutable33-compound DTI evaluation panel; this coverage database has not been evaluated as an independent DTI food panel'],
       'previous_frozen_database_path':str(previous),'previous_frozen_database_sha256':previous_hash,'previous_frozen_database_expected_sha256':expected,'previous_frozen_database_unchanged':previous_hash==expected,
       'outputs':outputs,'finalization_code_sha256':{p.name:digest(p) for p in sorted((root/'final_code' if (root/'final_code').exists() else root/'code').glob('*.py'))}}
    (derived/'extension_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps({'metrics':receipt['metrics'],'audit_passed':receipt['audit_passed'],'database_sha256':receipt['physical_database_sha256'],'previous_frozen_database_unchanged':receipt['previous_frozen_database_unchanged']}),flush=True)
    if not receipt['audit_passed'] or not receipt['previous_frozen_database_unchanged']:raise ValueError('final extension audit failed')
    return receipt


GENERIC_ANALYTE=re.compile(r'\b(total|equivalents?|mixture|polyphenols|flavonoids|tocopherols|folates|carbohydrates|proteins)\b|fatty acids|vitamin\s*[ade]\b|\b(?:pc|pe|tg|dg|sm|lpc|lpe|cer)\s*[\[(]|phosphatidyl|lysophosphatidyl|sphingomyelin|ceramide|triacylglycer|diacylglycer',re.I)

def normalized_analyte_name(value):
    import unicodedata
    # Retain +/-, positional and stereochemical labels. No fuzzy synonym matching.
    return ' '.join(unicodedata.normalize('NFKC',str(value or '')).casefold().split())

def strict_analyte_status(row,single_atom=False):
    quantitative=row.get('positive_value') is not None and row['positive_value']>0
    qualitative=str(row.get('eligibility') or '').startswith('eligible_qualitative_')
    if not quantitative and not qualitative:return 'neither_positive_nor_specific_qualitative_occurrence'
    if not row.get('standard_structure_verified') or not row.get('single_component'):return 'structure_unverified'
    if str(row.get('source_pmid') or '')=='30994344':return 'milk_assay_to_structure_mapping_unresolved'
    if GENERIC_ANALYTE.search(str(row.get('name') or '')) or GENERIC_ANALYTE.search(str(row.get('orig_source_name') or '')):return 'generic_class_or_lipid_species_assignment'
    if single_atom:return 'element_total_not_resolved_molecular_form'
    if row.get('standard_unspecified_stereo'):return 'standard_stereochemistry_unspecified'
    if not row.get('orig_source_id') or not row.get('orig_source_name'):return 'original_assayed_analyte_id_or_name_missing'
    if normalized_analyte_name(row['orig_source_name'])!=normalized_analyte_name(row['name']):return 'assayed_analyte_name_not_exact_standard_identity'
    if quantitative and not row.get('orig_unit'):return 'native_measurement_unit_missing'
    if re.search(r'hydroly|saponif|equivalent|predict|infer|in.silico',str(row.get('orig_method') or ''),re.I):return 'native_form_or_method_ambiguous'
    if row.get('citation')=='DUKE' and not row.get('published_edible_part_verified'):return 'edible_part_unverified'
    if qualitative:
        specific=row.get('orig_citation') if row.get('citation') in ('DUKE','USDA','DTU','PHENOL EXPLORER') else row.get('citation')
        text=str(specific or '').strip()
        if not text or text.upper() in ('MANUAL','UNKNOWN','DUKE','USDA','DTU','PHENOL EXPLORER','PHYTOHUB','DFC','KNAPSACK'):return 'qualitative_database_assignment_without_specific_source'
        if not re.search(r'10\.\d{4,9}/|\b(?:19|20)\d{2}\b|^\d{6,10}$',text):return 'qualitative_source_not_specific_publication'
        if re.search(r'corn silk|asparagus root|tissue.culture|insect.attractant',text,re.I):return 'reported_nonfood_tissue_or_purpose'
        return 'strict_source_reported_qualitative_single_analyte'
    return 'strict_positive_source_named_single_analyte'

def strict_audit(root,version=1):
    """Separate conservative positive-only audit; leave broader DB immutable."""
    from rdkit import Chem
    derived=root/'derived';destination=derived/('strict' if version==1 else 'strict_v2');destination.mkdir(exist_ok=True)
    broad=json.loads((derived/'extension_receipt.json').read_text());database=Path(broad['physical_database_path'])
    before=digest(database);con=duckdb.connect(str(database),read_only=True);con.execute('SET threads=2')
    fields=['id','food_id','public_food_id','source_id','public_compound_id','full_inchikey','name','source_inchi','canonical_isomeric_smiles','positive_value','source_pmid','citation','citation_type','orig_source_id','orig_source_name','orig_unit','orig_method','orig_citation','eligibility','orig_food_common_name','orig_food_part','published_edible_part_verified','standard_structure_verified','single_component','standard_unspecified_stereo']
    frame=con.execute('SELECT '+','.join(fields)+' FROM food_compound_extension_eligible_occurrences').fetchdf()
    broad_keys=con.execute("SELECT count(DISTINCT full_inchikey),count(DISTINCT full_inchikey) FILTER(WHERE source_pmid='30994344'),count(DISTINCT full_inchikey) FILTER(WHERE source_pmid IS DISTINCT FROM '30994344') FROM food_compound_extension_eligible_occurrences WHERE positive_value>0").fetchone()
    con.close()
    singles={}
    for key,smiles in frame[['full_inchikey','canonical_isomeric_smiles']].drop_duplicates().itertuples(index=False,name=None):
        molecule=Chem.MolFromSmiles(smiles);singles[key]=molecule is not None and molecule.GetNumAtoms()==1
    statuses=[]
    for native in frame.itertuples(index=False,name=None):
        row={name:None if pd.isna(value) else value for name,value in zip(fields,native)}
        statuses.append(strict_analyte_status(row,singles[row['full_inchikey']]))
    frame['strict_status']=statuses
    frame.to_parquet(destination/'strict_occurrence_audit.parquet',index=False)
    accepted=frame[frame.strict_status.str.startswith('strict_')]
    quantitative=accepted[accepted.strict_status=='strict_positive_source_named_single_analyte']
    qualitative=accepted[accepted.strict_status=='strict_source_reported_qualitative_single_analyte']
    standards=accepted.sort_values(['full_inchikey','id']).drop_duplicates('full_inchikey')
    standards.to_parquet(destination/'strict_verified_food_compounds.parquet',index=False)
    counts=frame.groupby('strict_status').agg(rows=('id','size'),unique_full_keys=('full_inchikey','nunique')).reset_index().to_dict('records')
    checks={'duplicate_strict_full_keys':len(standards)-standards.full_inchikey.nunique(),
        'milk_assignments_counted':int((accepted.source_pmid=='30994344').sum()),'unspecified_standard_stereo_counted':int(accepted.standard_unspecified_stereo.sum()),
        'missing_source_analyte_identity_counted':int((accepted.orig_source_id.isna()|accepted.orig_source_name.isna()).sum()),'eligible_rows_partition_difference':len(frame)-sum(x['rows'] for x in counts),
        'broader_database_changed':int(digest(database)!=before)}
    receipt={'schema_version':1,'protocol_id':f'food_compound_extension_source_analyte_audit_v{version}','frozen_utc':datetime.now(timezone.utc).isoformat(),
       'strict_unique_food_compounds_union':len(standards),'strict_positive_unique_food_compounds':quantitative.full_inchikey.nunique(),'strict_positive_food_occurrence_rows':len(quantitative),
       'strict_qualitative_unique_food_compounds':qualitative.full_inchikey.nunique(),'strict_qualitative_food_occurrence_rows':len(qualitative),'strict_qualitative_only_unique_food_compounds':len(set(qualitative.full_inchikey)-set(quantitative.full_inchikey)),
       'strict_quantitative_qualitative_overlap_keys':len(set(qualitative.full_inchikey)&set(quantitative.full_inchikey)),'strict_valid_foods':accepted.food_id.nunique(),
       'count_unit':'Source-reported single-analyte occurrence with original analyteID+exact chemical name matching its full-key-verified standard representation. Positive numeric and amount-free identified presence with specific publication are counted separately then unioned. Unresolved milk assignments, class/species totals, unspecified stereo and method ambiguities are excluded',
       'broader_standard_concept_count':broad['metrics']['eligible_unique_food_compounds'],'positive_numeric_standard_keys':broad_keys[0],'milk_assigned_positive_keys':broad_keys[1],'positive_nonmilk_keys_before_strict_guards':broad_keys[2],
       'milk_positive_keys_with_nonmilk_evidence':broad_keys[1]+broad_keys[2]-broad_keys[0],'strict_target_2000_passed':len(standards)>=2000 and all(x==0 for x in checks.values()),
       'original_assay_mapping_for_milk_available':False,'all_uncertain_milk_species_rows_excluded':True,'checks':checks,'audit_passed':all(x==0 for x in checks.values()),'rejection_categories':counts,
       'source_database_path':str(database),'source_database_sha256':before,'broader_receipt_path':str(derived/'extension_receipt.json'),'broader_receipt_sha256':digest(derived/'extension_receipt.json'),
       'source_license':broad['licenses'],'source_links':{'repository':broad['source_repository'],'paper':broad['source_paper_doi'],'milk_primary':'https://pubmed.ncbi.nlm.nih.gov/30994344/'},
       'limitations':['This conservative exact-name guard can exclude valid synonym-linked analytes; it is not an estimate of all food chemistry','Original chemical names support source analyte identity after native CompoundID and exact structure verification; names alone never create new structures','Missing amount does not exclude presence-only records with specific source evidence; blank amount plus only a database/taxonomy assignment is excluded','Native free-ligand exposure, dose, analytical stereoisomer resolution and health efficacy remain unestablished','The broader4630 standardized representations remain separate; they do not automatically pass the strict2000 verified-single-compound target'],
       'outputs':[{'path':str(p),'sha256':digest(p),'bytes':p.stat().st_size} for p in sorted(destination.glob('*.parquet'))],
       'executed_code_sha256':{p.name:digest(p) for p in sorted((root/('strict_code' if version==1 else 'strict_v2_code')).glob('*.py'))}}
    if version==2:
        receipt['superseded_strict_v1_receipt_sha256']=digest(derived/'strict/strict_receipt.json')
        receipt['correction']='Exclude original DTU nutrient0023 generic Vitamin D total: exact source-name equality cannot resolve native D2/D3 molecular identity. Earlier strict264 snapshot remains unchanged.'
    (destination/'strict_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_license','source_links','limitations','outputs','executed_code_sha256']}),flush=True)
    if not receipt['audit_passed']:raise ValueError('strict analyte audit failed')
    return receipt
