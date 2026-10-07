"""Strict, isolated public-data staging with explicit counting units.

Natural-product membership, edible-organism occurrence, and measured food
composition are deliberately separate evidence types. No production DB writer.
"""
from __future__ import annotations
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import shutil
import sqlite3
import zipfile
import duckdb
import pandas as pd

FULL_KEY = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")

def digest(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def verified_source(root: Path, sid: str) -> dict:
    receipt=json.loads((root/"receipts"/f"{sid}.json").read_text(encoding="utf-8"))
    if receipt.get("status")!="complete":raise ValueError(f"incomplete source {sid}")
    path=root/"raw"/receipt["path"].replace("\\","/").rsplit("/",1)[-1]
    if not path.exists() or digest(path)!=receipt["sha256"]:raise ValueError(f"source SHA mismatch {sid}")
    return receipt|{"actual_path":str(path)}

def extract_csv(root: Path, sid: str, requested: set[str]|None=None) -> list[Path]:
    receipt=verified_source(root,sid);archive_path=Path(receipt["actual_path"])
    destination=root/"extracted"/sid;destination.mkdir(parents=True,exist_ok=True)
    paths=[]
    with zipfile.ZipFile(archive_path) as archive:
        for info in archive.infolist():
            name=Path(info.filename).name
            if info.is_dir() or not name.endswith(".csv") or (requested and name not in requested):continue
            path=destination/name
            if not path.exists() or path.stat().st_size!=info.file_size:
                temporary=path.with_suffix(".csv.part")
                with archive.open(info) as source,temporary.open("wb") as output:shutil.copyfileobj(source,output,8*1024*1024)
                temporary.replace(path)
            paths.append(path)
    return paths

def sql_literal(value: str|Path) -> str:return "'"+str(value).replace("'","''")+"'"

def mass_exceeds_basis_sql(amount: str="amount", unit: str="unit") -> str:
    """Absolute physical bound for mass per100g; never sum overlapping nutrients."""
    return f"""CASE upper(trim({unit}))
      WHEN 'G' THEN {amount}>100
      WHEN 'MG' THEN {amount}>100000
      WHEN 'UG' THEN {amount}>100000000
      WHEN 'MCG' THEN {amount}>100000000
      ELSE false END"""

def csv_relation(path: Path) -> str:
    return f"read_csv({sql_literal(path)}, header=true, all_varchar=true, delim=',', quote='\"', escape='\"', strict_mode=true, null_padding=false, max_line_size=16777216)"

def materialize(con: duckdb.DuckDBPyConnection, derived: Path, name: str, query: str) -> int:
    target=derived/f"{name}.parquet"
    con.execute(f"COPY ({query}) TO {sql_literal(target)} (FORMAT PARQUET, COMPRESSION ZSTD)")
    con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet({sql_literal(target)})")
    return con.execute(f"SELECT count(*) FROM {name}").fetchone()[0]

def stage_fdc(con, root: Path, summary: dict, baseline: Path|None):
    paths={p.name:p for p in extract_csv(root,"usda_fdc_202604",{"food.csv","branded_food.csv","food_nutrient.csv","nutrient.csv"})}
    for name,file in [("fdc_native_food","food.csv"),("fdc_native_branded","branded_food.csv"),("fdc_native_nutrient","nutrient.csv"),("fdc_native_composition","food_nutrient.csv")]:
        con.execute(f"CREATE OR REPLACE TABLE {name} AS SELECT * FROM {csv_relation(paths[file])}")
    con.execute("""CREATE OR REPLACE TABLE food_candidates AS
    SELECT 'FDC:'||f.fdc_id AS food_id,'usda_fdc_202604' AS source_id,
      f.fdc_id AS source_food_id,f.description,f.data_type,f.publication_date,
      CASE WHEN f.data_type='branded_food' AND b.gtin_upc IS NOT NULL AND regexp_full_match(b.gtin_upc,'[0-9]{8,14}')
      THEN 'GTIN:'||lpad(b.gtin_upc,14,'0')||':'||coalesce(b.market_country,'unknown_market')
      ELSE 'FDC:'||f.fdc_id END AS identity_key,
      b.gtin_upc,b.brand_owner,b.market_country,
      row_number() OVER(PARTITION BY f.fdc_id ORDER BY f.publication_date DESC) AS source_id_rank
    FROM fdc_native_food f LEFT JOIN fdc_native_branded b USING(fdc_id)
    WHERE try_cast(f.fdc_id AS BIGINT)>0 AND length(trim(f.description))>0""")
    derived=root/"derived"
    native_n=materialize(con,derived,"food_records","SELECT * EXCLUDE(source_id_rank) FROM food_candidates WHERE source_id_rank=1")
    food_n=materialize(con,derived,"food","""SELECT * EXCLUDE(identity_rank) FROM
      (SELECT *,row_number() OVER(PARTITION BY identity_key ORDER BY try_cast(publication_date AS DATE) DESC NULLS LAST,try_cast(source_food_id AS BIGINT) DESC) identity_rank FROM food_records
       WHERE data_type IN ('branded_food','sr_legacy_food','foundation_food','survey_fndds_food'))
      WHERE identity_rank=1""")
    invalid=con.execute("""SELECT count(*) FROM fdc_native_composition c
      LEFT JOIN fdc_native_nutrient n ON c.nutrient_id=n.id
      WHERE try_cast(c.amount AS DOUBLE) IS NULL OR try_cast(c.amount AS DOUBLE)<0 OR n.id IS NULL OR n.unit_name IS NULL""").fetchone()[0]
    con.execute("""CREATE OR REPLACE TEMP VIEW typed_food_composition AS SELECT DISTINCT f.food_id,c.nutrient_id,n.name nutrient_name,
      try_cast(c.amount AS DOUBLE) amount,n.unit_name unit,'source100g' basis,c.id source_row_id,'usda_fdc_202604' source_id,
      c.derivation_id,c.data_points,c.min,c.max,c.median,c.loq,c.footnote,c.min_year_acquired
      FROM fdc_native_composition c JOIN food f ON f.source_food_id=c.fdc_id JOIN fdc_native_nutrient n ON c.nutrient_id=n.id
      WHERE try_cast(c.amount AS DOUBLE)>=0 AND isfinite(try_cast(c.amount AS DOUBLE)) AND n.unit_name IS NOT NULL""")
    rejected_n=materialize(con,derived,"food_composition_rejections","SELECT *, 'mass_exceeds100g_basis' rejection_reason FROM typed_food_composition WHERE "+mass_exceeds_basis_sql())
    composition_n=materialize(con,derived,"food_composition","SELECT * FROM typed_food_composition WHERE NOT "+mass_exceeds_basis_sql())
    checks={"duplicate_food_id":con.execute("SELECT count(*)-count(DISTINCT food_id) FROM food").fetchone()[0],
            "duplicate_identity_key":con.execute("SELECT count(*)-count(DISTINCT identity_key) FROM food").fetchone()[0],
            "composition_orphan_food":con.execute("SELECT count(*) FROM food_composition c ANTI JOIN food f USING(food_id)").fetchone()[0],
            "invalid_amount_or_unit_native_rows":invalid,"physically_impossible_mass_rows_rejected":rejected_n}
    delta={"baseline_reported_food_records":5845,"baseline_ids_available":False,"exact_added_source_records":None,"exact_added_canonical_food":None}
    if baseline is not None:
        db=sqlite3.connect(baseline.resolve().as_uri()+"?mode=ro",uri=True)
        if db.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='node'").fetchone()[0]:
            baseline_query="SELECT id,name,metadata_json FROM node WHERE kind='food'"
        else:
            baseline_query="SELECT id,label AS name,metadata_json FROM ng_entity WHERE kind='food'"
        snapshot=pd.read_sql_query(baseline_query,db);db.close()
        snapshot.to_parquet(derived/"baseline_food_ids.parquet",index=False)
        con.register("baseline_frame",snapshot)
        delta.update(baseline_ids_available=True,baseline_actual_food_records=len(snapshot),baseline_sha256=digest(baseline),
          exact_added_source_records=con.execute("SELECT count(*) FROM food_records f ANTI JOIN baseline_frame b ON f.food_id=b.id").fetchone()[0],
          exact_added_canonical_food=con.execute("SELECT count(*) FROM food f ANTI JOIN baseline_frame b ON f.food_id=b.id").fetchone()[0],
          delta_unit="new FDC source IDs relative to frozen existing graph; food identity dedup within new source separately")
    summary["food"]={"source_records":native_n,"canonical_product_or_food_identities":food_n,"composition_records":composition_n,"composition_rejected_mass_records":rejected_n,"mass_policy":"source100g mass <=100g,100000mg,100000000ug; IU/activity units retained natively without mass inference; overlapping nutrients never summed","food_type_counts":dict(con.execute("SELECT data_type,count(*) FROM food GROUP BY data_type").fetchall()),"food_registry_eligible_types":["branded_food","sr_legacy_food","foundation_food","survey_fndds_food"],"lab_sample_types_counted_as_food_entities":False,"deduplication":"branded GTIN padded14+market, latest publication; nonbranded nativeFDCID; no name-only entity merge","delta":delta,"checks":checks}
    con.execute("DROP VIEW typed_food_composition")
    for table in ["food_candidates","fdc_native_food","fdc_native_branded","fdc_native_nutrient","fdc_native_composition"]:con.execute(f"DROP TABLE {table}")

def structure_check(row: tuple) -> tuple:
    from rdkit.Chem import inchi
    row=tuple(v if isinstance(v,str) else None for v in row)
    identifier,smiles,standard_inchi,full_key,name,organisms,collections,dois=row
    valid=bool(full_key and FULL_KEY.fullmatch(full_key))
    computed=inchi.InchiToInchiKey(standard_inchi) if standard_inchi and standard_inchi.startswith("InChI=1S/") else None
    match=valid and computed==full_key
    single=bool(smiles and "." not in smiles)
    return (identifier,smiles,standard_inchi,full_key,name,organisms,collections,dois,valid,computed,match,single,
            "source_asserted_food_database_membership" if collections and "FooDB" in collections.split("|") else "natural_product_only")

def stage_coconut(con, root: Path, summary: dict):
    path=extract_csv(root,"coconut_202610")[0]
    frame=con.execute(f"SELECT identifier,canonical_smiles,standard_inchi,standard_inchi_key,name,organisms,collections,dois FROM {csv_relation(path)}").fetchdf()
    checked=[structure_check(row) for row in frame.itertuples(index=False,name=None)]
    frame=pd.DataFrame(checked,columns=["source_compound_id","smiles","standard_inchi","full_inchikey","name","organisms","collections","dois","full_key_format_valid","computed_inchikey","inchi_hash_verified","single_component","food_membership_kind"])
    frame["source_id"]="coconut_202610";frame.to_parquet(root/"derived"/"natural_product_records.parquet",index=False)
    con.execute(f"CREATE OR REPLACE VIEW natural_product_records AS SELECT * FROM read_parquet({sql_literal(root/'derived/natural_product_records.parquet')})")
    count=materialize(con,root/"derived","natural_products","""SELECT * EXCLUDE(rn) FROM
      (SELECT *,row_number() OVER(PARTITION BY full_inchikey ORDER BY source_compound_id) rn FROM natural_product_records WHERE inchi_hash_verified AND single_component) WHERE rn=1""")
    summary["natural_products"]={"native_records":len(frame),"distinct_verified_full_inchikey_single_structures":count,
      "inchi_hash_verified_records":int(frame.inchi_hash_verified.sum()),"single_component_records":int(frame.single_component.sum()),
      "food_database_membership_records":int((frame.food_membership_kind=="source_asserted_food_database_membership").sum()),
      "verification":"source standardInChI produces identical fullInChIKey; full stereo key retained; disconnectedSMILES excluded; source natural-product annotation retained",
      "experimental_food_occurrence_proven_by_collection_membership":False,"all_smiles_recomputed_and_stereochemistry_validated":False}

def stage_lotus(con, root: Path, summary: dict):
    receipt=verified_source(root,"lotus_202301")
    path=Path(receipt["actual_path"])
    rows=materialize(con,root/"derived","natural_product_occurrence","""SELECT DISTINCT structure_inchikey full_inchikey,organism_name,organism_wikidata,reference_doi,reference_wikidata,
      structure_wikidata,manual_validation,'lotus_202301' source_id,'natural_organism_occurrence' evidence_kind
      FROM """+csv_relation(path)+" WHERE regexp_full_match(structure_inchikey,'[A-Z]{14}-[A-Z]{10}-[A-Z]')")
    mapped=con.execute("SELECT count(*) FROM natural_product_occurrence o JOIN natural_products n USING(full_inchikey)").fetchone()[0]
    summary["lotus"]={"distinct_occurrence_records":rows,"exact_fullkey_mapped_occurrences":mapped,"food_occurrence_classified_without_edible_part_or_food_proof":False,"all_natural_occurrences_promoted_to_food":False}

def merge_ffq_cycle(frames: dict[str,pd.DataFrame], cycle: str) -> tuple[pd.DataFrame,pd.DataFrame]:
    demo=frames["DEMO"].copy();raw=frames["FFQRAW"].copy();bpx=frames["BPX"].copy();bpq=frames["BPQ"].copy()
    for name,frame in [("DEMO",demo),("FFQRAW",raw),("BPX",bpx),("BPQ",bpq)]:
        if frame.SEQN.duplicated().any():raise ValueError(f"duplicate person ID in {name}")
    required=["SEQN","SDMVPSU","SDMVSTRA","WTMEC2YR","WTINT2YR","RIDAGEYR","RIAGENDR","RIDRETH1"]
    people=demo[required].merge(raw[[c for c in ["SEQN","WTS_FFQ","FFQ_MISS"] if c in raw]],on="SEQN",how="left",validate="one_to_one")
    if "WTS_FFQ" not in people:raise ValueError("dedicated FFQ survey weights missing")
    bp_cols=[c for c in bpx if re.fullmatch(r"BPXSY[1-4]|BPXDI[1-4]",c)]
    people=people.merge(bpx[["SEQN"]+bp_cols],on="SEQN",how="left",validate="one_to_one")
    people=people.merge(bpq[[c for c in ["SEQN","BPQ020","BPQ040A","BPQ050A"] if c in bpq]],on="SEQN",how="left",validate="one_to_one")
    people["cycle"]=cycle;people["participant_id"]=cycle+":"+people.SEQN.astype("int64").astype(str)
    for outcome,prefix in [("observed_sbp_mean","BPXSY"),("observed_dbp_mean","BPXDI")]:
        people[outcome]=people[[c for c in bp_cols if c.startswith(prefix)]].where(lambda x:x>0).mean(axis=1)
    valid=people.observed_sbp_mean.notna()&people.observed_dbp_mean.notna()
    people["observed_bp_ge140_90"] = pd.Series(pd.NA,index=people.index,dtype="boolean")
    people.loc[valid,"observed_bp_ge140_90"]=(people.loc[valid,"observed_sbp_mean"]>=140)|(people.loc[valid,"observed_dbp_mean"]>=90)
    people["self_reported_htn_diagnosis"]=people.BPQ020.map({1.0:True,2.0:False}).astype("boolean")
    long=frames["FFQDC"].copy();long["cycle"]=cycle;long["participant_id"]=cycle+":"+long.SEQN.astype("int64").astype(str)
    keys=["participant_id","FFQ_VAR","FFQ_FOOD"]
    if long.duplicated(keys).any():raise ValueError("duplicated native FFQ person-variable-food keys")
    if not long.participant_id.isin(people.participant_id).all():raise ValueError("FFQ person missing DEMO")
    people["has_ffq_frequency"]=people.participant_id.isin(long.participant_id)
    long["frequency_unit"]="times_per_day";long["portion_size_available"]=False;long["absolute_nutrient_intake_derived"]=False
    long["source_id"]="nhanes_FFQDC_"+{"2003-2004":"C","2005-2006":"D"}[cycle]
    long["source_row_id"]=long.participant_id+":"+long.FFQ_VAR.astype("int64").astype(str)+":"+long.FFQ_FOOD.astype("int64").astype(str)
    return people,long

def stage_ffq(con,root:Path,summary:dict):
    people_all=[];long_all=[];cycles=[]
    for cycle,suffix in [("2003-2004","C"),("2005-2006","D")]:
        frames={}
        for name in ["DEMO","FFQRAW","FFQDC","BPX","BPQ","FOODLK","VARLK"]:
            receipt=verified_source(root,f"nhanes_{name}_{suffix}")
            frames[name]=pd.read_sas(receipt["actual_path"],format="xport",encoding="utf-8")
            frames[name].to_parquet(root/"derived"/f"nhanes_native_{name}_{suffix}.parquet",index=False)
        people,long=merge_ffq_cycle(frames,cycle)
        people_all.append(people);long_all.append(long)
        cycles.append({"cycle":cycle,"demographic_people":len(people),"ffq_people":long.participant_id.nunique(),"frequency_rows":len(long),"positive_ffq_weight_people":int((people.WTS_FFQ>0).sum()),"ffq_with_observed_bp_people":int(people.loc[people.participant_id.isin(long.participant_id),"observed_bp_ge140_90"].notna().sum())})
    for name,frame in [("ffq_participant",pd.concat(people_all,ignore_index=True)),("ffq_frequency",pd.concat(long_all,ignore_index=True))]:
        frame.to_parquet(root/"derived"/f"{name}.parquet",index=False)
        con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet({sql_literal(root/'derived'/f'{name}.parquet')})")
    summary["ffq"]={"cycles":cycles,"participant_mapping":"cycle+SEQN within same NHANES cycle only","weights":"native WTS_FFQ for FFQ; WTMEC2YR/WTINT2YR retained; SDMVPSU+SDMVSTRA retained; no pooled or survey-inference analysis performed","bp_label":"arithmetic mean of positive observed BPXSY1-4/BPXDI1-4; observed>=140/90 separate from BPQ020 diagnosis; not clinical repeated-visit diagnosis","portion_size_available":False,"absolute_nutrient_intake_derived":False,"paired_omics_claim":False}

def stage_ctd(con,root:Path,summary:dict):
    path=root/"raw"/"CTD_chem_gene_ixns.csv.gz"
    if not path.exists():summary["ctd"]={"status":"not_staged"};return
    fields=None;skip=0
    with gzip.open(path,"rt",encoding="utf-8",newline="") as stream:
        for line in stream:
            if line.startswith("# ChemicalName,"):fields=line[2:].strip().split(",")
            if line.strip() and not line.startswith("#"):break
            skip+=1
    if not fields:raise ValueError("CTD field schema absent")
    columns="{"+",".join(sql_literal(f)+":'VARCHAR'" for f in fields)+"}"
    relation=f"read_csv({sql_literal(path)},header=false,skip={skip},columns={columns},delim=',',strict_mode=true,null_padding=false,max_line_size=16777216)"
    key="md5(to_json(["+",".join('coalesce("'+f+'",\'\')' for f in fields)+"]))"
    rows=materialize(con,root/"derived","chemical_gene_interaction",f"SELECT DISTINCT {key} source_row_id,*, 'ctd_20260929' source_id,'curated_molecular_response' evidence_kind FROM {relation} WHERE ChemicalID IS NOT NULL AND GeneID IS NOT NULL AND PubMedIDs IS NOT NULL")
    summary["ctd"]={"source_sha256":digest(path),"distinct_curated_interaction_rows":rows,"human_rows":con.execute("SELECT count(*) FROM chemical_gene_interaction WHERE OrganismID='9606'").fetchone()[0],"native_missing_required_fields":con.execute(f"SELECT count(*) FROM {relation} WHERE ChemicalID IS NULL OR GeneID IS NULL OR PubMedIDs IS NULL").fetchone()[0],"unit":"distinct original ChemicalID/GeneID/organism/action/PMID-rich source record; all native columns in dedup hash","physical_binding_assumed":False,"inference_assumed":False}

def sqlite_lookup(con,root:Path,summary:dict):
    path=root/"derived"/"expansion_lookup.sqlite"
    db=sqlite3.connect(path)
    db.execute("CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value_json TEXT NOT NULL)")
    db.execute("INSERT OR REPLACE INTO metadata VALUES('summary',?)",(json.dumps(summary,ensure_ascii=False),))
    for name in ["food","ffq_participant"]:
        columns=con.execute(f"DESCRIBE {name}").fetchall()
        db.execute(f"DROP TABLE IF EXISTS {name}")
        definitions=','.join('"'+c[0]+'" TEXT' for c in columns)
        db.execute(f"CREATE TABLE {name}({definitions})")
        cursor=con.execute(f"SELECT * FROM {name}")
        while batch:=cursor.fetchmany(20000):db.executemany(f"INSERT INTO {name} VALUES({','.join(['?']*len(columns))})",batch)
        db.execute(f"CREATE INDEX {name}_id ON {name}({ 'food_id' if name=='food' else 'participant_id'})")
    db.commit();summary["sqlite_integrity"]=db.execute("PRAGMA integrity_check").fetchone()[0];db.close()

def register_genotypes(con,root:Path,summary:dict):
    receipt_path=root/"genotype_reference_receipt.json"
    if not receipt_path.exists():return
    receipt=json.loads(receipt_path.read_text())
    if receipt.get("status")!="complete":summary["genomics_reference"]={"status":"failed","error":receipt.get("error")};return
    for name in ["genotype_reference_call","genotype_reference_sample"]:
        path=root/"derived"/f"{name}.parquet"
        expected=next(item["sha256"] for item in receipt["outputs"] if Path(item["path"]).name==path.name)
        if digest(path)!=expected:raise ValueError("reference genotype output hash mismatch")
        con.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet({sql_literal(path)})")
    summary["genomics_reference"]={"status":"complete","variants":receipt["variants"],"individuals":receipt["vcf_samples"],"genotype_calls":receipt["genotype_calls"],"build":"GRCh37","region":receipt["region"],"reference_context_only":True,"individual_genotype_not_gwas_summary":True,"diet_person_linkage":False,"source_receipt_sha256":digest(receipt_path)}

def run(root:Path,threads:int=2,baseline:Path|None=None)->dict:
    if threads not in range(1,5):raise ValueError("threads must be1..4")
    (root/"derived").mkdir(parents=True,exist_ok=True)
    output=root/"derived"/"expansion.duckdb"
    con=duckdb.connect(str(output));con.execute(f"SET threads={threads}");con.execute("SET memory_limit='12GB'");con.execute("SET enable_progress_bar=false")
    summary={"schema_version":1,"generated_utc":datetime.now(timezone.utc).isoformat(),"production_db_modified":False,"counting_policy":"typed source records counted once within each native namespace; heterogeneous units not silently summed"}
    for name,function in [("food",lambda:stage_fdc(con,root,summary,baseline)),("natural_products",lambda:stage_coconut(con,root,summary)),("lotus",lambda:stage_lotus(con,root,summary)),("ffq",lambda:stage_ffq(con,root,summary)),("ctd",lambda:stage_ctd(con,root,summary))]:
        print(json.dumps({"stage":name,"status":"running"}),flush=True)
        function();(root/"derived"/"summary.partial.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
        print(json.dumps({"stage":name,"status":"complete"}),flush=True)
    # A membership flag is not a validated food occurrence or a measured amount.
    summary["verified_food_single_compounds"]={"count":0,"target":2000,"status":"not_met","reason":"FooDB bulk403; COCONUT FooDB membership and LOTUS organism occurrence do not establish edible part/measured food occurrence"}
    register_genotypes(con,root,summary)
    sqlite_lookup(con,root,summary)
    con.execute("CHECKPOINT");con.close()
    files=[{"path":str(p),"sha256":digest(p),"bytes":p.stat().st_size} for p in sorted((root/"derived").glob("*")) if p.suffix in {".parquet",".duckdb",".sqlite"}]
    summary["outputs"]=files
    (root/"derived"/"summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding="utf-8")
    return summary
