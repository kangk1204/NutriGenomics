"""Source-specific composition adapters; nutrients are not exact isolated molecules."""
from __future__ import annotations
import csv
import hashlib
import json
import math
import zipfile
from decimal import Decimal, InvalidOperation
from pathlib import Path
from .graph import canonical, connect, edge, node, source, validate
from .sources import file_record, registry_hook, sha256

NUTRIENTS = {
 'ENERC':('1008','Energy','kcal'),'WATER':('1051','Water','g'),'PROT':('1003','Protein','g'),
 'FATCE':('1004','Total lipid','g'),'ASH':('1007','Ash','g'),'CHOCDF':('1005','Carbohydrate','g'),
 'SUGAR':('2000','Total sugars','g'),'FIBTG':('1079','Fiber','g'),'CA':('1087','Calcium','mg'),
 'FE':('1089','Iron','mg'),'P':('1091','Phosphorus','mg'),'K':('1092','Potassium','mg'),
 'NAT':('1093','Sodium','mg'),'VITA_RAE':('1106','Vitamin A RAE','ug'),'RETOL':('1105','Retinol','ug'),
 'CARTB':('1107','Beta-carotene','ug'),'THIA':('1165','Thiamin','mg'),'RIBF':('1166','Riboflavin','mg'),
 'NIA':('1167','Niacin','mg'),'VITC':('1162','Total vitamin C','mg'),'VITD':('1114','Vitamin D D2+D3','ug'),
 'CHOLE':('1253','Cholesterol','mg'),'FASAT':('1258','Saturated fatty acids','g'),'FATRN':('1257','Trans fatty acids','g')}


def parse_value(value) -> tuple[float | None,str]:
    if value is None or str(value).strip().casefold() in {'','na','n/a','-','미측정','미분석'}: return None,'missing'
    token=str(value).strip()
    if token.casefold() in {'tr','trace','nd','미량','불검출','미검출','<lod','<loq'} or token.startswith('<'): return None,'below_detection'
    try: amount=Decimal(token.replace(',',''))
    except InvalidOperation: return None,'unparsed'
    if not amount.is_finite() or amount<0: raise ValueError('invalid source nutrient value')
    return float(amount),'measured_zero' if amount==0 else 'measured'


def measurement(db, food: str, nutrient: str, name: str, amount, unit, basis, status, record, raw) -> None:
    node(db,nutrient,'nutrient',name)
    key=hashlib.sha256(canonical([food,nutrient,record['id'],basis]).encode()).hexdigest()
    db.execute('INSERT OR IGNORE INTO measurement VALUES (?,?,?,?,?,?,?,?,?)',(key,food,nutrient,amount,unit,basis,status,record['id'],canonical(raw)))
    edge(db,food,nutrient,'contains_nutrient','observed_composition',None,{'measurement_id':key,'basis':basis,'status':status,'source':record['source'],'amount_is_not_blood_concentration':True},record['id'])


def fdc_rows(path: Path, audit: dict | None = None):
    def unpack(payload):
        if isinstance(payload,list):foods=payload
        elif isinstance(payload,dict):
            wrappers=[key for key in ('FoundationFoods','SRLegacyFoods','SurveyFoods','FNDDSFoods') if key in payload]
            if len(wrappers)>1:raise ValueError('ambiguous FDC dataset wrapper')
            foods=payload[wrappers[0]] if wrappers else [payload]
        else:raise ValueError('invalid FDC JSON')
        if not isinstance(foods,list):raise ValueError('FDC foods must be a list')
        for food in foods:
            if food is None:
                if audit is not None:audit['null_food_entries']=audit.get('null_food_entries',0)+1
                continue
            if not isinstance(food,dict) or 'fdcId' not in food: raise ValueError('not an FDC food-detail snapshot')
            if food.get('dataType') not in {'SR Legacy','Foundation','Survey (FNDDS)'}: raise ValueError('unsupported FDC nutrient basis/data type')
            yield food
    if path.suffix.lower()=='.zip':
        # Read official downloadable JSON in place. No archive extraction or filename guessing.
        with zipfile.ZipFile(path) as archive:
            members=[member for member in archive.infolist() if not member.is_dir() and member.filename.lower().endswith('.json')]
            if not members:raise ValueError('FDC ZIP has no JSON')
            for member in members:
                with archive.open(member) as handle:
                    yield from unpack(json.load(handle))
    else:
        yield from unpack(json.loads(path.read_text(encoding='utf-8')))


def ingest_fdc(db,path,registry=None):
    record=file_record(path,'USDA-FDC','https://fdc.nal.usda.gov/',None,'CC0/public domain source; retain attribution and exact preparation')
    source(db,record); foods=0;audit={'null_food_entries':0,'invalid_quantity_entries':0}
    for food in fdc_rows(path,audit):
        identifier='FDC:'+str(food['fdcId']); foods+=1
        node(db,identifier,'food',food['description'],{'active_source':record['id'],'data_type':food['dataType'],'publication_date':food.get('publicationDate'),'food_attributes':food.get('foodAttributes',[]),'preparation':food['description']})
        for attribute in food.get('foodAttributes',[]):
            value=str(attribute.get('value',''))
            if '/FOODON_' in value:
                ontology='FOODON:'+value.rsplit('_',1)[-1];node(db,ontology,'food_ontology',ontology)
                edge(db,identifier,ontology,'source_foodon_annotation','curated',None,{'source':'USDA-FDC','source_attribute':attribute,'food_equivalence':False},record['id'])
        for row in food.get('foodNutrients',[]):
            nutrient=row.get('nutrient',{})
            try:
                amount,status=parse_value(row.get('amount'))
            except ValueError:
                # Preserve the official anomalous record without using it as a quantity.
                amount,status=None,'invalid_source_value'
                audit['invalid_quantity_entries']+=1
            derivation=row.get('foodNutrientDerivation',{});code=derivation.get('code')
            if amount is not None and code!='A':status='imputed' if code=='Z' else 'derived'
            unit=nutrient.get('unitName');unit={'µg':'ug','μg':'ug'}.get(unit,unit)
            measurement(db,identifier,'FDCNutrient:'+str(nutrient['id']),nutrient['name'],amount,unit,'source_100g',status,record,row)
    record['food_quality_audit']={**audit,'valid_food_records':foods}
    db.execute('UPDATE source SET metadata_json=? WHERE id=?',(canonical(record),record['id']))
    registry_hook(registry,record)
    return foods


def kfind_release(directory: Path):
    manifest_path=directory/'manifest.json';manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
    if str(manifest.get('dataset_id')) not in {'15100065','15100070'}:raise ValueError('unsupported K-FIND dataset')
    hashes=manifest.get('sha256',{})
    pages=sorted(directory.glob('page-*.json'))
    if not hashes or 'columns.json' not in hashes or len(pages)!=manifest.get('page_count'):raise ValueError('incomplete K-FIND source manifest')
    for filename,digest in hashes.items():
        if Path(filename).name!=filename or sha256(directory/filename)!=digest:raise ValueError('K-FIND hash mismatch or unsafe path')
    columns=json.loads((directory/'columns.json').read_text(encoding='utf-8')).get('tableVO',{}).get('colNmList',[])
    if not {'FOOD_CD','FOOD_NM','NUT_CON_SRTR_QUA',*NUTRIENTS}.issubset(columns):raise ValueError('missing K-FIND nutrient columns')
    rows=[];counts=[]
    for page in pages:
        if page.name not in hashes:raise ValueError('unhashed K-FIND page')
        batch=json.loads(page.read_text(encoding='utf-8'))
        if not isinstance(batch,list) or any(not isinstance(r,dict) for r in batch):raise ValueError('invalid K-FIND page')
        rows.extend(batch);counts.append(len(batch))
    if len(rows)!=manifest.get('total_count') or counts!=manifest.get('page_record_counts'):raise ValueError('K-FIND page/row count mismatch')
    return manifest,rows


def ingest_kfind(db,inputs,selection_path,registry=None):
    selection=json.loads(selection_path.read_text(encoding='utf-8'))['foods']
    if not selection:raise ValueError('empty selection')
    selected_keys=[(str(s['dataset_id']),s['food_code'],s['food_name']) for s in selection]
    if len(set(selected_keys))!=len(selected_keys):raise ValueError('duplicate selected food')
    releases={}
    for directory in inputs:
        manifest,rows=kfind_release(directory);dataset_id=str(manifest['dataset_id'])
        if dataset_id in releases:raise ValueError('duplicate K-FIND release dataset')
        record=file_record(directory/'manifest.json','K-FIND','https://www.data.go.kr/data/'+dataset_id+'/standard.do',manifest.get('provider_reference_date_max'),'Source-specific reuse conditions unresolved; private research')
        record['input_page_sha256']=manifest['sha256'];source(db,record)
        releases[dataset_id]=(manifest,rows,record)
    for item in selection:
        dataset_id=str(item['dataset_id'])
        if dataset_id not in releases:raise ValueError('selected K-FIND release not supplied')
        _,rows,record=releases[dataset_id]
        matches=[r for r in rows if str(r['FOOD_CD'])==item['food_code'] and str(r['FOOD_NM'])==item['food_name']]
        if len(matches)!=1:raise ValueError('selected food code/name not uniquely matched')
        row=matches[0]
        basis=str(row.get('NUT_CON_SRTR_QUA','')).replace(' ','').lower()
        if basis not in {'100g','100.0g'} or row.get('DATA_PROD_NM')!='분석':raise ValueError('only verified 100g analyzed K-FIND foods are admitted')
        identifier=f'K-FIND:{dataset_id}:{item["food_code"]}:'+hashlib.sha256(item['food_name'].encode()).hexdigest()[:12]
        node(db,identifier,'food',item['name'],{'active_source':record['id'],'food_code':item['food_code'],'source_name':item['food_name'],'preparation':item['preparation_state'],'REFUSE':row.get('REFUSE'),'refuse_conversion_applied':False})
        for field,(nutrient,name,unit) in NUTRIENTS.items():
            amount,status=parse_value(row.get(field))
            measurement(db,identifier,'FDCNutrient:'+nutrient,name,amount,unit,'source_100g',status,record,{'field':field,'raw_value':row.get(field),'food_code':item['food_code']})
    for _,_,record in releases.values():registry_hook(registry,record)
    return len(selection)


def ingest_rda(db,path,mapping_path,registry=None):
    """Explicitly reviewed column mapping; never guess Excel headers or nutrient units."""
    mapping=json.loads(mapping_path.read_text(encoding='utf-8'))
    if mapping.get('basis')!='source_100g' or mapping.get('basis_verified') is not True:raise ValueError('RDA 100g basis must be reviewed')
    record=file_record(path,'RDA',mapping['source_url'],mapping['version'],'RDA source-specific use conditions must be retained')
    source(db,record)
    if path.suffix.lower()=='.xlsx':
        import openpyxl
        workbook=openpyxl.load_workbook(path,read_only=True,data_only=True)
        sheet=workbook[mapping['sheet']];iterator=sheet.iter_rows(values_only=True)
        for _ in range(mapping.get('header_row',1)-1):next(iterator)
        columns=[str(x) for x in next(iterator)];rows=(dict(zip(columns,r)) for r in iterator)
    elif path.suffix.lower()=='.csv':
        handle=path.open(encoding=mapping.get('encoding','utf-8-sig'),newline='');rows=csv.DictReader(handle)
    else:raise ValueError('RDA input must be csv or xlsx')
    count=0
    try:
        for row in rows:
            if not row.get(mapping['food_code']):continue
            identifier='RDA:'+mapping['version']+':'+str(row[mapping['food_code']]);count+=1
            node(db,identifier,'food',str(row[mapping['food_name']]),{'active_source':record['id'],'preparation':row.get(mapping.get('preparation')),'version':mapping['version']})
            for column,n in mapping['nutrients'].items():
                if column not in row:raise ValueError('mapped RDA nutrient column missing')
                amount,status=parse_value(row[column]);measurement(db,identifier,n['id'],n['name'],amount,n['unit'],'source_100g',status,record,{'field':column,'raw_value':row[column]})
    finally:
        if path.suffix.lower()=='.xlsx':workbook.close()
        else:handle.close()
    registry_hook(registry,record);return count


def ingest_food(database:Path,format:str,inputs:list[Path],selection=None,mapping=None,registry=None):
    db=connect(database)
    try:
        with db:
            if format=='fdc':count=sum(ingest_fdc(db,p,registry) for p in inputs)
            elif format=='kfind':
                if selection is None:raise ValueError('K-FIND requires --selection')
                count=ingest_kfind(db,inputs,selection,registry)
            elif format=='rda':
                if mapping is None:raise ValueError('RDA requires --mapping')
                count=sum(ingest_rda(db,p,mapping,registry) for p in inputs)
            else:raise ValueError('unknown food adapter')
            result=validate(db)
            if result['errors']:raise ValueError(result['errors'])
        result['foods_ingested']=count;return result
    finally:db.close()


def meal(db,foods:list[dict]) -> dict:
    totals={};missing=[]
    for food in foods:
        grams=float(food['grams'])
        if not math.isfinite(grams) or grams<=0:raise ValueError('grams must be finite and positive')
        metadata=db.execute("SELECT metadata_json FROM node WHERE id=? AND kind='food'",(food['food_id'],)).fetchone()
        if not metadata:raise ValueError('unknown food')
        active=json.loads(metadata[0]).get('active_source')
        if not active:raise ValueError('food active release unresolved')
        for row in db.execute('SELECT * FROM measurement WHERE food=? AND source_id=?',(food['food_id'],active)):
            if row['basis']!='source_100g' or row['status'] not in {'measured','measured_zero'} or row['amount'] is None:
                missing.append({'food':food['food_id'],'nutrient':row['nutrient'],'status':row['status']});continue
            key=(row['nutrient'],row['unit']);totals[key]=totals.get(key,0)+row['amount']*grams/100
    return {'totals':[{'nutrient':key[0],'unit':key[1],'amount':amount} for key,amount in sorted(totals.items())], 'excluded_or_unknown':missing,'blood_concentration_predicted':False,'causal_health_effect_predicted':False}
