import json
import math
import zipfile
import pytest
from nutriomics_evidence.food import ingest_food,parse_value,meal,NUTRIENTS
from nutriomics_evidence.graph import connect
from nutriomics_evidence.identity import validate_mapping
from nutriomics_evidence.sources import sha256

def test_missing_detection_zero_separated():
    assert parse_value(None)==(None,'missing')
    assert parse_value('0')==(0.,'measured_zero')
    assert parse_value('<0.01')==(None,'below_detection')
    assert parse_value('Tr')==(None,'below_detection')
    with pytest.raises(ValueError):parse_value('-2')
    with pytest.raises(ValueError):parse_value('NaN')

def test_fdc_imputed_zero_is_excluded_and_foodon_not_structure(tmp_path):
    path=tmp_path/'fdc.json';path.write_text(json.dumps({'fdcId':1,'dataType':'SR Legacy','description':'Raw food','foodAttributes':[{'value':'http://purl.obolibrary.org/obo/FOODON_03301710'}],'foodNutrients':[
       {'nutrient':{'id':1,'name':'A','unitName':'g'},'amount':2,'foodNutrientDerivation':{'code':'A'}},
       {'nutrient':{'id':2,'name':'B','unitName':'mg'},'amount':0,'foodNutrientDerivation':{'code':'Z'}},
       {'nutrient':{'id':3,'name':'C','unitName':'mg'}}]}))
    database=tmp_path/'db.sqlite';ingest_food(database,'fdc',[path]);db=connect(database,readonly=True)
    result=meal(db,[{'food_id':'FDC:1','grams':200}]);assert result['totals'][0]['amount']==4
    assert len(result['excluded_or_unknown'])==2
    assert db.execute("SELECT relation FROM edge WHERE object='FOODON:03301710'").fetchone()[0]=='source_foodon_annotation'
    assert db.execute('SELECT count(*) FROM identifier_mapping').fetchone()[0]==0
    db.close()

def test_name_mapping_cannot_become_exact_structure():
    row={'subject':'MESH:C1','namespace':'PubChem','identifier':'1','relation':'exact_structure','basis':'same name','method':'same_name','reviewed':True}
    with pytest.raises(ValueError):validate_mapping(row)
    row['relation']='candidate_name';assert validate_mapping(row)==row

def test_exact_protonation_key_mismatch_rejected():
    row={'subject':'MESH:C1','namespace':'PubChem','identifier':'1','relation':'exact_structure','basis':'keys','method':'same_standard_inchikey','reviewed':True,'subject_inchikey':'AAAAAAAAAAAAAA-BBBBBBBBBB-C','object_inchikey':'AAAAAAAAAAAAAA-BBBBBBBBBB-D'}
    with pytest.raises(ValueError):validate_mapping(row)

def test_new_fdc_source_updates_active_release_without_double_count(tmp_path):
    database=tmp_path/'db.sqlite'
    for amount in [2,3]:
        path=tmp_path/f'{amount}.json';path.write_text(json.dumps({'fdcId':1,'dataType':'SR Legacy','description':'Food','foodNutrients':[{'nutrient':{'id':1,'name':'N','unitName':'g'},'amount':amount,'foodNutrientDerivation':{'code':'A'}}]}));ingest_food(database,'fdc',[path])
    db=connect(database,readonly=True);assert meal(db,[{'food_id':'FDC:1','grams':100}])['totals'][0]['amount']==3;db.close()


def test_official_fdc_zip_wrappers(tmp_path):
    path=tmp_path/'survey.zip'
    with zipfile.ZipFile(path,'w') as archive:
        archive.writestr('survey.json',json.dumps({'SurveyFoods':[{'fdcId':1,'dataType':'Survey (FNDDS)','description':'Cooked food','foodNutrients':[]}]}))
    result=ingest_food(tmp_path/'db.sqlite','fdc',[path]);assert result['foods_ingested']==1


@pytest.fixture
def kfind_fixture(tmp_path):
    directory=tmp_path/'release';directory.mkdir()
    row={key:'0' for key in NUTRIENTS}
    row.update({'FOOD_CD':'R1','FOOD_NM':'쌀밥','NUT_CON_SRTR_QUA':'100g','DATA_PROD_NM':'분석','NAT':'3.5','VITC':'Tr'})
    (directory/'page-0001.json').write_text(json.dumps([row],ensure_ascii=False),encoding='utf-8')
    (directory/'columns.json').write_text(json.dumps({'tableVO':{'colNmList':list(row)}}))
    manifest={'dataset_id':'15100065','page_count':1,'page_record_counts':[1],'total_count':1,'sha256':{name:sha256(directory/name) for name in ['page-0001.json','columns.json']}}
    (directory/'manifest.json').write_text(json.dumps(manifest))
    selection=tmp_path/'selection.json';selection.write_text(json.dumps({'foods':[{'dataset_id':'15100065','food_code':'R1','food_name':'쌀밥','name':'밥','preparation_state':'조리 후'}]},ensure_ascii=False),encoding='utf-8')
    return directory,selection,row,manifest


def test_kfind_identity_basis_zero_detection_and_source_hash(kfind_fixture,tmp_path):
    directory,selection,_,_=kfind_fixture;database=tmp_path/'db.sqlite'
    result=ingest_food(database,'kfind',[directory],selection);assert result['foods_ingested']==1
    db=connect(database,readonly=True)
    assert db.execute('SELECT count(*) FROM measurement').fetchone()[0]==24
    assert db.execute("SELECT amount,status FROM measurement WHERE nutrient='FDCNutrient:1162'").fetchone()[:]==(None,'below_detection')
    assert db.execute("SELECT amount,status FROM measurement WHERE nutrient='FDCNutrient:1093'").fetchone()[:]==(3.5,'measured')
    db.close()


@pytest.mark.parametrize('damage',['hash','basis','name'])
def test_kfind_rejects_unverified_input(kfind_fixture,tmp_path,damage):
    directory,selection,row,manifest=kfind_fixture
    if damage=='name':
        payload=json.loads(selection.read_text(encoding='utf-8'));payload['foods'][0]['food_name']='different name';selection.write_text(json.dumps(payload))
    else:
        row['NUT_CON_SRTR_QUA']='1 serving';page=directory/'page-0001.json';page.write_text(json.dumps([row],ensure_ascii=False),encoding='utf-8')
        if damage=='basis':
            manifest['sha256'][page.name]=sha256(page);(directory/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError):ingest_food(tmp_path/'db.sqlite','kfind',[directory],selection)
