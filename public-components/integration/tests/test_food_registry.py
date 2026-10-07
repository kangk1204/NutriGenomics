import pytest
duckdb=pytest.importorskip('duckdb')
from fastapi.testclient import TestClient
from nutriomics_final.api import create_app
from nutriomics_final.food_registry import FoodRegistry


def make_database(path):
    with duckdb.connect(str(path)) as c:
        c.execute('CREATE TABLE food(food_id VARCHAR,source_id VARCHAR,source_food_id VARCHAR,description VARCHAR,data_type VARCHAR,publication_date VARCHAR,identity_key VARCHAR)')
        c.execute("INSERT INTO food VALUES ('FDC:1','USDA','1','Green tea','Foundation','2026-04','USDA:1')")
        c.execute('CREATE TABLE food_composition(food_id VARCHAR,nutrient_id VARCHAR,nutrient_name VARCHAR,amount DOUBLE,unit VARCHAR,basis VARCHAR,source_row_id VARCHAR,source_id VARCHAR)')
        c.execute("INSERT INTO food_composition VALUES ('FDC:1','1','Sodium',4,'mg','source100g','A','USDA'),('FDC:1','2','Protein',5,'g','portion','B','USDA'),('FDC:1','3','Unknown',NULL,'mg','source100g','C','USDA')")


def test_actual_portion_scaling_and_incompatible_units(tmp_path):
    path=tmp_path/'food.duckdb';make_database(path)
    with TestClient(create_app(tmp_path,food_database=path)) as client:
        assert client.get('/foods',params={'q':"' OR 1=1 --"}).json()['items']==[]
        assert client.get('/foods',params={'q':'tea'}).json()['items'][0]['food_id']=='FDC:1'
        r=client.get('/foods/FDC:1/nutrients',params={'grams':150})
        assert r.json()['nutrients'][0]['portion_amount']==6
        assert r.json()['nutrients'][0]['unit']=='mg'
        assert r.json()['excluded_incompatible_rows']==2
        assert r.json()['clinical_effect_inferred'] is False
        assert client.get('/foods/FDC:1/nutrients',params={'grams':0}).status_code==422
        assert client.get('/foods/unknown/nutrients',params={'grams':100}).status_code==404


def test_no_missing_data_as_zero_and_root_boundary(tmp_path):
    with TestClient(create_app(tmp_path)) as client:
        assert client.get('/foods').status_code==503
    with pytest.raises(ValueError,match='inside'):
        FoodRegistry(tmp_path,tmp_path.parent/'outside.duckdb')


def test_impossible_source_mass_is_never_served(tmp_path):
    path=tmp_path/'food.duckdb';make_database(path)
    with duckdb.connect(str(path)) as c:
        c.execute("INSERT INTO food_composition VALUES ('FDC:1','4','Caffeine',233333,'mg','source100g','D','USDA')")
    response=FoodRegistry(tmp_path,path).composition('FDC:1',150)
    assert response['excluded_incompatible_rows']==3
    assert len(response['nutrients'])==1
