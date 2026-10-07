import json
from nutriomics_evidence.food import ingest_food,meal
from nutriomics_evidence.graph import connect

def test_official_negative_quantity_retained_but_never_summed(tmp_path):
    path=tmp_path/'fdc.json'
    path.write_text(json.dumps({'FoundationFoods':[{'fdcId':1,'description':'test','dataType':'Foundation','foodNutrients':[
      {'nutrient':{'id':1,'name':'Anomalous nutrient','unitName':'mg'},'amount':-1},
      {'nutrient':{'id':2,'name':'Measured nutrient','unitName':'mg'},'amount':2,'foodNutrientDerivation':{'code':'A'}}]},None]}))
    database=tmp_path/'graph.sqlite'
    ingest_food(database,'fdc',[path])
    db=connect(database)
    row=db.execute("SELECT * FROM measurement WHERE nutrient='FDCNutrient:1'").fetchone()
    assert row['amount'] is None and row['status']=='invalid_source_value'
    assert json.loads(row['raw_json'])['amount']==-1
    quality=json.loads(db.execute('SELECT metadata_json FROM source').fetchone()[0])['food_quality_audit']
    assert quality=={'null_food_entries':1,'invalid_quantity_entries':1,'valid_food_records':1}
    result=meal(db,[{'food_id':'FDC:1','grams':100}])
    assert result['totals']==[{'nutrient':'FDCNutrient:2','unit':'mg','amount':2.0}]
    assert result['excluded_or_unknown'][0]['status']=='invalid_source_value'
    db.close()
