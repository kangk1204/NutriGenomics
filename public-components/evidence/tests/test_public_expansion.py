import csv
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import zipfile
import duckdb
import pandas as pd
import pytest
from nutriomics_evidence.public_expansion import merge_ffq_cycle,stage_fdc,structure_check,verified_source,mass_exceeds_basis_sql

def fixture_archive(root, tables):
    (root/'raw').mkdir();(root/'receipts').mkdir();(root/'derived').mkdir()
    path=root/'raw'/'fdc.zip'
    with zipfile.ZipFile(path,'w') as archive:
        for name,rows in tables.items():
            text=io.StringIO();writer=csv.writer(text);writer.writerows(rows)
            archive.writestr('release/'+name,text.getvalue())
    (root/'receipts'/'usda_fdc_202604.json').write_text(json.dumps({'status':'complete','sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'path':str(path)}))

def test_food_version_dedup_composition_and_exact_baseline(tmp_path):
    fixture_archive(tmp_path,{
        'food.csv':[['fdc_id','data_type','description','food_category_id','publication_date'],['1','branded_food','Same product old','','2024-01-01'],['2','branded_food','Same product new','','2026-01-01'],['3','foundation_food','Measured apple','','2026-01-01'],['4','sample_food','Laboratory sample','','2026-01-01']],
        'branded_food.csv':[['fdc_id','gtin_upc','market_country','brand_owner'],['1','123456789012','US','Brand'],['2','00123456789012','US','Brand']],
        'nutrient.csv':[['id','name','unit_name'],['10','Vitamin C','MG']],
        'food_nutrient.csv':[['id','fdc_id','nutrient_id','amount','derivation_id','data_points','min','max','median','loq','footnote','min_year_acquired'],['a','1','10','3','','','','','','','',''],['b','2','10','4','','','','','','','',''],['c','3','10','0','','','','','','','',''],['d','3','10','-1','','','','','','','',''],['e','3','99','5','','','','','','','','']]
    })
    baseline=tmp_path/'baseline.sqlite';db=sqlite3.connect(baseline);db.execute('CREATE TABLE node(id TEXT,kind TEXT,name TEXT,metadata_json TEXT)');db.execute("INSERT INTO node VALUES('FDC:3','food','Apple','{}')");db.commit();db.close()
    con=duckdb.connect();summary={};stage_fdc(con,tmp_path,summary,baseline)
    assert summary['food']['source_records']==4
    assert summary['food']['canonical_product_or_food_identities']==2
    assert summary['food']['composition_records']==2
    assert summary['food']['checks']['invalid_amount_or_unit_native_rows']==2
    assert summary['food']['delta']['exact_added_canonical_food']==1
    assert con.execute("SELECT count(*) FROM food WHERE food_id='FDC:4'").fetchone()[0]==0
    assert con.execute('SELECT food_id,amount,basis FROM food_composition ORDER BY food_id').fetchall()==[('FDC:2',4.0,'source100g'),('FDC:3',0.0,'source100g')]

def test_source_hash_corruption_fails_before_load(tmp_path):
    fixture_archive(tmp_path,{'food.csv':[['fdc_id'],['1']]})
    with (tmp_path/'raw/fdc.zip').open('ab') as output:output.write(b'corruption')
    with pytest.raises(ValueError,match='SHA mismatch'):verified_source(tmp_path,'usda_fdc_202604')

def test_per100g_mass_bounds_keep_exact_limits_and_nonmass_units():
    con=duckdb.connect();con.execute('CREATE TABLE measurements(amount DOUBLE,unit VARCHAR)')
    con.executemany('INSERT INTO measurements VALUES(?,?)',[(100,'G'),(100.1,'g'),(100000,'MG'),(100001,'MG'),(100000000,'UG'),(100000001,'UG'),(13333333,'IU'),(200,'KCAL')])
    assert con.execute('SELECT unit,amount FROM measurements WHERE '+mass_exceeds_basis_sql()).fetchall()==[('g',100.1),('MG',100001.0),('UG',100000001.0)]
    con.close()

def ffq_fixture():
    return {
      'DEMO':pd.DataFrame({'SEQN':[7.,8.],'SDMVPSU':[1.,2.],'SDMVSTRA':[11.,11.],'WTMEC2YR':[40.,60.],'WTINT2YR':[42.,62.],'RIDAGEYR':[50.,60.],'RIAGENDR':[1.,2.],'RIDRETH1':[1.,2.]}),
      'FFQRAW':pd.DataFrame({'SEQN':[7.,8.],'WTS_FFQ':[50.,70.],'FFQ_MISS':[0.,1.]}),
      'BPX':pd.DataFrame({'SEQN':[7.,8.],'BPXSY1':[0.,150.],'BPXSY2':[120.,140.],'BPXDI1':[80.,90.],'BPXDI2':[80.,100.]}),
      'BPQ':pd.DataFrame({'SEQN':[7.,8.],'BPQ020':[1.,2.]}),
      'FFQDC':pd.DataFrame({'SEQN':[7.,8.],'FFQ_VAR':[1.,1.],'FFQ_FOOD':[1.,1.],'FFQ_FREQ':[.5,1.]})
    }

def test_ffq_preserves_survey_design_and_diagnosis_distinction():
    people,long=merge_ffq_cycle(ffq_fixture(),'2003-2004')
    assert people.WTS_FFQ.tolist()==[50.,70.]
    assert people.observed_sbp_mean.tolist()==[120.,145.]
    assert people.observed_bp_ge140_90.tolist()==[False,True]
    assert people.self_reported_htn_diagnosis.tolist()==[True,False]
    assert long.frequency_unit.tolist()==['times_per_day']*2
    assert not long.portion_size_available.any()
    p2,_=merge_ffq_cycle(ffq_fixture(),'2005-2006')
    assert not set(people.participant_id)&set(p2.participant_id)

def test_ffq_duplicate_person_and_absent_special_weights_rejected():
    frames=ffq_fixture();frames['BPX']=pd.concat([frames['BPX'],frames['BPX'].iloc[[0]]])
    with pytest.raises(ValueError,match='duplicate person'):merge_ffq_cycle(frames,'2003-2004')
    frames=ffq_fixture();frames['FFQRAW']=frames['FFQRAW'].drop(columns='WTS_FFQ')
    with pytest.raises(ValueError,match='FFQ survey weights'):merge_ffq_cycle(frames,'2003-2004')

def test_full_inchikey_and_mixture_guards():
    ethanol=('id','CCO','InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3','LFQSCWFLJHTTHZ-UHFFFAOYSA-N','ethanol','foodplant','FooDB','doi')
    checked=structure_check(ethanol)
    assert checked[10] is True and checked[11] is True
    assert checked[12]=='source_asserted_food_database_membership'
    wrong=structure_check(ethanol[:3]+('AAAAAAAAAAAAAA-UHFFFAOYSA-N',)+ethanol[4:])
    assert wrong[10] is False
    mixture=structure_check((ethanol[0],'CCO.O')+ethanol[2:])
    assert mixture[10] is True and mixture[11] is False
