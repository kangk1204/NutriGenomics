import json
import sqlite3
import sys
from pathlib import Path
from fastapi.testclient import TestClient
from nutriomics_final.api import create_app

# Install the sibling evidence package for deployments; source checkout works for this integration test.
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'evidence/src'))


def test_literature_snapshot_api_preserves_unreviewed_status_and_pagination(tmp_path):
    db=tmp_path/'literature.sqlite'
    with sqlite3.connect(db) as con:
        con.executescript('CREATE TABLE metadata(key,value_json); CREATE TABLE paper(pmid,payload_json); CREATE TABLE assertion(id,pmid,status,payload_json);')
        con.execute('INSERT INTO metadata VALUES(?,?)',('snapshot',json.dumps({'source_report_sha256':'sha','papers':1,'approved_assertions':0})))
        con.execute('INSERT INTO paper VALUES(?,?)',('123',json.dumps({'pmid':'123','summary_current':True,'summary':{'results':['DASH reduced BP']}})))
        con.execute('INSERT INTO assertion VALUES(?,?,?,?)',('a1','123','candidate',json.dumps({'tier':'extracted','binding_label':False})))
    with TestClient(create_app(tmp_path,literature_database=db)) as client:
        assert client.get('/literature/status').json()['approved_assertions']==0
        assert client.get('/literature/papers/123').json()['pmid']=='123'
        assert client.get('/literature/papers/missing').status_code==404
        assert client.get('/literature/papers?limit=0').status_code==422
        assert client.get('/literature/papers?offset=1').json()['items']==[]
        assert client.get('/literature/evidence?reviewed=true').json()['total']==0
        assert client.get('/literature/evidence?pmid=123').json()['items'][0]['tier']=='extracted'
        assert client.get('/literature/evidence?q=%25').json()['total']==0


def test_missing_or_outside_snapshot_returns_unavailable(tmp_path):
    with TestClient(create_app(tmp_path,literature_database=tmp_path.parent/'outside.sqlite')) as client:
        assert client.get('/literature/status').status_code==503
