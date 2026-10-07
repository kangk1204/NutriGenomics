import json
import sys
import threading
import time
from pathlib import Path
import pytest
from fastapi.testclient import TestClient
from nutriomics_final.api import create_app, command_for
from nutriomics_final.jobs import JobStore,Runner

def synthetic_manifest_command(job):
    # Tiny API fixture satisfies the same artifact contract as a real dispatcher.
    program = '''from pathlib import Path
import json, hashlib
p=Path(OUTPUT)
v={'protocol_id':PROTOCOL,'input_hash':'a'*64,'model_version':'fixture-v1','sources':[],'output_manifest':'output_manifest.json'}
(p/'provenance.json').write_text(json.dumps(v))
files=[{'path':f.name,'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(p.iterdir()) if f.is_file() and f.name not in {'run.log','output_manifest.json'}]
(p/'output_manifest.json').write_text(json.dumps({**v,'files':files}))
'''
    return [sys.executable, '-c', program.replace('OUTPUT', repr(job['output'])).replace('PROTOCOL', repr(job['protocol_id']))]

def test_protocol_migration_preserves_all_legacy_records_and_result_bytes(tmp_path):
    import sqlite3
    import hashlib
    database=tmp_path/'legacy.sqlite';outputs=tmp_path/'outputs';outputs.mkdir()
    with sqlite3.connect(database) as db:
        db.execute('CREATE TABLE jobs (id TEXT PRIMARY KEY, algorithm TEXT NOT NULL, action TEXT NOT NULL, parameters TEXT NOT NULL, state TEXT NOT NULL, created TEXT NOT NULL, updated TEXT NOT NULL, output TEXT NOT NULL, error TEXT, code_version TEXT NOT NULL, model_version TEXT, source_manifest TEXT, parent_id TEXT)')
        for index in range(13):
            folder=outputs/str(index);folder.mkdir();(folder/'result.txt').write_text('frozen '+str(index))
            db.execute('INSERT INTO jobs VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)',(str(index),'evidence','query','{"limit":10}','succeeded','old-created','old-updated',str(folder),None,'old-code','old-model','[]',None))
        before=db.execute('SELECT * FROM jobs ORDER BY id').fetchall()
    result_hashes={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in outputs.rglob('result.txt')}
    store=JobStore(database,outputs)
    with sqlite3.connect(database) as db:
        after=db.execute('SELECT id,algorithm,action,parameters,state,created,updated,output,error,code_version,model_version,source_manifest,parent_id FROM jobs ORDER BY id').fetchall()
        assert before==after
        assert db.execute('SELECT count(*) FROM jobs WHERE protocol_id IS NULL AND input_hash IS NULL AND output_manifest IS NULL').fetchone()[0]==13
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==value for p,value in result_hashes.items())
    assert len(store.list())==13
    assert store.create('evidence','query',{},'new-code',protocol_id='phase1-final-20261003-v1')['protocol_id']=='phase1-final-20261003-v1'

def test_success_result_hash_and_provenance(tmp_path):
    store=JobStore(tmp_path/'jobs.sqlite',tmp_path/'outputs')
    job=store.create('evidence','query',{},'abc123')
    def command(job):
        return [sys.executable,'-c',"from pathlib import Path;import json; p=Path("+repr(job['output'])+");(p/'value.txt').write_text('measured');(p/'provenance.json').write_text(json.dumps({'model_version':'v1','code_version':'source-bytes-v2','sources':[{'id':'CTD','kind':'direct'}]}))"]
    assert Runner(store,command).run(job['id'])['state']=='succeeded'
    result=store.result(job['id'])
    assert result['job']['model_version']=='v1'
    assert result['job']['code_version']=='source-bytes-v2'
    assert result['job']['source_manifest'][0]['id']=='CTD'
    assert any(x['path']=='value.txt' and len(x['sha256'])==64 for x in result['files'])

def test_failure_cancel_and_explicit_restart(tmp_path):
    store=JobStore(tmp_path/'jobs.sqlite',tmp_path/'outputs')
    failed=store.create('evidence','validate',{},'v1')
    assert Runner(store,lambda j:[sys.executable,'-c','raise SystemExit(3)']).run(failed['id'])['state']=='failed'
    job=store.create('evidence','validate',{},'v1')
    runner=Runner(store,lambda j:[sys.executable,'-c','import time;time.sleep(30)'])
    thread=threading.Thread(target=lambda:runner.run(job['id']));thread.start()
    for _ in range(100):
        if job['id'] in runner.processes:break
        time.sleep(.01)
    assert runner.cancel(job['id'])['state']=='cancelled'
    thread.join(3); assert not thread.is_alive()
    recovering=store.create('evidence','validate',{},'v1');store.transition(recovering['id'],'queued','running')
    assert store.recover()==1
    assert store.get(recovering['id'])['state']=='interrupted'
    with pytest.raises(ValueError):store.result(job['id'])

def test_path_escape_and_shell_text_never_executed(tmp_path):
    outside=tmp_path.parent/'outside.sqlite';outside.write_text('')
    with pytest.raises(ValueError):
        command_for(tmp_path,sys.executable,{'algorithm':'evidence','action':'validate','parameters':{'database':'../outside.sqlite'},'output':str(tmp_path/'result')})
    db=tmp_path/'data.sqlite';db.write_text('')
    command=command_for(tmp_path,sys.executable,{'algorithm':'evidence','action':'query','parameters':{'database':'data.sqlite','disease':'MESH:D006973; touch /tmp/pwned'},'output':str(tmp_path/'result')})
    assert command[-1]=='MESH:D006973; touch /tmp/pwned'
    assert isinstance(command,list)

def test_api_invalid_jobs_cancel_retry_and_result(tmp_path):
    app=create_app(tmp_path,command_builder=synthetic_manifest_command)
    with TestClient(app) as client:
        assert client.post('/jobs',json={'algorithm':'evidence','action':'train'}).status_code==422
        assert client.get('/jobs/missing').status_code==404
        response=client.post('/jobs',json={'algorithm':'evidence','action':'validate'})
        assert response.status_code==202
        job_id=response.json()['id']
        for _ in range(100):
            state=client.get('/jobs/'+job_id).json()['state']
            if state=='succeeded':break
            time.sleep(.02)
        assert state=='succeeded'
        assert client.get('/jobs/'+job_id+'/result').status_code==200
        assert client.get('/jobs/'+job_id+'/files/run.log').status_code==200
        assert client.get('/jobs/'+job_id+'/files/unlisted.txt').status_code==404
        assert client.post('/jobs/'+job_id+'/retry').status_code==409

def test_pending_includes_older_jobs_beyond_recent_listing(tmp_path):
    store=JobStore(tmp_path/'jobs.sqlite',tmp_path/'outputs')
    first=store.create('evidence','validate',{},'v1')
    for _ in range(101):
        store.create('evidence','validate',{},'v1')
    assert first['id'] not in {item['id'] for item in store.list()}
    assert store.pending()[0]['id']==first['id']

def test_transient_metadata_lock_does_not_kill_worker(tmp_path):
    import sqlite3
    app=create_app(tmp_path,command_builder=synthetic_manifest_command)
    original=app.state.store.pending
    calls=[]
    def transient():
        calls.append(1)
        if len(calls)==1:raise sqlite3.OperationalError('database is locked')
        return original()
    app.state.store.pending=transient
    with TestClient(app) as client:
        job=client.post('/jobs',json={'algorithm':'evidence','action':'validate'}).json()
        for _ in range(100):
            state=client.get('/jobs/'+job['id']).json()['state']
            if state=='succeeded':break
            time.sleep(.02)
        assert state=='succeeded'
        assert client.get('/health').json()['worker_alive']
