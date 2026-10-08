"""Tiny real default-dispatch regressions; no model training or downloads."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import pytest
from nutriomics_final.execution_sources import source_layout,verify_module_origins

CHILD = r"""
import hashlib,json,os,sys,time
from pathlib import Path
root=Path(sys.argv[1]).resolve();legacy=sys.argv[2]=='legacy';hostile=sys.argv[3]=='hostile'
final='nutriomics-final-report' if legacy else 'integration';evidence='nutriomics-evidence-engine' if legacy else 'evidence'
sys.path[:0]=[str(root/final/'src'),str(root/evidence/'src')]
from fastapi.testclient import TestClient
from nutriomics_final.api import create_app
from nutriomics_evidence.graph import connect,node,source,edge
p=root/'tiny.sqlite';db=connect(p)
with db:
    source(db,{'id':'fixture-source','source':'synthetic','sha256':'a'*64,'version':'v1'})
    node(db,'fixture:one','chemical','synthetic');node(db,'fixture:two','disease','synthetic')
    edge(db,'fixture:one','fixture:two','fixture_relation','curated',None,{'synthetic':True},'fixture-source','123')
db.close();before=hashlib.sha256(p.read_bytes()).hexdigest()
if hostile:
    shadow=root/'shadow';shadow.mkdir();package=shadow/'nutriomics_evidence';package.mkdir()
    (package/'__init__.py').write_text('')
    (package/'__main__.py').write_text('raise RuntimeError("harmless conflicting module selected")')
    os.environ['PYTHONPATH']=str(shadow);os.chdir(shadow)
app=create_app(root)
with TestClient(app) as client:
    response=client.post('/jobs',json={'algorithm':'evidence','action':'query','parameters':{'database':'tiny.sqlite','limit':1}})
    assert response.status_code==202,response.text
    job=response.json();identifier=job['id']
    for _ in range(200):
        job=client.get('/jobs/'+identifier).json()
        if job['state'] in {'succeeded','failed'}:break
        time.sleep(.02)
    assert job['state']=='succeeded',(job, (Path(job['output'])/'run.log').read_text())
    result=client.get('/jobs/'+identifier+'/result');assert result.status_code==200,result.text
    assert result.json()['artifact_verification']['status']=='frozen_verified'
    data=json.loads(client.get('/jobs/'+identifier+'/files/evidence.json').text)
    assert len(data['evidence'])==1 and data['evidence'][0]['relation']=='fixture_relation'
    output=Path(job['output']);provenance=json.loads((output/'provenance.json').read_text())
    declared={str(Path(item['path']).resolve()):item['sha256'] for snapshot in provenance['code_snapshots'] for item in snapshot['files']}
    assert provenance['module_origins']
    assert any(item['module']=='nutriomics_evidence.graph' for item in provenance['module_origins'])
    assert any(item['module']=='nutriomics_final.dispatch' for item in provenance['module_origins'])
    for item in provenance['module_origins']:
        actual=Path(item['path']).resolve()
        assert actual.is_relative_to(root) and 'shadow' not in actual.parts,item
        assert declared[str(actual)]==item['sha256']==hashlib.sha256(actual.read_bytes()).hexdigest(),item
    manifest=json.loads((output/'output_manifest.json').read_text())
    for item in manifest['files']:
        artifact=output/item['path'];assert item['sha256']==hashlib.sha256(artifact.read_bytes()).hexdigest()
    assert provenance['sources'][0]['sha256']==before==hashlib.sha256(p.read_bytes()).hexdigest()
    assert provenance['protocol_id']=='evidence-query-20261003-v1'
    assert (output/'execution_modules.json').name in {item['path'] for item in manifest['files']}
    (output/'evidence.json').write_text('changed fixture')
    assert client.get('/jobs/'+identifier+'/result').status_code==409
    print(json.dumps({'layout':'legacy' if legacy else 'public','hostile_import_environment':hostile,'status':'passed','origins':len(provenance['module_origins']),'artifacts':len(manifest['files']),'database_unchanged':True,'tampered_result_rejected':True}))
"""

@pytest.mark.parametrize('legacy',[False,True])
@pytest.mark.parametrize('hostile',[False,True])
def test_default_evidence_dispatch_binds_real_code_and_results(tmp_path,legacy,hostile):
    component=Path(__file__).resolve().parents[1]
    for source,target in [(component,'nutriomics-final-report' if legacy else 'integration'),
                          (component.parent/'evidence','nutriomics-evidence-engine' if legacy else 'evidence')]:
        shutil.copytree(source/'src',tmp_path/target/'src',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    process=subprocess.run([sys.executable,'-I','-c',CHILD,str(tmp_path),'legacy' if legacy else 'public','hostile' if hostile else 'plain'],
                           text=True,capture_output=True,timeout=15)
    assert process.returncode==0,process.stdout+process.stderr
    receipt=json.loads(process.stdout)
    assert receipt['status']=='passed' and receipt['origins']>=2


def test_ambiguous_layout_and_missing_layout_fail_closed(tmp_path):
    with pytest.raises(ValueError,match='missing'):source_layout(tmp_path,'evidence')
    for directory,module in [('integration','nutriomics_final'),('evidence','nutriomics_evidence'),
                             ('nutriomics-final-report','nutriomics_final'),('nutriomics-evidence-engine','nutriomics_evidence')]:
        package=tmp_path/directory/'src'/module;package.mkdir(parents=True);(package/'__init__.py').write_text('')
    with pytest.raises(ValueError,match='Ambiguous'):source_layout(tmp_path,'evidence')


def test_module_origin_receipt_must_match_before_snapshot(tmp_path):
    path=tmp_path/'fixture.py';path.write_text('synthetic = True')
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    snapshots=[{'files':[{'path':str(path),'sha256':digest}]}]
    valid=[{'module':'fixture','path':str(path),'sha256':digest}]
    verify_module_origins(snapshots,valid)
    for invalid in [[],[{'module':'fixture','path':str(path),'sha256':'b'*64}],
                    [{'module':'fixture','path':str(tmp_path/'other.py'),'sha256':digest}]]:
        with pytest.raises(ValueError):verify_module_origins(snapshots,invalid)


def test_source_layout_rejects_outside_root(tmp_path):
    outside=tmp_path/'outside';outside.mkdir();root=tmp_path/'root';root.mkdir()
    (root/'nutriomics-final-report').symlink_to(outside,target_is_directory=True)
    with pytest.raises(ValueError,match='escaped'):source_layout(root,'evidence')
