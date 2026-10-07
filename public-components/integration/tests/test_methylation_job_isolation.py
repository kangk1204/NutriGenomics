import json
import sys
from pathlib import Path
import pytest
from nutriomics_final.api import command_for


def test_evaluate_copies_predictions_instead_of_writing_historical_root(tmp_path):
    study=tmp_path/'study';out=tmp_path/'newjob';out.mkdir()
    contents={'results/training_config.json':json.dumps({'spaces':['common'],'methods':['nested5']}),
        'data/prepared/external_metadata.tsv':'id\tlabel\n1\t0\n',
        'results/common_ridge/nested5/oof_predictions.tsv':'sample_id\ty\tprobability\n1\t0\t0.2\n',
        'results/common_ridge/final_training_audit.json':'{}',
        'docs/validation/metrics_all.json':'{"stale":true}'}
    for name,content in contents.items():
        path=study/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(content)
    before={name:(study/name).read_bytes() for name in contents}
    job={'algorithm':'methylation','action':'evaluate','parameters':{'root':'study'},'output':str(out)}
    command=command_for(tmp_path,sys.executable,job)
    assert command[command.index('--root')+1]==str(out)
    assert (out/'results/common_ridge/nested5/oof_predictions.tsv').read_bytes()==before['results/common_ridge/nested5/oof_predictions.tsv']
    assert not (out/'docs/validation/metrics_all.json').exists()
    assert before=={name:(study/name).read_bytes() for name in contents}
    assert len(json.loads((out/'evaluation_inputs.json').read_text()))==4


def test_evaluate_refuses_ignored_scope_overrides_and_dry_run_does_not_copy(tmp_path):
    study=tmp_path/'study'
    for name in ('results/training_config.json','data/prepared/external_metadata.tsv'):
        path=study/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('{}')
    job={'algorithm':'methylation','action':'evaluate','parameters':{'root':'study','spaces':['full']},'output':str(tmp_path/'job')}
    with pytest.raises(ValueError,match='scope is fixed'):
        command_for(tmp_path,sys.executable,job,prepare=False)
    job['parameters']={'root':'study'}
    command_for(tmp_path,sys.executable,job,prepare=False)
    assert not (tmp_path/'job').exists()
