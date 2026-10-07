"""Synthetic fixtures test artifact contracts, never scientific performance."""
import hashlib
import pandas as pd
import pytest
from nutriomics_dti.data import structure
from nutriomics_dti.io import sha256,write_json
from nutriomics_dti.models import train,evaluate,predict
from nutriomics_dti.splits import load_locked


def fixture(smiles,offset):
    records = []
    sequence = "ACDEFGHIKLMNPQRSTVWY"*2
    target = hashlib.sha256(sequence.encode()).hexdigest()[:24]
    for index,text in enumerate(smiles):
        canonical,compound,scaffold = structure(text)
        records.append({"smiles":canonical,"protein_sequence":sequence,"compound_id":compound,"target_id":target,
                        "pair_id":compound+":"+target,"measurement_id":str(index+offset),"source_group":"synthetic:"+str(index+offset),
                        "pkd":5+(index+offset)*0.05,"scaffold":scaffold})
    return pd.DataFrame(records)


def test_locked_training_evaluation_prediction_and_tamper_detection(tmp_path):
    directory = tmp_path/"split"
    directory.mkdir()
    frames = {"train":fixture(["CCO","CCN","CC","CCC","CCCC","c1ccccc1","c1ccncc1","CC(=O)O"],0),
              "validation":fixture(["CCCl","CCBr"],20),"test":fixture(["CCCO","CCCN"],40)}
    files = {}
    for name,frame in frames.items():
        frame.to_csv(directory/f"{name}.csv.gz",index=False)
        files[name] = {"rows":len(frame),"sha256":sha256(directory/f"{name}.csv.gz")}
    write_json(directory/"lock.json",{"mode":"random","seed":17,"files":files,"overlap_counts":{"pair_id":[0,0,0],"source_group":[0,0,0]}})
    record = train(directory,tmp_path/"model",trees=3)
    assert record["locked_test_sha256"]==files["test"]["sha256"]
    result = evaluate(directory,tmp_path/"model",tmp_path/"result")
    assert result["test_metrics"]["n"]==2
    output = predict(tmp_path/"model",frames["test"][["smiles","protein_sequence"]])
    assert len(output)==2
    assert output.evidence_status.str.contains("not an experimentally measured interaction").all()
    assert (output.predicted_kd_nm>0).all()
    with open(directory/"test.csv.gz","ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError,match="Locked test file changed"):
        load_locked(directory)
