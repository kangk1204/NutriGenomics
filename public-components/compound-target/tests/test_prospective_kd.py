import hashlib
import json
from pathlib import Path
import pandas as pd
import pytest
from nutriomics_dti.data import structure, to_pkd
from nutriomics_dti.source_ledger_v2 import canonical_aliases, close_aliases, native_guards, food_external_guard
from nutriomics_dti.prospective_kd import freeze, development, consume_test


def row(index=0,**changes):
    smi,key,scaffold = structure("C"*(index%12+2)+"O")
    sequence = "ACDEFGHIKLMNPQRSTVWY"*2+"A"*index
    tid = hashlib.sha256(sequence.encode()).hexdigest()[:24]
    value = float(index+1)
    return {"measurement_id":f"m{index}","pair_id":key+":"+tid,"compound_id":key,"target_id":tid,"smiles":smi,"scaffold":scaffold,"protein_sequence":sequence,"source_group":f"doi:10.0/{index}","doi":f"10.0/{index}","pmid":"","patent":"","assay_id":"","endpoint":"Kd","relation":"=","unit":"nM","value_nm":value,"pkd":to_pkd(value),**changes}


def test_alias_normalization_connects_doi_pmid_and_multiple_assays():
    aliases = canonical_aliases(row(doi="https://doi.org/10.0/A",assay_id="one;two",patent="us-123 a",pubchem_aid="aid1963802"))
    assert {"doi:10.0/a","assay:one","assay:two","patent:US123A","pubchem_aid:1963802"}<=set(aliases)
    frame = close_aliases(pd.DataFrame([row(doi="10.0/A",pmid="123"),row(1,doi="",pmid="123")]))
    assert frame.source_group.nunique()==1


@pytest.mark.parametrize("changes",[{"endpoint":"Ki"},{"relation":"<"},{"value_nm":0},{"unit":"uM"},{"pkd":-4},{"compound_id":"name-match"},{"target_id":"wrong"}])
def test_native_kd_guards_reject_label_and_identity_errors(changes):
    with pytest.raises(ValueError):
        native_guards(pd.DataFrame([row(**changes)]))


def test_food_guard_never_accepts_incomplete_historical_ledger():
    reference = pd.DataFrame([row()])
    external = pd.DataFrame([row(1)])
    retained,result = food_external_guard(reference,external,set(external.compound_id),False)
    assert retained.empty and result["status"]=="unevaluated"
    assert result["minimum"]=={"n":30,"compounds":5,"targets":3,"source_groups":3}


def test_training_loader_does_not_read_test_and_failed_test_consumes_attempt(tmp_path):
    frame = pd.DataFrame([row(i) for i in range(200)])
    input_ = tmp_path/"prepared.csv.gz"; frame.to_csv(input_,index=False)
    out = tmp_path/"locks"
    freeze(input_,out,modes=("source",))
    directory = out/"source"/"20261003"
    (directory/"test.csv.gz").write_bytes(b"broken test")
    (training,validation),_ = development(directory)
    assert len(training)>0 and len(validation)>0
    with pytest.raises(ValueError,match="Sealed test"):
        consume_test(directory,"ridge",tmp_path/"result")
    with pytest.raises(FileExistsError):
        consume_test(directory,"ridge",tmp_path/"second")
    assert json.loads((directory/"test_attempt_ridge.json").read_text())["status"]=="attempt_consumed_before_label_access"


def test_baseline_training_reads_only_development(tmp_path):
    from nutriomics_dti.native_kd_runner import train
    frame = pd.DataFrame([row(i) for i in range(150)])
    path = tmp_path/"prepared.csv.gz"; frame.to_csv(path,index=False)
    freeze(path,tmp_path/"locks",modes=("source",))
    directory = tmp_path/"locks/source/20261003"
    (directory/"test.csv.gz").unlink()
    result = train(directory,tmp_path/"ridge","ridge",threads=1)
    assert result["validation_metrics"]["n"]==22
    assert result["fresh_weights"]
    assert not (directory/"test_attempt_ridge.json").exists()
