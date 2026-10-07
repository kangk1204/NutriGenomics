import pytest
from nutriomics_dti.data import to_pkd,structure,prepare,uncertain_construct


def test_pkd_units_and_exact_endpoint():
    assert to_pkd("10",unit="nM")==pytest.approx(8)
    assert to_pkd("= 10",unit="nM")==pytest.approx(8)
    assert to_pkd("0.01",unit="uM")==pytest.approx(8)
    assert to_pkd("1e-8",unit="M")==pytest.approx(8)
    for value in (">10","<5","~2","0","-1","nan","inf"):
        with pytest.raises(ValueError):
            to_pkd(value)
    with pytest.raises(ValueError):
        to_pkd("10",endpoint="IC50")
    with pytest.raises(ValueError):
        to_pkd("10",relation=">")
    with pytest.raises(ValueError):
        to_pkd("10",unit="mg/L")


def test_stereochemical_identity_is_retained():
    assert structure("N[C@@H](C)C(=O)O")[1]!=structure("N[C@H](C)C(=O)O")[1]
    with pytest.raises(ValueError):
        structure("invalid!")


def test_native_parser_filters_censoring_and_endpoint(tmp_path):
    raw = tmp_path/"native.tsv"
    raw.write_text("BindingDB Reactant_set_id\tLigand SMILES\tKd (nM)\tIC50 (nM)\tBindingDB Target Chain Sequence\tArticle DOI\tPMID\tNumber of Protein Chains in Target (>1 implies a multichain complex)\n"
                   "1\tCCO\t10\t\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\t10.1/a\t123\t1\n"
                   "2\tOCC\t10\t\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\t\t123\t1\n"
                   "3\tCCC\t>10\t\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\t10.1/b\t\t1\n"
                   "4\tCCCC\t\t1\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\t10.1/c\t\t1\n",encoding="utf-8")
    audit = prepare(raw,tmp_path/"prepared")
    assert audit["counts"]["retained_measurements"]==1
    assert audit["counts"]["duplicate_measurements"]==1
    assert audit["counts"]["censored_or_invalid_kd"]==1
    assert audit["counts"]["no_kd"]==1


def test_numbered_native_chain_schema_and_construct_guard(tmp_path):
    raw = tmp_path/"native202609.tsv"
    raw.write_text("BindingDB Reactant_set_id\tLigand SMILES\tKd (nM)\tBindingDB Target Chain Sequence 1\tTarget Name\tArticle DOI\tUniProt (SwissProt) Primary ID of Target Chain 1\tNumber of Protein Chains in Target (>1 implies a multichain complex)\tBindingDB Target Chain Sequence 2\n"
                   "1\tCCO\t10\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\tKinase\t10.1/a\tP12345\t1\t\n"
                   "2\tCCC\t10\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\tKinase [T315I]\t10.1/b\tP12345\t1\t\n"
                   "3\tCCN\t10\tACDEFGHIKLMNPQRSTVWYACDEFGHIK\tKinase fusion\t10.1/c\tP12345\t1\t\n",encoding="utf-8")
    audit = prepare(raw,tmp_path/"prepared")
    assert audit["counts"]["retained_measurements"]==1
    assert audit["counts"]["explicit_mutant_truncation_or_fusion"]==2
    import pandas as pd
    data = pd.read_csv(tmp_path/"prepared"/"measurements.csv.gz")
    assert data.uniprot.tolist()==["P12345"]
    assert uncertain_construct("ABL1 T315I")
    assert uncertain_construct("protein [1-200]")
    assert not uncertain_construct("Cyclin-dependent kinase 2")
    assert not uncertain_construct("Adenosine receptor A2A")
