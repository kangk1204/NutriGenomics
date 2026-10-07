import pandas as pd
import pytest
from nutriomics_dti.data import structure
from nutriomics_dti.candidate_cards import DISEASES, cas_checksum, disease_protocol, map_food_chemicals, build_cards


def inputs():
    smiles, key, _ = structure("CCO")
    panel = pd.DataFrame([{"compound_id": key, "smiles": smiles, "pubchem_cid": 702,
                           "food_verified": True, "name": "test-analyte", "positive_food_records": 1}])
    registry = pd.DataFrame([{"ChemicalName": "test-chemical", "ChemicalID": "C1", "PubChemCID": "702",
                              "CasRN": "64-17-5", "InChIKey": ""}])
    diseases = pd.DataFrame([{"DiseaseID": item["disease_id"], "DiseaseName": item["expected_name"]}
                             for item in DISEASES.values()])
    associations = pd.DataFrame([{"ChemicalID": "C1", "DiseaseID": "MESH:D065626",
                                   "DirectEvidence": "marker/mechanism", "PubMedIDs": "123"}])
    return panel, registry, diseases, associations


def test_cas_checksum_and_ambiguous_mapping_fail_closed():
    assert cas_checksum("64-17-5") and not cas_checksum("64-17-6")
    panel, registry, _, _ = inputs()
    duplicated = pd.concat([registry, registry.assign(ChemicalID="C2")], ignore_index=True)
    mappings, rejected = map_food_chemicals(panel, duplicated, {})
    assert not mappings and len(rejected) == 1


def test_fullkey_conflict_never_degrades_to_name_cas_or_cid():
    panel, registry, _, _ = inputs()
    registry["InChIKey"] = structure("CCCO")[1]
    mappings, rejected = map_food_chemicals(panel, registry, {"702": ["64-17-5"]})
    assert not mappings and len(rejected) == 1


def test_marker_association_stays_nonbenefit_with_actual_shortfall():
    result = build_cards(*inputs(), {}, {})
    card = result["cards"][0]
    assert card["unique_therapeutic_pmids"] == 0
    assert not card["efficacy_verified"] and not card["full_stereochemical_identity_verified"]
    assert result["summaries"]["nafld"]["shortfall"] == 16
    assert result["summaries"]["immune_related_inflammation"]["shortfall"] == 6
    assert card["evidence"][0]["species"] == "not reported in this CTD chemical-disease download"


def test_exact_disease_id_name_mismatch_rejected():
    _, _, diseases, _ = inputs()
    diseases.loc[0, "DiseaseName"] = "Fatty Liver, Alcoholic"
    with pytest.raises(ValueError, match="Exact disease"):
        disease_protocol(diseases)
