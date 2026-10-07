from nutriomics_evidence.supervised_relations import candidates, targets, export, tuple_scores, LABELS
import numpy as np


def test_candidate_preserves_source_ids_and_offsets():
    doc = {"pmid": "1", "text": "A inhibits B.", "entities": [
        {"id": "T1", "type": "CHEMICAL", "start": 0, "end": 1, "text": "A"},
        {"id": "T2", "type": "GENE", "start": 11, "end": 12, "text": "B"}]}
    row, = candidates(doc)
    assert row["evidence_text"] == doc["text"][row["evidence_start"]:row["evidence_end"]]
    assert "SELECTED_CHEMICAL" in row["features"] and "SELECTED_GENE" in row["features"]
    assert (row["arg1"], row["arg2"]) == ("T1", "T2")


def test_multilabel_and_outside_candidate_gold_stay_in_denominator():
    row = {"pmid": "1", "arg1": "T1", "arg2": "T2", "evidence_start": 0, "evidence_end": 1, "evidence_text": "x"}
    gold = {("1", "INHIBITOR", "T1", "T2"), ("2", "ACTIVATOR", "T3", "T4")}
    y = targets([row], gold)
    assert y.sum() == 1
    probability = np.zeros((1, len(LABELS)))
    probability[0, LABELS.index("INHIBITOR")] = .9
    predicted = export([row], probability, .5)
    assert tuple_scores(gold, predicted) == {"tp": 1, "fp": 0, "fn": 1, "precision": 1., "recall": .5, "f1": 2/3}
    assert predicted[0]["binding_affinity"] is None
