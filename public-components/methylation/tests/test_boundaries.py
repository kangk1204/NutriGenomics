import gzip

import numpy as np
import pandas as pd
import pytest

from nutriomics_methylation.io import common_probe_mask, covariate_report, phenotype, read_matrix
from nutriomics_methylation.model import FoldLocalProbeSelector, nested_predictions
from nutriomics_methylation.evaluation import match_groups, performance


def test_common_mask_accepts_only_measurement_ids_and_is_order_invariant():
    expected = np.array(["cg1", "cg3"])
    np.testing.assert_array_equal(common_probe_mask(["cg3", "cg2", "cg1"], ["cg4", "cg1", "cg3"]), expected)
    with pytest.raises(ValueError, match="Duplicate"):
        common_probe_mask(["cg1", "cg1"], ["cg1"])


def test_test_values_cannot_change_fitted_imputer_or_feature_selection():
    training = np.array([[0.1, 0.2, np.nan, 0.5], [0.3, 0.2, 0.1, 0.5], [0.2, 0.2, 0.3, 0.5], [0.4, 0.2, 0.2, 0.5]])
    fitted = FoldLocalProbeSelector(k=2, max_missing=0.3).fit(training)
    medians = fitted.medians_.copy()
    indices = fitted.selected_indices_.copy()
    fitted.transform(np.array([[999, -100, np.nan, 10]]))
    np.testing.assert_array_equal(fitted.medians_, medians)
    np.testing.assert_array_equal(fitted.selected_indices_, indices)
    assert fitted.n_fit_samples_ == 4
    assert set(indices) == {0, 2}


def test_declared_metadata_labels_not_accession_ranges():
    meta = pd.DataFrame({"characteristic_subject_status": ["healthy control", "hypertensive patient", "prehypertensive patient"]}, index=["GSM3", "GSM1", "GSM2"])
    assert phenotype(meta, "GSE193795").to_dict() == {"GSM3": "control", "GSM1": "HTN", "GSM2": "preHT"}
    with pytest.raises(ValueError, match="Unsupported"):
        phenotype(meta, "invented_external")
    with pytest.raises(ValueError, match="Missing"):
        phenotype(pd.DataFrame({"title": ["HTN"]}), "GSE193795")
    report = covariate_report(meta)
    assert report["individual_age_field"] == []
    assert report["individual_sex_field"] == []
    assert not report["age_adjusted"] and not report["sex_adjusted"]


def test_matrix_label_alignment_checked(tmp_path):
    path = tmp_path / "wrong.gz"
    with gzip.open(path, "wt") as f:
        f.write('!Sample_geo_accession\t"GSM1"\t"GSM2"\n!series_matrix_table_begin\n"ID_REF"\t"GSM2"\t"GSM1"\n"cg1"\t0.1\t0.2\n!series_matrix_table_end\n')
    with pytest.raises(ValueError, match="order"):
        read_matrix(path)


def test_nested_audit_inner_and_outer_samples_never_overlap():
    rng = np.random.default_rng(31)
    X = rng.uniform(0, 1, (30, 20))
    y = np.tile([0, 1], 15)
    samples = np.array([f"S{i}" for i in range(30)])
    p, folds, audit = nested_predictions(X, y, samples, "ridge", "nested5", cores=1, inner_folds=3, ks=(5,), Cs=(0.1,))
    assert len(p) == 30 and np.isfinite(p).all()
    assert len(np.unique(folds)) == 5
    for entry in audit:
        assert not set(entry["train_sample_ids"]).intersection(entry["test_sample_ids"])
        for inner in entry["inner_splits"]:
            assert not set(inner["train_indices"]).intersection(inner["validation_indices"])
            assert inner["preprocessing_fit_n"] == len(inner["train_indices"])


def test_metric_probability_alignment_and_intervals():
    y = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    p = np.array([.1, .2, .4, .7, .3, .6, .8, .9])
    result = performance(y, p, bootstraps=30)
    assert result["sensitivity"] == result["specificity"] == .75
    assert result["ci95"]["sensitivity"][0] < .75 < result["ci95"]["sensitivity"][1]
    assert "does not refit CV" in result["interval_method"]["scope"]
    with pytest.raises(ValueError, match="aligned"):
        performance(y, p[:-1], bootstraps=3)


def test_external_match_pairs_verified():
    metadata = pd.DataFrame({"title": [f"Sample {i}" for i in range(1, 17)],
                             "characteristic_match": [f"match with sample {i+1 if i%2 else i-1}" for i in range(1, 17)]})
    groups = match_groups(metadata)
    assert len(np.unique(groups)) == 8
    metadata.loc[0, "characteristic_match"] = "unknown"
    with pytest.raises(ValueError, match="Cannot map"):
        match_groups(metadata)
