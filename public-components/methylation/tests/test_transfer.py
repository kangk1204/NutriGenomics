import numpy as np
import pytest

from nutriomics_methylation.transfer import grouped_splits, source_only_tune


def example_source():
    rng = np.random.default_rng(8)
    X = rng.random((16, 30), dtype=np.float32)
    X[0, 0] = np.nan
    y = np.tile([0, 1], 8)
    samples = np.array([f"source{i}" for i in range(16)])
    groups = np.repeat([f"pair{i}" for i in range(8)], 2)
    return X, y, samples, groups


def test_transfer_matching_partner_never_crosses_train_validation():
    _, y, _, groups = example_source()
    splits = grouped_splits(y, groups, 4, 17)
    seen = []
    for train, validation in splits:
        assert not set(groups[train]) & set(groups[validation])
        assert set(y[train]) == set(y[validation]) == {0, 1}
        seen.extend(validation)
    assert sorted(seen) == list(range(16))


def test_transfer_source_only_fit_cannot_learn_target_scaling_or_imputation():
    X, y, samples, groups = example_source()
    fitted, audit = source_only_tune(X, y, samples, groups, "ridge", 17, ks=(5,), Cs=(.1,))
    medians = fitted["probes"].medians_.copy()
    means = fitted["scale"].mean_.copy()
    fitted.predict_proba(np.full((88, 30), 1000, dtype=np.float32))
    fitted.predict_proba(np.full((88, 30), np.nan, dtype=np.float32))
    np.testing.assert_equal(fitted["probes"].medians_, medians)
    np.testing.assert_equal(fitted["scale"].mean_, means)
    assert fitted["probes"].n_fit_samples_ == 16
    assert fitted["scale"].n_samples_seen_ == 16
    assert audit["target_used_in_tuning"] is False
    for fold in audit["inner_splits"]:
        assert fold["preprocessing_fit_n"] == 12
        assert set(fold["training_sample_ids"]) | set(fold["validation_sample_ids"]) == set(samples)
        assert not set(fold["training_pair_ids"]) & set(fold["validation_pair_ids"])


def test_transfer_rejects_missing_groups_or_single_class_folds():
    _, y, _, groups = example_source()
    with pytest.raises(ValueError, match="Too few aligned"):
        grouped_splits(y, groups[:-1], 4, 17)
    with pytest.raises(ValueError, match="both outcomes"):
        grouped_splits(np.zeros(16), groups, 4, 17)
