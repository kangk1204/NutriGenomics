import gzip

import numpy as np
import pandas as pd
import pytest
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from nutriomics_methylation.genoa import bp_proxy, predict_frozen_source, signal_to_beta, stream_beta, validate_metadata_qc
from nutriomics_methylation.model import FoldLocalProbeSelector, classifier


def test_bp_proxy_prespecified_boundaries_and_missingness():
    y = bp_proxy([140, 130, 119, 120, 110, np.nan, 0], [70, 90, 79, 79, 80, 70, 70])
    assert y.tolist() == [1, 1, 0, -1, -1, -1, -1]


def test_signal_beta_offset_and_measurement_detection_qc():
    beta = signal_to_beta([100, 100, -1, 100, 100], [100, 100, 100, 100, 100], [.001, .011, .001, np.nan, .01])
    assert beta[0] == pytest.approx(1/3)
    assert beta[4] == pytest.approx(1/3)
    assert np.isnan(beta[1:4]).all()


def test_genoa_native_signal_sample_and_probe_order_not_names_or_position(tmp_path):
    path = tmp_path / "signals.gz"
    with gzip.open(path, "wt") as out:
        out.write('ID_REF B.Methylated.signal B.Unmethylated.signal B.Detection.Pval A.Methylated.signal A.Unmethylated.signal A.Detection.Pval\n')
        out.write('cg2 100 100 <1E-16 200 100 .005\n')
        out.write('cg1 200 100 .011 100 100 .001\n')
    meta = pd.DataFrame({"raw_sample_id": ["A", "B"]})
    beta, qc = stream_beta(path, meta, np.array(["cg1", "cg2"]))
    assert beta.shape == (2, 2)
    assert beta[0, 0] == pytest.approx(1/3) and beta[0, 1] == .5
    assert np.isnan(beta[1, 0]) and beta[1, 1] == pytest.approx(1/3)
    assert qc["raw_metadata_sample_identity_join"] == "exact and bijective"
    with pytest.raises(ValueError, match="misses"):
        stream_beta(path, meta, np.array(["cg_absent"]))
    with pytest.raises(ValueError, match="exactly match"):
        stream_beta(path, pd.DataFrame({"raw_sample_id": ["A", "C"]}), np.array(["cg1"]))
    partial, partial_qc = stream_beta(path, meta, np.array(["cg1", "cg_absent"]), allow_structural_missing=True)
    assert np.isnan(partial[:, 1]).all()
    assert partial_qc["unprovided_source_probes"] == ["cg_absent"]


def test_genoa_prediction_cannot_fit_target_or_change_source_state():
    rng = np.random.default_rng(18)
    X = rng.random((88, 20), dtype=np.float32)
    y = np.tile([0, 1], 44)
    model = Pipeline([("probes", FoldLocalProbeSelector(5)), ("scale", StandardScaler()), ("model", classifier("ridge", .1, None, 17))]).fit(X, y)
    bundle = {"feature_space": "common", "probes": np.array([f"cg{i}" for i in range(20)]),
              "training_samples": np.array([f"source{i}" for i in range(88)]), "model": model}
    medians, mean = model["probes"].medians_.copy(), model["scale"].mean_.copy()
    result = predict_frozen_source(bundle, np.full((272, 20), np.nan), bundle["probes"])
    assert len(result) == 272 and np.isfinite(result).all()
    np.testing.assert_equal(model["probes"].medians_, medians)
    np.testing.assert_equal(model["scale"].mean_, mean)
    assert model["probes"].n_fit_samples_ == model["scale"].n_samples_seen_ == 88
    with pytest.raises(ValueError, match="identity/order"):
        predict_frozen_source(bundle, X, bundle["probes"][::-1])


def test_genoa_phenotype_cannot_move_to_another_native_sample():
    frame = pd.DataFrame({"sample_id": ["GSM1"], "raw_sample_id": ["native1"], "sex": ["F"], "SBP": [140.], "DBP": [70.], "age": [60.]})
    qc = {"samples": [{"gsm": "GSM1", "description": "native1", "sex": "F", "sbp_mmhg": 140., "dbp_mmhg": 70., "age_years": 60.}]}
    validate_metadata_qc(frame, qc)
    with pytest.raises(ValueError, match="phenotype"):
        validate_metadata_qc(frame.assign(SBP=110.), qc)
    with pytest.raises(ValueError, match="alignment"):
        validate_metadata_qc(frame.assign(raw_sample_id="native2"), qc)
