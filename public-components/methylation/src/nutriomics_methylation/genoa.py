"""Frozen source88 classifiers on independently deposited GENOA BP proxies.

No target-cohort fit, normalization, batch correction or hyperparameter tuning.
GEO measured-BP proxy is a different outcome from diagnosed essential HTN.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from .evaluation import performance, plots
from .io import sha256, versions, write_json


def bp_proxy(sbp, dbp):
    sbp, dbp = np.asarray(sbp, dtype=float), np.asarray(dbp, dtype=float)
    if sbp.shape != dbp.shape:
        raise ValueError("SBP/DBP measurements are not aligned")
    valid = np.isfinite(sbp) & np.isfinite(dbp) & (sbp > 0) & (dbp > 0)
    y = np.full(sbp.shape, -1, dtype=np.int8)
    y[valid & ((sbp >= 140) | (dbp >= 90))] = 1
    y[valid & (sbp < 120) & (dbp < 80)] = 0
    return y


def signal_to_beta(methylated, unmethylated, detection_p):
    m, u, p = np.asarray(methylated, dtype=float), np.asarray(unmethylated, dtype=float), np.asarray(detection_p, dtype=float)
    if m.shape != u.shape or m.shape != p.shape:
        raise ValueError("M/U/P shapes differ")
    valid = np.isfinite(m) & np.isfinite(u) & np.isfinite(p) & (m >= 0) & (u >= 0) & (p >= 0) & (p <= .01)
    beta = np.full(m.shape, np.nan, dtype=np.float32)
    beta[valid] = m[valid] / (m[valid] + u[valid] + 100.0)
    return beta


def parse_number(value, detection=False):
    value = value.strip('"')
    if detection and value.startswith("<"):
        # A published upper bound below threshold is conservative measurement QC.
        value = value[1:]
    try:
        return float(value)
    except ValueError:
        return np.nan


def metadata_frame(path):
    frame = pd.read_csv(path, keep_default_na=False)
    required = {"sample_id", "raw_sample_id", "SBP", "DBP", "age", "sex"}
    if not required.issubset(frame):
        raise ValueError(f"GENOA metadata requires exact columns {sorted(required)}")
    if frame.sample_id.duplicated().any() or frame.raw_sample_id.duplicated().any() or frame.sample_id.eq("").any() or frame.raw_sample_id.eq("").any():
        raise ValueError("GENOA individual/sample identities must be complete and unique")
    for column in ("SBP", "DBP", "age"):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["y_bp_proxy"] = bp_proxy(frame.SBP, frame.DBP)
    frame["endpoint"] = "GEO measured-BP proxy; medication correction unknown; not diagnosed essential HTN"
    frame["age_50_65_secondary"] = frame.age.between(50, 65, inclusive="both")
    return frame


def validate_metadata_qc(frame, source_qc):
    records = source_qc.get("samples", [])
    if len(records) != len(frame):
        raise ValueError("GENOA CSV does not match the independently checked sample count")
    reference = {record["gsm"]: record for record in records}
    if set(reference) != set(frame.sample_id):
        raise ValueError("GENOA GSM identities differ from acquisition QC")
    for row in frame.itertuples(index=False):
        record = reference[row.sample_id]
        if row.raw_sample_id != record["description"] or row.sex != record["sex"]:
            raise ValueError("GENOA GSM/native sample/sex alignment differs from acquisition QC")
        if not np.allclose([row.SBP, row.DBP, row.age], [record["sbp_mmhg"], record["dbp_mmhg"], record["age_years"]], rtol=0, atol=1e-8, equal_nan=True):
            raise ValueError("GENOA individual phenotype differs from acquisition QC")


def stream_beta(raw_signals, metadata, probes, allow_structural_missing=False):
    """Read native signal rows; retain fixed source probes without any target fit."""
    probes = np.asarray(probes)
    if len(set(probes)) != len(probes):
        raise ValueError("Frozen source probe IDs duplicated")
    positions = {probe: i for i, probe in enumerate(probes)}
    beta = np.full((len(metadata), len(probes)), np.nan, dtype=np.float32)
    seen = set()
    scanned = 0
    with gzip.open(raw_signals, "rt", encoding="utf-8") as stream:
        header = [value.strip('"') for value in next(stream).split()]
        if header[0] != "ID_REF" or len(set(header)) != len(header):
            raise ValueError("Native signal header ID_REF/column uniqueness invalid")
        columns = {column: i for i, column in enumerate(header)}
        native_samples = {column.removesuffix(".Methylated.signal") for column in header[1:] if column.endswith(".Methylated.signal")}
        if native_samples != set(metadata.raw_sample_id):
            raise ValueError("Native raw signal triplets do not exactly match metadata sample identities")
        ordered = []
        for sample in metadata.raw_sample_id:
            names = [sample + "." + suffix for suffix in ("Methylated.signal", "Unmethylated.signal", "Detection.Pval")]
            if not all(name in columns for name in names):
                raise ValueError(f"Incomplete M/U/P triplet for {sample}")
            ordered.append([columns[name] for name in names])
        ordered = np.asarray(ordered)
        if len(header) != 1 + 3 * len(metadata):
            raise ValueError("Raw signals have unknown orphan columns")
        for line in stream:
            scanned += 1
            probe = line.split(maxsplit=1)[0].strip('"')
            if probe not in positions:
                continue
            if probe in seen:
                raise ValueError(f"Native source probe duplicated: {probe}")
            tokens = line.split()
            if len(tokens) != len(header):
                raise ValueError(f"Native row width mismatch for {probe}")
            m = np.array([parse_number(tokens[index]) for index in ordered[:, 0]])
            u = np.array([parse_number(tokens[index]) for index in ordered[:, 1]])
            p = np.array([parse_number(tokens[index], detection=True) for index in ordered[:, 2]])
            beta[:, positions[probe]] = signal_to_beta(m, u, p)
            seen.add(probe)
            if len(seen) % 2000 == 0:
                print(f"GENOA frozen probes: {len(seen)}/{len(probes)}", flush=True)
        # Full stream consumption also checks gzip CRC/truncation.
    missing = sorted(set(probes) - seen)
    if missing and not allow_structural_missing:
        raise ValueError(f"GENOA misses {len(missing)} fixed source measurement IDs; no platform fabrication: {missing[:10]}")
    return beta, {"raw_probe_rows_scanned": scanned, "fixed_source_probes_found": len(seen),
                  "unprovided_source_probes": [str(p) for p in missing], "unprovided_source_probe_count": len(missing),
                  "coverage_policy": "Full frozen source feature width; unprovided target probes remain NaN and use original source-fitted medians only" if allow_structural_missing else "require full source probe coverage",
                  "raw_metadata_sample_identity_join": "exact and bijective", "gzip_full_stream_read": True}


def predict_frozen_source(bundle, target_beta, target_probes):
    if bundle.get("feature_space") != "common" or not np.array_equal(bundle["probes"], target_probes):
        raise ValueError("Target probe identity/order differs from frozen common source model")
    if len(bundle["training_samples"]) != 88 or bundle["model"]["probes"].n_fit_samples_ != 88 or bundle["model"]["scale"].n_samples_seen_ != 88:
        raise ValueError("GENOA evaluation requires the original source88-fitted preprocessing")
    before = joblib.hash(bundle["model"])
    with threadpool_limits(limits=1):
        prediction = bundle["model"].predict_proba(target_beta)[:, 1]
    if before != joblib.hash(bundle["model"]):
        raise AssertionError("Prediction changed source-fitted model state")
    if not np.isfinite(prediction).all():
        raise ValueError("GENOA model probabilities nonfinite")
    return prediction


def prepared_measurements(root, raw_signals, metadata_csv, sample_qc_json):
    """Fixed, nonlearned M/U/P transformation may precede model availability."""
    root, raw_signals, metadata_csv, sample_qc_json = map(Path, (root, raw_signals, metadata_csv, sample_qc_json))
    frame = metadata_frame(metadata_csv)
    source_qc = json.loads(sample_qc_json.read_text(encoding="utf-8"))
    validate_metadata_qc(frame, source_qc)
    probes = np.load(root / "data" / "prepared" / "matrices.npz", allow_pickle=False)["common_probes"]
    cache = root / "data" / "genoa_fixed_common"
    cache.mkdir(parents=True, exist_ok=True)
    contract = {"raw_signals_sha256": sha256(raw_signals), "metadata_csv_sha256": sha256(metadata_csv),
                "sample_qc_json_sha256": sha256(sample_qc_json),
                "frozen_common_probe_sha256": sha256(root / "data" / "prepared" / "common_probes.tsv"),
                "formula": "M/(M+U+100)", "detection_p_threshold": .01, "target_fit": False,
                "structural_missing_policy": "unprovided source probes remain NaN; original source-trained imputer only"}
    manifest = cache / "measurement_cache.json"
    matrix = cache / "target_beta.npz"
    if manifest.exists():
        saved = json.loads(manifest.read_text(encoding="utf-8"))
        if saved["contract"] != contract or saved["matrix_sha256"] != sha256(matrix):
            raise ValueError("GENOA measurement cache input or bytes changed")
        data = np.load(matrix, allow_pickle=False)
        if not np.array_equal(data["sample_ids"], frame.sample_id.to_numpy()) or not np.array_equal(data["probe_ids"], probes):
            raise ValueError("GENOA measurement cache sample/probe identities changed")
        return frame, data["beta"], probes, saved["raw_qc"]
    beta, raw_qc = stream_beta(raw_signals, frame, probes, allow_structural_missing=True)
    np.savez_compressed(matrix, beta=beta, sample_ids=frame.sample_id.to_numpy(dtype=str), probe_ids=probes)
    write_json(manifest, {"contract": contract, "raw_qc": raw_qc, "matrix_sha256": sha256(matrix),
                          "shape": list(beta.shape), "target_learned_preprocessing": False})
    return frame, beta, probes, raw_qc


def run(root, raw_signals, metadata_csv, sample_qc_json, bootstraps=1000):
    root, raw_signals, metadata_csv, sample_qc_json = map(Path, (root, raw_signals, metadata_csv, sample_qc_json))
    original_config_path = root / "results" / "training_config.json"
    original = json.loads(original_config_path.read_text(encoding="utf-8"))
    source_qc = json.loads(sample_qc_json.read_text(encoding="utf-8"))
    if source_qc.get("data_accession") != "GSE157131" or source_qc.get("platform") != "GPL13534" or not source_qc.get("exact_sample_id_join"):
        raise ValueError("Actual independently verified GENOA450k sample mapping QC required")
    if source_qc.get("family_id_present") or source_qc.get("medication_present") or source_qc.get("clinical_htn_label_present"):
        raise ValueError("GENOA proxy protocol metadata contract changed; record a new protocol rather than silently reinterpret")
    frame = metadata_frame(metadata_csv)
    validate_metadata_qc(frame, source_qc)
    probes = np.load(root / "data" / "prepared" / "matrices.npz", allow_pickle=False)["common_probes"]
    if sha256(root / "data" / "prepared" / "common_probes.tsv") != original["common_mask_sha256"]:
        raise ValueError("Frozen source common probe mask changed")
    output = root / "results_genoa" / "GSE157131_bp_proxy"
    output.mkdir(parents=True, exist_ok=True)
    config = {"accession": "GSE157131", "platform": "GPL13534", "source_accession": "GSE193795", "source_n": 88,
              "target_endpoint": "Measured-BP proxy; not diagnosed essential hypertension", "target_n_metadata": len(frame),
              "case_definition": "SBP>=140 OR DBP>=90", "control_definition": "SBP<120 AND DBP<80",
              "middle_bp_excluded": True, "secondary_age_band": "50<=age<=65, predefined before target predictions",
              "signal_beta_formula": "M/(M+U+100)", "per_measurement_qc": "Detection P>.01 or invalid M/U/P -> missing",
              "protocol_amendment_before_target_predictions": "Source10749 features retained. Target-file-unprovided probes stay NaN and use existing source88-fitted medians only. No feature dropping or target median fitting.",
              "target_used_for_fit_or_tuning": False, "target_normalization_or_batch_fit": False,
              "age_or_sex_adjusted": False, "family_ids_available": False, "individual_medication_available": False,
              "blood_pressure_medication_correction": "GEO individual values' correction status unknown",
              "threshold": .5, "bootstrap_replicates": bootstraps, "versions": versions(),
              "raw_signals_sha256": sha256(raw_signals), "metadata_csv_sha256": sha256(metadata_csv),
              "sample_qc_json_sha256": sha256(sample_qc_json), "source_training_config_sha256": sha256(original_config_path),
              "common_probe_mask_sha256": original["common_mask_sha256"], "frozen_model_sha256": {},
              "limitations": ["Unknown family dependence is not modeled by participant bootstrap intervals.",
                              "Treated hypertension may appear in measured normal-BP controls; clinical diagnoses/individual medications are unavailable.",
                              "Raw signal beta differs from source deposited processed beta; no cross-cohort normalization is fitted.",
                              "Blood leukocyte fractions/preprocessing/population differences remain; age/sex are descriptive only."]}
    bundles = {}
    for family in ("ridge", "elasticnet"):
        path = root / "results" / f"common_{family}" / "final_model.joblib"
        config["frozen_model_sha256"][family] = sha256(path)
        bundles[family] = joblib.load(path)
        if bundles[family]["common_mask_sha256"] != original["common_mask_sha256"] or set(bundles[family]["training_samples"]) & set(frame.sample_id):
            raise ValueError("GENOA source model mask/sample overlap contract invalid")
    lock_path = output / "analysis_config.json"
    if lock_path.exists() and json.loads(lock_path.read_text(encoding="utf-8")) != config:
        raise ValueError("GENOA evaluation input/model/config changed; preserve fixed outputs and declare another analysis")
    write_json(lock_path, config)
    if (output / "complete.json").exists():
        return json.loads((output / "metrics_all.json").read_text(encoding="utf-8"))
    cached_frame, beta, cached_probes, raw_qc = prepared_measurements(root, raw_signals, metadata_csv, sample_qc_json)
    if not np.array_equal(cached_frame.sample_id, frame.sample_id) or not np.array_equal(cached_probes, probes):
        raise ValueError("Cached GENOA measurement/sample ordering differs from model input")
    frame["fixed_probe_missing_fraction"] = np.mean(~np.isfinite(beta), axis=1)
    measured = ~np.isin(probes, raw_qc["unprovided_source_probes"])
    frame["unprovided_probe_fraction"] = 1 - float(np.mean(measured))
    frame["measured_probe_detection_or_signal_failure_fraction"] = np.mean(~np.isfinite(beta[:, measured]), axis=1)
    # No target-derived sample cutoff; report every declared proxy participant.
    frame.to_csv(output / "target_sample_qc.tsv", sep="\t", index=False)
    qc = {**raw_qc, "sample_missing_fraction_quantiles": {str(q): float(frame.fixed_probe_missing_fraction.quantile(q)) for q in (0, .25, .5, .75, 1)},
          "measured_probe_detection_or_signal_failure_fraction_quantiles": {str(q): float(frame.measured_probe_detection_or_signal_failure_fraction.quantile(q)) for q in (0, .25, .5, .75, 1)},
          "frozen_selected_probe_coverage": {family: {"selected_probes": len(bundle["model"]["probes"].selected_indices_),
               "unprovided_selected_probes": [str(p) for p in bundle["probes"][bundle["model"]["probes"].selected_indices_] if p in set(raw_qc["unprovided_source_probes"])]}
               for family, bundle in bundles.items()},
          "beta_mean_descriptive": float(np.nanmean(beta)), "target_missingness_sample_exclusions": 0,
          "bp_counts": {"case": int((frame.y_bp_proxy == 1).sum()), "control": int((frame.y_bp_proxy == 0).sum()),
                        "middle_or_missing": int((frame.y_bp_proxy == -1).sum())},
          "age_quantiles": {str(q): float(frame.age.quantile(q)) for q in (0, .25, .5, .75, 1)},
          "sex_counts": frame.sex.value_counts().to_dict(), "metadata_mapping": source_qc}
    write_json(output / "measurement_qc.json", qc)
    results = {}
    for family in ("ridge", "elasticnet"):
        prediction = predict_frozen_source(bundles[family], beta, probes)
        predicted = frame.copy()
        predicted["probability"] = prediction
        predicted.to_csv(output / f"{family}_all_predictions.tsv", sep="\t", index=False)
        for subset_name, include in [("all_proxy", frame.y_bp_proxy >= 0),
                                     ("age_50_65_secondary", (frame.y_bp_proxy >= 0) & frame.age_50_65_secondary)]:
            subset = predicted.loc[include].rename(columns={"y_bp_proxy": "y"}).copy()
            subset.to_csv(output / f"{family}_{subset_name}_predictions.tsv", sep="\t", index=False)
            if set(subset.y) != {0, 1}:
                raise ValueError("Declared GENOA proxy evaluation requires both outcomes")
            values = performance(subset.y.to_numpy(), subset.probability.to_numpy(), bootstraps=bootstraps)
            values["endpoint"] = config["target_endpoint"]
            values["family_dependence_captured"] = False
            values["interpretation"] = "Conditional fixed-classifier association with independently deposited measured-BP proxy; not clinical HTN validation or proof of improved disease generalization."
            results[f"{family}/{subset_name}"] = values
            plots(subset, output / f"{family}_{subset_name}_roc_calibration.png")
            (output / "calibration_bins.tsv").replace(output / f"{family}_{subset_name}_calibration_bins.tsv")
    write_json(output / "metrics_all.json", results)
    write_json(output / "complete.json", {"status": "real_data_completed", "proxy_evaluation_n": qc["bp_counts"]["case"] + qc["bp_counts"]["control"],
                                          "clinical_htn_validation": False, "source_models_refitted": False})
    target = root / "docs" / "genoa"
    target.mkdir(parents=True, exist_ok=True)
    import shutil
    for p in output.glob("*"):
        if p.suffix in (".json", ".tsv", ".svg") and p.is_file():
            shutil.copy2(p, target / p.name)
    lines = ["# Independently deposited GENOA measured-BP proxy evaluation", "",
             "Source88 common-probe ridge and elastic-net models remain frozen. Target endpoint is measured-BP proxy, not clinical essential hypertension.",
             "No target-cohort fitting/tuning, covariate adjustment, threshold selection or prediction inversion. Unknown family dependence and medication correction remain.",
             f"Coverage-limited deployment: {raw_qc['fixed_source_probes_found']}/{len(probes)} source probe IDs measured; {raw_qc['unprovided_source_probe_count']} unprovided probes stay NaN and use source-fitted medians. No complete platform-coverage claim.", "",
             "| Analysis | N | AUROC (95% conditional interval) | Brier | Sensitivity | Specificity |",
             "|---|---:|---|---:|---:|---:|"]
    for name, value in results.items():
        lo, hi = value["ci95"]["auroc"]
        lines.append(f"| {name} | {value['n']} | {value['auroc']:.3f} ({lo:.3f}–{hi:.3f}) | {value['brier']:.3f} | {value['sensitivity']:.3f} | {value['specificity']:.3f} |")
    lines += ["", "Case SBP≥140 or DBP≥90; control SBP<120 and DBP<80; middle BP excluded. Age50–65 is a prespecified secondary descriptive population, not matching/covariate adjustment.",
              "Participant-bootstrap intervals condition on fixed predictions and do not account for unknown family dependence or model-selection variability. All original clinical-HTN and reciprocal-transfer results must remain reported separately.", ""]
    (target / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")
    write_json(target / "artifact_manifest.json", {"files": [{"file": p.name, "sha256": sha256(p)} for p in target.iterdir() if p.is_file() and p.name != "artifact_manifest.json"],
               "code_sha256": sha256(Path(__file__)), "versions": versions(), "clinical_endpoint_equivalent": False,
               "source_model_refit": False, "original_results_overwritten": False})
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--raw-signals", type=Path, required=True)
    parser.add_argument("--metadata-csv", type=Path, required=True)
    parser.add_argument("--sample-qc-json", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    args = parser.parse_args(argv)
    run(args.root.resolve(), args.raw_signals.resolve(), args.metadata_csv.resolve(), args.sample_qc_json.resolve(), args.bootstrap_replicates)


if __name__ == "__main__":
    main()
