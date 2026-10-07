"""Exploratory reciprocal transfer; frozen primary analyses are never overwritten."""
from __future__ import annotations

import json
import shutil
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .evaluation import match_groups, performance, plots
from .io import sha256, versions, write_json
from .model import FoldLocalProbeSelector, candidates, classifier, marker_stability


def grouped_splits(y, groups, n_splits, seed):
    y, groups = np.asarray(y), np.asarray(groups)
    if len(y) != len(groups) or len(np.unique(groups)) < n_splits:
        raise ValueError("Too few aligned participant groups for source CV")
    splits = list(GroupKFold(n_splits=n_splits, shuffle=True, random_state=seed).split(np.zeros(len(y)), y, groups))
    for train, validation in splits:
        if set(groups[train]) & set(groups[validation]):
            raise AssertionError("Matched participant groups leaked across source folds")
        if set(y[train]) != {0, 1} or set(y[validation]) != {0, 1}:
            raise ValueError("Every source training/validation fold requires both outcomes")
    return splits


def source_only_tune(X, y, samples, groups, family, seed, folds=4,
                     ks=(100, 500), Cs=(0.01, 0.1, 1.0), max_missing=0.1):
    """The interface has no target features/outcomes: all fitting uses source only."""
    grid = candidates(family, ks, Cs)
    scores, convergence = [[] for _ in grid], [0 for _ in grid]
    split_audit = []
    for inner_id, (train, validation) in enumerate(grouped_splits(y, groups, folds, seed)):
        selector = FoldLocalProbeSelector(max(ks), max_missing).fit(X[train])
        a, b = selector.transform(X[train]), selector.transform(X[validation])
        scaled = {}
        for k in ks:
            scaler = StandardScaler().fit(a[:, :min(k, a.shape[1])])
            scaled[k] = (scaler.transform(a[:, :min(k, a.shape[1])]), scaler.transform(b[:, :min(k, b.shape[1])]))
        for index, params in enumerate(grid):
            train_scaled, validation_scaled = scaled[params["k"]]
            fitted = classifier(family, params["C"], params["l1_ratio"], seed)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", ConvergenceWarning)
                fitted.fit(train_scaled, y[train])
            scores[index].append(float(roc_auc_score(y[validation], fitted.predict_proba(validation_scaled)[:, 1])))
            convergence[index] += sum(isinstance(w.message, ConvergenceWarning) for w in caught)
        split_audit.append({"inner_fold": inner_id, "training_sample_ids": samples[train].tolist(),
                            "validation_sample_ids": samples[validation].tolist(),
                            "training_pair_ids": sorted(set(groups[train])), "validation_pair_ids": sorted(set(groups[validation])),
                            "preprocessing_fit_n": selector.n_fit_samples_})
    means = [float(np.mean(s)) for s in scores]
    best_index = int(np.argmax(means))
    best = grid[best_index]
    fitted = Pipeline([("probes", FoldLocalProbeSelector(best["k"], max_missing)), ("scale", StandardScaler()),
                       ("model", classifier(family, best["C"], best["l1_ratio"], seed))])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        fitted.fit(X, y)
    audit = {"source_sample_ids": samples.tolist(), "best_params": best, "inner_mean_auc": means[best_index],
             "inner_splits": split_audit, "target_used_in_tuning": False,
             "candidate_results": [{**g, "fold_auc": s, "mean_auc": mean, "convergence_warnings": count}
                                   for g, s, mean, count in zip(grid, scores, means, convergence)],
             "refit_convergence_warnings": sum(isinstance(w.message, ConvergenceWarning) for w in caught)}
    return fitted, audit


def one_source_fold(X, y, samples, groups, train, test, family, seed, fold_id):
    with threadpool_limits(limits=1):
        fitted, audit = source_only_tune(X[train], y[train], samples[train], groups[train], family, seed + fold_id, folds=3)
        prediction = fitted.predict_proba(X[test])[:, 1]
    audit.update({"outer_fold": fold_id, "train_sample_ids": samples[train].tolist(), "test_sample_ids": samples[test].tolist(),
                  "selected_indices": fitted["probes"].selected_indices_.tolist(), "coefficients": fitted["model"].coef_[0].tolist()})
    return test, prediction, audit


def reciprocal(root, cores=12, seed=20261001, bootstraps=1000):
    if not 1 <= cores <= 12:
        raise ValueError("Transfer worker CPU limit is 1–12 cores")
    if bootstraps < 100:
        raise ValueError("Transfer evaluation requires at least 100 bootstrap replicates")
    root = Path(root)
    prepared = root / "data" / "prepared"
    original_config_path = root / "results" / "training_config.json"
    if not original_config_path.exists():
        original_config_path = root / "docs" / "validation" / "training_config.json"
    original_config = json.loads(original_config_path.read_text(encoding="utf-8"))
    matrix_path = prepared / "matrices.npz"
    mask_path = prepared / "common_probes.tsv"
    if sha256(matrix_path) != original_config["prepared_matrix_sha256"] or sha256(mask_path) != original_config["common_mask_sha256"]:
        raise ValueError("Original frozen measurement inputs/mask changed; reciprocal transfer refused")
    data = np.load(matrix_path, allow_pickle=False)
    source_X, source_y, source_samples = data["external"], data["external_y"], data["external_samples"]
    include = data["primary_phenotypes"] != "preHT"
    target_samples = data["primary_samples"][include]
    target_y = (data["primary_phenotypes"][include] == "HTN").astype(int)
    probes = data["common_probes"]
    positions = {p: i for i, p in enumerate(data["probe_ids"])}
    indices = np.array([positions[p] for p in probes])
    target_X = data["primary"][include][:, indices]
    if (len(source_y), int(np.sum(source_y)), len(target_y), int(np.sum(target_y))) != (16, 8, 88, 44):
        raise ValueError("Reciprocal protocol requires GSE42774 8+8 and GSE193795 44+44 only")
    if set(source_samples) & set(target_samples) or source_X.shape[1] != len(probes) or target_X.shape[1] != len(probes):
        raise ValueError("Reciprocal sample or frozen common-probe alignment invalid")
    metadata = pd.read_csv(prepared / "external_metadata.tsv", sep="\t", index_col=0)
    if metadata.index.tolist() != source_samples.tolist():
        raise ValueError("Reciprocal source matched-pair metadata misaligned")
    groups = match_groups(metadata)
    if groups is None:
        raise ValueError("Published matched source pairs unavailable; do not invent unpaired CV")
    output = root / "results_transfer" / "42774_to193795"
    output.mkdir(parents=True, exist_ok=True)
    config = {"analysis": "exploratory_reciprocal_transfer", "source": "GSE42774", "source_n": 16,
              "target": "GSE193795", "target_n": 88, "target_previously_observed": True,
              "source_pair_groups": 8, "source_outer_group_folds": 4, "source_inner_group_folds": 3,
              "final_tuning_source_group_folds": 4, "k_grid": [100, 500], "C_grid": [0.01, 0.1, 1.0],
              "elasticnet_l1_ratios": [0.1, 0.5], "max_training_missing_fraction": 0.1, "decision_threshold": 0.5,
              "seed": seed, "cores": cores, "common_probe_count": len(probes), "preHT_used_for_training": False,
              "common_mask_sha256": sha256(mask_path), "prepared_matrix_sha256": sha256(matrix_path),
              "primary_training_config_sha256": sha256(original_config_path),
              "primary_training_config_path": original_config_path.relative_to(root).as_posix(),
              "bootstrap_replicates": bootstraps, "versions": versions(),
              "interpretation": "Supplementary direction-of-transfer experiment. Previously evaluated target outcomes preclude a new untouched external validation claim; no improvement claim from selecting the better direction."}
    existing = output / "transfer_config.json"
    if existing.exists() and json.loads(existing.read_text(encoding="utf-8")) != config:
        raise ValueError("Transfer config changed; preserve this run and use a separately declared analysis")
    write_json(existing, config)
    metrics = {}
    for family in ("ridge", "elasticnet"):
        run = output / family
        run.mkdir(parents=True, exist_ok=True)
        if (run / "complete.json").exists():
            saved = json.loads((run / "metrics.json").read_text(encoding="utf-8"))
            metrics.update(saved)
            continue
        source_oof = np.full(len(source_y), np.nan)
        outer_fold = np.full(len(source_y), -1)
        splits = grouped_splits(source_y, groups, 4, seed)
        results = Parallel(n_jobs=min(cores, 4))(
            delayed(one_source_fold)(source_X, source_y, source_samples, groups, train, test, family, seed, fold)
            for fold, (train, test) in enumerate(splits))
        audits = []
        for test, prediction, audit in results:
            if np.isfinite(source_oof[test]).any():
                raise AssertionError("Source participant assigned to multiple OOF folds")
            source_oof[test], outer_fold[test] = prediction, audit["outer_fold"]
            audits.append(audit)
        if not np.isfinite(source_oof).all():
            raise AssertionError("Missing source OOF predictions")
        pd.DataFrame({"sample_id": source_samples, "y": source_y, "probability": source_oof,
                      "outer_fold": outer_fold, "match_pair": groups}).to_csv(run / "source_oof_predictions.tsv", sep="\t", index=False)
        write_json(run / "source_fold_audit.json", audits)
        marker_stability(audits, probes).to_csv(run / "source_marker_stability.tsv", sep="\t", index=False)
        with threadpool_limits(limits=1):
            fitted, final_audit = source_only_tune(source_X, source_y, source_samples, groups, family, seed)
            # Target outcome is used only after this source-only fit and prediction.
            prediction = fitted.predict_proba(target_X)[:, 1]
        joblib.dump({"model": fitted, "training_samples": source_samples, "probes": probes,
                     "common_mask_sha256": config["common_mask_sha256"]}, run / "source_final_model.joblib")
        write_json(run / "source_final_tuning_audit.json", final_audit)
        target_frame = pd.DataFrame({"sample_id": target_samples, "y": target_y, "probability": prediction,
                                     "cohort": "GSE193795", "previously_evaluated_target": True})
        target_frame.to_csv(run / "target_predictions.tsv", sep="\t", index=False)
        family_metrics = {f"{family}/source_nested4_pair_cv": performance(source_y, source_oof, seed, bootstraps, groups),
                          f"{family}/target_reciprocal_transfer": performance(target_y, prediction, seed, bootstraps)}
        for name, values in family_metrics.items():
            values["interpretation"] = config["interpretation"]
            values["target_used_for_tuning"] = False
        write_json(run / "metrics.json", family_metrics)
        plots(target_frame, run / "target_roc_calibration.png")
        metrics.update(family_metrics)
        write_json(run / "complete.json", {"status": "real_data_completed", "source_n": 16, "target_n": 88,
                                            "target_previously_observed": True, "clinical_generalization_established": False})
    write_json(output / "metrics_all.json", metrics)
    target = root / "docs" / "transfer"
    target.mkdir(parents=True, exist_ok=True)
    for path in [existing, output / "metrics_all.json"]:
        shutil.copy2(path, target / path.name)
    for path in output.glob("*/*.tsv"):
        table = pd.read_csv(path, sep="\t")
        if "marker_stability" in path.name:
            table = table.head(100)
        table.to_csv(target / f"{path.parent.name}_{path.name}", sep="\t", index=False)
    for path in output.glob("*/*audit.json"):
        shutil.copy2(path, target / f"{path.parent.name}_{path.name}")
    lines = ["# Exploratory reciprocal transfer", "", config["interpretation"], "",
             "Source: GSE42774 n16 (8 published matched pairs). Target: previously evaluated GSE193795 HTN44/control44; preHT44 excluded.",
             "All preprocessing and hyperparameters are learned from source training groups only. No target normalization, batch correction, feature selection, threshold selection or coefficient flipping.", "",
             "| Analysis | N | AUROC (95% conditional interval) | Brier | Sensitivity | Specificity |",
             "|---|---:|---|---:|---:|---:|"]
    for name, values in metrics.items():
        lo, hi = values["ci95"]["auroc"]
        lines.append(f"| {name} | {values['n']} | {values['auroc']:.3f} ({lo:.3f}–{hi:.3f}) | {values['brier']:.3f} | {values['sensitivity']:.3f} | {values['specificity']:.3f} |")
    lines += ["", "Intervals condition on fixed realized predictions. Tiny source size and cross-platform/population differences remain. No individual age/sex adjustment was invented.",
              "Original forward-transfer and primary nested-CV results remain frozen and must be reported beside these supplementary results. A more favorable point estimate does not establish improved generalization.", ""]
    (target / "SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")
    write_json(target / "artifact_manifest.json", {"files": [{"file": p.name, "sha256": sha256(p)} for p in target.iterdir() if p.is_file() and p.name != "artifact_manifest.json"],
               "code_sha256": sha256(Path(__file__)), "target_previously_observed": True,
               "code_files": [{"path": p.relative_to(root).as_posix(), "sha256": sha256(p)} for p in sorted((root / "src").rglob("*.py"))],
               "original_primary_results_overwritten": False, "clinical_generalization_claim": False})
    return metrics
