from __future__ import annotations

import itertools
import json
import warnings
from collections import Counter
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from .io import sha256, versions, write_json


class FoldLocalProbeSelector(TransformerMixin, BaseEstimator):
    """Training-only missingness filter, median imputation and variance ordering.

    The returned columns retain variance rank, allowing tuning k without fitting
    preprocessing more than once per inner training fold. No labels are used.
    """
    def __init__(self, k=500, max_missing=0.1):
        self.k = k
        self.max_missing = max_missing

    def fit(self, X, y=None):
        X = np.asarray(X, dtype=np.float32)
        self.n_features_in_ = X.shape[1]
        self.n_fit_samples_ = X.shape[0]
        observed = np.isfinite(X)
        keep = observed.mean(axis=0) >= 1 - self.max_missing
        # Avoid all-missing warnings and preserve original feature positions.
        candidates = np.flatnonzero(keep & observed.any(axis=0))
        if not len(candidates):
            raise ValueError("No probes pass training-fold missingness QC")
        clean = X[:, candidates]
        medians = np.nanmedian(clean, axis=0)
        clean = np.where(np.isfinite(clean), clean, medians)
        variances = np.var(clean, axis=0)
        nonconstant = variances > 0
        candidates, medians, variances = candidates[nonconstant], medians[nonconstant], variances[nonconstant]
        if not len(candidates):
            raise ValueError("No nonconstant probes in training fold")
        # Stable tie ordering depends only on frozen measurement ID order.
        order = np.argsort(-variances, kind="stable")[:min(self.k, len(candidates))]
        self.selected_indices_ = candidates[order]
        self.medians_ = medians[order]
        self.variances_ = variances[order]
        return self

    def transform(self, X):
        X = np.asarray(X, dtype=np.float32)
        if X.shape[1] != self.n_features_in_:
            raise ValueError("Prediction probe width differs from fitted probe mask")
        selected = X[:, self.selected_indices_]
        return np.where(np.isfinite(selected), selected, self.medians_)


def classifier(family, C, l1_ratio, seed):
    if family == "ridge":
        return LogisticRegression(C=C, penalty="l2", solver="liblinear", max_iter=3000, tol=1e-4, random_state=seed)
    if family == "elasticnet":
        return LogisticRegression(C=C, penalty="elasticnet", solver="saga", l1_ratio=l1_ratio,
                                  max_iter=5000, tol=1e-3, random_state=seed)
    raise ValueError(f"Unknown family {family}")


def candidates(family, ks=(100, 500), Cs=(0.01, 0.1, 1.0)):
    ratios = (0.1, 0.5) if family == "elasticnet" else (None,)
    return [{"k": k, "C": C, "l1_ratio": ratio} for k, C, ratio in itertools.product(ks, Cs, ratios)]


def tune(X, y, family, seed, inner_folds=5, ks=(100, 500), Cs=(0.01, 0.1, 1.0), max_missing=0.1):
    grid = candidates(family, ks, Cs)
    scores = [[] for _ in grid]
    convergence = [0 for _ in grid]
    split = StratifiedKFold(n_splits=inner_folds, shuffle=True, random_state=seed)
    split_audit = []
    for inner_id, (train, test) in enumerate(split.split(X, y)):
        assert not set(train).intersection(test)
        selector = FoldLocalProbeSelector(max(ks), max_missing).fit(X[train])
        training = selector.transform(X[train])
        validation = selector.transform(X[test])
        scaled = {}
        for k in ks:
            width = min(k, training.shape[1])
            scaler = StandardScaler().fit(training[:, :width])
            scaled[k] = (scaler.transform(training[:, :width]), scaler.transform(validation[:, :width]))
        for i, params in enumerate(grid):
            a, b = scaled[params["k"]]
            model = classifier(family, params["C"], params["l1_ratio"], seed)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", ConvergenceWarning)
                model.fit(a, y[train])
            convergence[i] += sum(isinstance(w.message, ConvergenceWarning) for w in caught)
            scores[i].append(float(roc_auc_score(y[test], model.predict_proba(b)[:, 1])))
        split_audit.append({"inner_fold": inner_id, "train_indices": train.tolist(), "validation_indices": test.tolist(),
                            "preprocessing_fit_n": selector.n_fit_samples_})
    means = np.array([np.mean(s) for s in scores])
    # Stable tie-break: smaller k, stronger regularization, listed l1 ratio.
    best_index = int(np.argmax(means))
    best = grid[best_index]
    fitted = Pipeline([("probes", FoldLocalProbeSelector(best["k"], max_missing)), ("scale", StandardScaler()),
                       ("model", classifier(family, best["C"], best["l1_ratio"], seed))])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        fitted.fit(X, y)
    audit = {"best_params": best, "inner_mean_auc": float(means[best_index]), "inner_splits": split_audit,
             "candidate_results": [{**g, "mean_auc": float(np.mean(s)), "fold_auc": s, "convergence_warnings": c}
                                   for g, s, c in zip(grid, scores, convergence)],
             "refit_convergence_warnings": sum(isinstance(w.message, ConvergenceWarning) for w in caught)}
    return fitted, audit


def one_outer(X, y, train, test, family, seed, inner_folds, ks, Cs, max_missing, fold_id, samples):
    with threadpool_limits(limits=1):
        fitted, audit = tune(X[train], y[train], family, seed + fold_id, inner_folds, ks, Cs, max_missing)
        predictions = fitted.predict_proba(X[test])[:, 1]
    audit.update({"outer_fold": fold_id, "train_sample_ids": samples[train].tolist(), "test_sample_ids": samples[test].tolist(),
                  "all_preprocessing_training_only": True, "selected_indices": fitted["probes"].selected_indices_.tolist(),
                  "coefficients": fitted["model"].coef_[0].tolist()})
    return test, predictions, audit


def nested_predictions(X, y, samples, family, method, seed=20261001, cores=12, inner_folds=5,
                       ks=(100, 500), Cs=(0.01, 0.1, 1.0), max_missing=0.1, temp_folder=None):
    if cores < 1 or cores > 12:
        raise ValueError("Worker CPU limit is 1–12 cores")
    outer = (StratifiedKFold(5, shuffle=True, random_state=seed) if method == "nested5" else LeaveOneOut())
    if method not in ("nested5", "loocv"):
        raise ValueError("Unknown outer CV method")
    splits = list(outer.split(X, y))
    outputs = Parallel(n_jobs=cores, max_nbytes="10M", temp_folder=temp_folder)(
        delayed(one_outer)(X, y, a, b, family, seed, inner_folds, ks, Cs, max_missing, i, samples)
        for i, (a, b) in enumerate(splits))
    prediction = np.full(len(y), np.nan)
    fold = np.full(len(y), -1)
    audits = []
    for test, values, audit in outputs:
        if np.isfinite(prediction[test]).any():
            raise AssertionError("Repeated test assignment")
        prediction[test], fold[test] = values, audit["outer_fold"]
        audits.append(audit)
    if not np.isfinite(prediction).all():
        raise AssertionError("Not every participant has exactly one OOF prediction")
    return prediction, fold, audits


def marker_stability(audits, probes):
    selected = Counter()
    nonzero = Counter()
    coefficients = {}
    for audit in audits:
        for index, coef in zip(audit["selected_indices"], audit["coefficients"]):
            probe = probes[index]
            selected[probe] += 1
            if abs(coef) > 1e-8:
                nonzero[probe] += 1
            coefficients.setdefault(probe, []).append(coef)
    rows = [{"probe_id": p, "selected_folds": n, "selection_frequency": n / len(audits),
             "nonzero_frequency": nonzero[p] / len(audits), "median_standardized_coefficient_when_selected": float(np.median(coefficients[p])),
             "positive_sign_fraction_when_selected": float(np.mean(np.array(coefficients[p]) > 0)),
             "interpretation": "exploratory stability; not an independently validated biomarker"} for p, n in selected.items()]
    return pd.DataFrame(rows).sort_values(["selection_frequency", "nonzero_frequency"], ascending=False)


def train(root, spaces=("common", "full"), methods=("nested5", "loocv"), cores=12, seed=20261001,
          ks=(100, 500), Cs=(0.01, 0.1, 1.0), max_missing=0.1):
    prepared = root / "data" / "prepared"
    data = np.load(prepared / "matrices.npz", allow_pickle=False)
    phenotype = data["primary_phenotypes"]
    include = phenotype != "preHT"
    y = (phenotype[include] == "HTN").astype(np.int8)
    samples = data["primary_samples"][include]
    probe_ids = data["probe_ids"]
    common = data["common_probes"]
    if len(samples) != 88 or set(y) != {0, 1} or np.sum(y) != 44:
        raise ValueError("Primary model must contain HTN44/control44 only")
    output = root / "results"
    output.mkdir(parents=True, exist_ok=True)
    temporary = root / "cache" / "joblib"
    temporary.mkdir(parents=True, exist_ok=True)
    config = {"spaces": list(spaces), "methods": list(methods), "cores": cores, "seed": seed, "k_grid": list(ks),
              "C_grid": list(Cs), "inner_folds": 5, "max_training_missing_fraction": max_missing, "decision_threshold": 0.5,
              "threshold_selection": "fixed before outcomes, never chosen using external results", "preHT_used_for_training": False,
              "common_mask_sha256": sha256(prepared / "common_probes.tsv"), "versions": versions()}
    config["prepared_matrix_sha256"] = sha256(prepared / "matrices.npz")
    existing = output / "training_config.json"
    if existing.exists() and json.loads(existing.read_text(encoding="utf-8")) != config:
        raise ValueError("Existing results belong to a different config/input/version. Use a separate result directory; do not silently resume stale outputs.")
    write_json(output / "training_config.json", config)
    annotations = pd.read_csv(prepared / "probe_annotation.tsv", sep="\t", index_col=0, low_memory=False)
    for space in spaces:
        if space not in ("full", "common"):
            raise ValueError("Feature space must be full or common")
        indices = np.arange(len(probe_ids)) if space == "full" else np.array([np.flatnonzero(probe_ids == p)[0] for p in common])
        probes = probe_ids[indices]
        X = data["primary"][include][:, indices]
        for family in ("ridge", "elasticnet"):
            run = output / f"{space}_{family}"
            run.mkdir(parents=True, exist_ok=True)
            for method in methods:
                target = run / method
                target.mkdir(parents=True, exist_ok=True)
                if (target / "complete.json").exists():
                    print(f"Skipping complete {space}/{family}/{method}", flush=True)
                    continue
                print(f"Fitting {space}/{family}/{method}: n={len(y)} p={X.shape[1]} cores={cores}", flush=True)
                pred, folds, audit = nested_predictions(X, y, samples, family, method, seed, cores, ks=ks, Cs=Cs, max_missing=max_missing, temp_folder=str(temporary))
                pd.DataFrame({"sample_id": samples, "y": y, "probability": pred, "outer_fold": folds}).to_csv(target / "oof_predictions.tsv", sep="\t", index=False)
                write_json(target / "fold_audit.json", audit)
                marker_stability(audit, probes).join(annotations, on="probe_id", rsuffix="_annotation").to_csv(target / "marker_stability.tsv", sep="\t", index=False)
                write_json(target / "complete.json", {"n": len(y), "p": X.shape[1], "method": method, "family": family, "status": "real_data_completed"})
            if not (run / "final_model.joblib").exists():
                print(f"Final training-only refit {space}/{family}", flush=True)
                with threadpool_limits(limits=1):
                    model, final_audit = tune(X, y, family, seed, ks=ks, Cs=Cs, max_missing=max_missing)
                joblib.dump({"model": model, "probes": probes, "training_samples": samples, "feature_space": space,
                             "common_mask_sha256": config["common_mask_sha256"]}, run / "final_model.joblib")
                write_json(run / "final_training_audit.json", final_audit)
                prediction = model.predict_proba(data["primary"][~include][:, indices])[:, 1]
                pd.DataFrame({"sample_id": data["primary_samples"][~include], "phenotype": "preHT", "probability": prediction,
                              "used_in_fit_or_tuning": False}).to_csv(run / "preHT_predictions.tsv", sep="\t", index=False)
                if space == "common":
                    if not np.array_equal(probes, common):
                        raise AssertionError("External prediction probe order was changed")
                    external_prediction = model.predict_proba(data["external"])[:, 1]
                    pd.DataFrame({"sample_id": data["external_samples"], "y": data["external_y"], "probability": external_prediction}).to_csv(run / "external_predictions.tsv", sep="\t", index=False)
                else:
                    write_json(run / "external_not_applicable.json", {"reason": "Full450k feature space cannot be applied to 27k; no missing-probe fabrication."})
    return config
