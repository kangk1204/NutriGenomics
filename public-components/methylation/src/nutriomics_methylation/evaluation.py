from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logit
from scipy.stats import beta
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (auc, average_precision_score, brier_score_loss,
                             confusion_matrix, precision_recall_curve, roc_auc_score, roc_curve)

from .io import sha256, versions, write_json


def point_metrics(y, p, threshold=0.5):
    y, p = np.asarray(y), np.asarray(p, dtype=float)
    if (y.ndim != 1 or p.ndim != 1 or y.shape != p.shape or y.size == 0
            or y.dtype.kind not in 'biuf' or not np.isfinite(y).all()
            or set(y) != {0, 1} or not np.isfinite(p).all()
            or ((p < 0) | (p > 1)).any()
            or not isinstance(threshold, (int, float)) or isinstance(threshold, bool)
            or not np.isfinite(threshold) or not 0 <= threshold <= 1):
        raise ValueError("Metrics require aligned finite binary outcome probabilities")
    y = y.astype(int)
    tn, fp, fn, tp = confusion_matrix(y, p >= threshold, labels=[0, 1]).ravel()
    precision, recall, _ = precision_recall_curve(y, p)
    return {"n": len(y), "n_positive": int(y.sum()), "auroc": float(roc_auc_score(y, p)),
            "average_precision": float(average_precision_score(y, p)), "auprc_trapezoid": float(auc(recall, precision)),
            "brier": float(brier_score_loss(y, p)), "sensitivity": float(tp / (tp + fn)),
            "specificity": float(tn / (tn + fp)), "threshold": threshold,
            "tp": int(tp), "tn": int(tn), "fp": int(fp), "fn": int(fn)}


def exact_binomial_interval(successes, trials, alpha=0.05):
    return [0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, trials - successes + 1)),
            1.0 if successes == trials else float(beta.ppf(1 - alpha / 2, successes + 1, trials - successes))]


def performance(y, p, seed=20261001, bootstraps=1000, groups=None):
    y, p = np.asarray(y), np.asarray(p)
    output = point_metrics(y, p)
    rng = np.random.default_rng(seed)
    draws = []
    group_ids = np.unique(groups) if groups is not None else None
    for _ in range(bootstraps):
        if groups is None:
            indices = np.concatenate([rng.choice(np.flatnonzero(y == c), int((y == c).sum()), replace=True) for c in (0, 1)])
        else:
            indices = np.concatenate([np.flatnonzero(groups == g) for g in rng.choice(group_ids, len(group_ids), replace=True)])
        if len(np.unique(y[indices])) == 2:
            draws.append(point_metrics(y[indices], p[indices]))
    output["ci95"] = {key: np.quantile([d[key] for d in draws], [0.025, 0.975]).tolist()
                      for key in ["auroc", "average_precision", "auprc_trapezoid", "brier"]}
    output["ci95"]["sensitivity"] = exact_binomial_interval(output["tp"], output["tp"] + output["fn"])
    output["ci95"]["specificity"] = exact_binomial_interval(output["tn"], output["tn"] + output["fp"])
    output["interval_method"] = {"discrimination_brier": "matched-pair bootstrap" if groups is not None else "stratified participant bootstrap",
                                  "sensitivity_specificity": "Clopper-Pearson exact binomial; external pairing not modeled for these intervals",
                                  "bootstrap_replicates": len(draws), "seed": seed,
                                  "scope": "Conditional on realized predictions, does not refit CV and does not capture model-selection variability."}
    z = logit(np.clip(p, 1e-6, 1 - 1e-6)).reshape(-1, 1)
    calibration = LogisticRegression(C=1e6, solver="lbfgs", max_iter=3000).fit(z, y)
    output["calibration_slope_descriptive"] = float(calibration.coef_[0, 0])
    output["calibration_intercept_descriptive"] = float(calibration.intercept_[0])
    return output


def match_groups(metadata):
    column = "characteristic_match"
    if column not in metadata:
        return None
    import re
    group = []
    for _, row in metadata.iterrows():
        own = re.search(r"Sample\s*(\d+)", row["title"], re.I)
        partner = re.search(r"sample\s*(\d+)", row[column], re.I)
        if own is None or partner is None:
            raise ValueError("Cannot map published matched sample IDs")
        group.append("-".join(str(i) for i in sorted([int(own.group(1)), int(partner.group(1))])))
    counts = pd.Series(group).value_counts()
    if len(counts) != 8 or not (counts == 2).all():
        raise ValueError("External matched pairs are inconsistent")
    return np.array(group)


def plots(predictions, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    y, p = predictions.y.to_numpy(), predictions.probability.to_numpy()
    fpr, tpr, _ = roc_curve(y, p)
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), constrained_layout=True)
    axes[0].plot(fpr, tpr, color="#156082", label=f"AUROC {roc_auc_score(y, p):.3f}")
    axes[0].plot([0, 1], [0, 1], linestyle="--", color="grey")
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", xlim=(0, 1), ylim=(0, 1))
    axes[0].legend(loc="lower right")
    bins = pd.DataFrame({"p": p, "y": y, "bin": np.minimum((p * 5).astype(int), 4)}).groupby("bin").agg(mean_probability=("p", "mean"), event_fraction=("y", "mean"), n=("y", "size"))
    axes[1].plot(bins.mean_probability, bins.event_fraction, "o-", color="#E97132")
    axes[1].plot([0, 1], [0, 1], linestyle="--", color="grey")
    axes[1].set(xlabel="Mean predicted probability", ylabel="Observed HTN fraction", xlim=(0, 1), ylim=(0, 1))
    bins.to_csv(path.parent / "calibration_bins.tsv", sep="\t")
    fig.savefig(path, dpi=180)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def evaluate(root, bootstraps=1000):
    rows, all_metrics = [], {}
    config = json.loads((root / "results" / "training_config.json").read_text(encoding="utf-8"))
    required = []
    for space in config["spaces"]:
        for family in ("ridge", "elasticnet"):
            directory = root / "results" / f"{space}_{family}"
            required.extend(directory / method / "oof_predictions.tsv" for method in config["methods"])
            required.append(directory / "preHT_predictions.tsv")
            if space == "common":
                required.append(directory / "external_predictions.tsv")
    missing = [str(p.relative_to(root)) for p in required if not p.exists()]
    if missing:
        raise ValueError(f"Training contract is incomplete; evaluation/export cannot present partial runs as complete: {missing}")
    emeta = pd.read_csv(root / "data" / "prepared" / "external_metadata.tsv", sep="\t", index_col=0)
    pairs = match_groups(emeta)
    for path in sorted((root / "results").glob("**/*predictions.tsv")):
        frame = pd.read_csv(path, sep="\t")
        key = path.parent.relative_to(root / "results").as_posix() + ("/external" if path.name == "external_predictions.tsv" else "")
        if path.name == "preHT_predictions.tsv":
            write_json(path.parent / "preHT_summary.json", {"n": len(frame), "used_for_training": False,
                       "probability_median": float(frame.probability.median()), "probability_iqr": frame.probability.quantile([0.25, 0.75]).tolist(),
                       "fraction_at_or_above_prespecified_threshold": float((frame.probability >= 0.5).mean()),
                       "interpretation": "No binary disease outcome assigned to preHT; probabilities are exploratory, not calibrated risk."})
            continue
        aligned_pairs = None
        if path.name == "external_predictions.tsv":
            if frame.sample_id.tolist() != emeta.index.tolist():
                raise ValueError("External labels/pairs do not align with prediction sample order")
            aligned_pairs = pairs
        metrics = performance(frame.y.to_numpy(), frame.probability.to_numpy(), bootstraps=bootstraps, groups=aligned_pairs)
        metrics["source"] = path.relative_to(root).as_posix()
        metrics["prediction_sha256"] = sha256(path)
        if aligned_pairs is not None:
            metrics["interpretation"] = "Exploratory independent cross-platform transfer, n16 young African American males; cannot establish clinical generalization."
        write_json(path.parent / ("external_metrics.json" if aligned_pairs is not None else "metrics.json"), metrics)
        plots(frame, path.parent / ("external_roc_calibration.png" if aligned_pairs is not None else "roc_calibration.png"))
        for measure in ["auroc", "average_precision", "auprc_trapezoid", "brier", "sensitivity", "specificity"]:
            rows.append({"analysis": key, "n": metrics["n"], "metric": measure, "estimate": metrics[measure],
                         "ci95_lower": metrics["ci95"][measure][0], "ci95_upper": metrics["ci95"][measure][1],
                         "threshold": metrics["threshold"], "interval_scope": "conditional on realized predictions"})
        all_metrics[key] = metrics
    if not rows:
        raise ValueError("No completed real predictions to evaluate")
    pd.DataFrame(rows).to_csv(root / "results" / "performance.tsv", sep="\t", index=False)
    write_json(root / "results" / "metrics_all.json", all_metrics)
    return all_metrics


def export(root):
    """Small auditable result tables for version control; raw inputs remain ignored."""
    import shutil
    import subprocess
    target = root / "docs" / "validation"
    target.mkdir(parents=True, exist_ok=True)
    written_paths = set()
    for path in [root / "results" / "performance.tsv", root / "results" / "metrics_all.json",
                 root / "results" / "training_config.json", root / "data" / "source_manifest.json",
                 root / "data" / "prepared" / "qc_manifest.json"]:
        if not path.exists():
            raise FileNotFoundError(f"Required real validation artifact missing: {path}")
        destination = target / path.name
        shutil.copy2(path, destination)
        written_paths.add(destination)
    for run in sorted((root / "results").glob("*_*")):
        if not run.is_dir():
            continue
        for path in run.glob("**/marker_stability.tsv"):
            table = pd.read_csv(path, sep="\t", low_memory=False).head(100)
            relative_name = "_".join(path.relative_to(root / "results").parts)
            destination = target / relative_name
            table.to_csv(destination, sep="\t", index=False)
            written_paths.add(destination)
        for path in run.glob("**/*predictions.tsv"):
            destination = target / "_".join(path.relative_to(root / "results").parts)
            shutil.copy2(path, destination)
            written_paths.add(destination)
        for path in run.glob("**/*roc_calibration.svg"):
            destination = target / "_".join(path.relative_to(root / "results").parts)
            shutil.copy2(path, destination)
            written_paths.add(destination)
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, stderr=subprocess.DEVNULL, text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit = "uncommitted_initial_implementation"
    summary = pd.read_csv(target / "performance.tsv", sep="\t")
    lines = ["# Observed validation results", "", "Research-only results on deposited public samples. No guaranteed performance, clinical risk or causal dietary response.",
             "", "All intervals condition on the realized predictions. They do not include CV/model-selection uncertainty.", "",
             "| Analysis | N | AUROC (95% conditional interval) | Average precision | Brier |", "|---|---:|---|---:|---:|"]
    for analysis, group in summary.groupby("analysis", sort=True):
        group = group.set_index("metric")
        a = group.loc["auroc"]
        lines.append(f"| {analysis} | {int(a['n'])} | {a.estimate:.3f} ({a.ci95_lower:.3f}–{a.ci95_upper:.3f}) | {group.loc['average_precision','estimate']:.3f} | {group.loc['brier','estimate']:.3f} |")
    lines += ["", "The common probe mask was frozen from measurement IDs before outcome mapping. All learned QC, median imputation, variance selection, scaling and tuning used training folds only.",
              "", "preHT44 was excluded from fit/tuning. GSE42774 is an exploratory n16 transfer study across 450k/27k, age, ancestry and sex distributions. Individual age/sex adjustment was not fabricated.", ""]
    summary_path = target / "SUMMARY.md"
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    written_paths.add(summary_path)
    files = [{"path": p.relative_to(root).as_posix(), "sha256": sha256(p), "bytes": p.stat().st_size} for p in sorted(written_paths)]
    code = [{"path": p.relative_to(root).as_posix(), "sha256": sha256(p)} for p in sorted((root / "src").rglob("*.py"))]
    write_json(target / "artifact_manifest.json", {"files": files, "inventory_scope": "files_written_by_current_export", "code_files": code, "code_revision_at_export": commit, "versions": versions(),
               "not_clinical": True, "no_causal_dietary_claim": True, "analysis_unit": "participant/GSM", "data_access": "official public GEO HTTPS"})
    return target
