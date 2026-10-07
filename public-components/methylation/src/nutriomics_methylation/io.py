from __future__ import annotations

import csv
import gzip
import hashlib
import json
import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

SOURCES = {
    "GSE193795_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193795/matrix/GSE193795_series_matrix.txt.gz",
    "GSE193795_family": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193795/soft/GSE193795_family.soft.gz",
    "GSE42774_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42774/matrix/GSE42774_series_matrix.txt.gz",
    "GSE42774_family": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42774/soft/GSE42774_family.soft.gz",
}


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def versions():
    from importlib.metadata import version
    return {"python": platform.python_version(), **{p: version(p) for p in ["numpy", "pandas", "scipy", "scikit-learn", "joblib"]}}


def fetch(root: Path):
    raw = root / "data" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    manifest = []
    for key, url in SOURCES.items():
        path = raw / url.rsplit("/", 1)[1]
        if not path.exists():
            partial = path.with_suffix(path.suffix + ".partial")
            print(f"Downloading {key}: {url}", flush=True)
            with urlopen(Request(url, headers={"User-Agent": "NutriOmicsResearch/0.1 public GEO reproducibility"}), timeout=600) as response, partial.open("wb") as out:
                for chunk in iter(lambda: response.read(1024 * 1024), b""):
                    out.write(chunk)
            partial.replace(path)
        # Validate complete gzip, not just existence. Detection P extraction happens separately.
        with gzip.open(path, "rb") as f:
            for _ in iter(lambda: f.read(1024 * 1024), b""):
                pass
        manifest.append({"key": key, "url": url, "path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size,
                         "sha256": sha256(path), "checked_at_utc": datetime.now(timezone.utc).isoformat()})
    write_json(root / "data" / "source_manifest.json", {"sources": manifest, "versions": versions()})
    return manifest


def read_matrix(path):
    headers = {}
    rows = []
    ids = []
    columns = None
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as f:
        in_table = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            values = next(csv.reader([line], delimiter="\t"))
            if in_table:
                if columns is None:
                    columns = values[1:]
                else:
                    ids.append(values[0])
                    rows.append(np.array([float(v) if v not in ("null", "NULL", "NA", "", "NaN") else np.nan for v in values[1:]], dtype=np.float32))
            elif line.startswith("!Sample_"):
                headers.setdefault(values[0], []).append(values[1:])
    if not rows or columns is None:
        raise ValueError(f"No matrix values in {path}")
    table = pd.DataFrame(np.stack(rows), index=ids, columns=columns)
    if table.index.has_duplicates or table.columns.has_duplicates:
        raise ValueError("Duplicate CpG or sample accession")
    accessions = headers["!Sample_geo_accession"][0]
    if columns != accessions:
        raise ValueError("Matrix columns differ from sample metadata order")
    metadata = pd.DataFrame(index=accessions)
    metadata.index.name = "sample_id"
    for key, entries in headers.items():
        for j, values in enumerate(entries):
            if len(values) != len(accessions):
                raise ValueError(f"Wrong metadata width: {key}")
            column = key.removeprefix("!Sample_") + (f"_{j+1}" if len(entries) > 1 else "")
            metadata[column] = values
            if key.startswith("!Sample_characteristics"):
                for sample_id, value in zip(accessions, values):
                    if ":" in value:
                        field, val = value.split(":", 1)
                        field = "characteristic_" + re.sub(r"\W+", "_", field.strip().lower())
                        metadata.loc[sample_id, field] = val.strip()
    return table, metadata


def phenotype(metadata, dataset):
    """Map declared phenotypes only, never accession ranges or beta profiles."""
    if dataset == "GSE193795":
        candidates = [c for c in metadata if c.startswith("characteristic_subject_status")]
    elif dataset == "GSE42774":
        candidates = [c for c in metadata if c.startswith("characteristic_disease_state")]
    else:
        raise ValueError("Unsupported external dataset: outcome mapping requires a reviewed adapter")
    if len(candidates) != 1:
        raise ValueError(f"Missing/ambiguous declared phenotype in {dataset}: {candidates}")
    def label(text):
        text = str(text).strip().lower()
        if "pre" in text and ("hyper" in text or "ht" in text):
            return "preHT"
        if "healthy" in text or "control" in text or "normotensive" in text or "normortensive" in text:
            return "control"
        if "hyper" in text:
            return "HTN"
        raise ValueError(f"Unknown phenotype: {text}")
    return metadata[candidates[0]].map(label)


def covariate_report(metadata):
    age = [c for c in metadata if c.startswith("characteristic_age")]
    sex = [c for c in metadata if c.startswith(("characteristic_sex", "characteristic_gender"))]
    return {"individual_age_field": age, "individual_sex_field": sex,
            "age_adjusted": False, "sex_adjusted": False,
            "reason": "No covariate adjustment prespecified; group-level publication demographics are not individual covariates."}


def common_probe_mask(primary_probes, external_probes):
    """Only IDs enter this function. Neither outcomes nor methylation values can enter."""
    if len(set(primary_probes)) != len(primary_probes) or len(set(external_probes)) != len(external_probes):
        raise ValueError("Duplicate probe IDs")
    return np.array(sorted(set(primary_probes).intersection(external_probes)), dtype=str)


def read_family(path, primary_probes, primary_samples):
    """Stream official SOFT: annotation and declared per-sample Detection P values."""
    lookup = {p: i for i, p in enumerate(primary_probes)}
    samples = {s: i for i, s in enumerate(primary_samples)}
    detection = np.full((len(primary_probes), len(primary_samples)), np.nan, dtype=np.float32)
    annotation = []
    table = None
    header = None
    sample = None
    detected_samples = set()
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("^SAMPLE ="):
                sample = line.split("=", 1)[1].strip()
            if line.startswith("!platform_table_begin"):
                table, header = "platform", None
                continue
            if line.startswith("!sample_table_begin"):
                table, header = "sample", None
                continue
            if line.startswith(("!platform_table_end", "!sample_table_end")):
                table = None
                continue
            if table is None:
                continue
            values = line.rstrip("\r\n").split("\t")
            if header is None:
                header = [x.strip('"') for x in values]
                if table == "sample":
                    dcols = [i for i, x in enumerate(header) if "detect" in x.lower()]
                    dcol = dcols[0] if len(dcols) == 1 else None
                continue
            probe = values[0].strip('"')
            if probe not in lookup:
                continue
            if table == "platform":
                kept = {k: values[i].strip('"') for i, k in enumerate(header) if i < len(values) and (i == 0 or any(t in k.lower() for t in ["gene", "chr", "coordinate", "position", "mapinfo", "island", "relation"]))}
                kept["probe_id"] = probe
                annotation.append(kept)
            elif sample in samples and dcol is not None:
                value = values[dcol].strip('"')
                detection[lookup[probe], samples[sample]] = float(value) if value not in ("null", "NA", "", "NaN") else np.nan
                detected_samples.add(sample)
    if detected_samples != set(primary_samples):
        raise ValueError(f"Detection P unavailable for samples: {set(primary_samples) - detected_samples}")
    ann = pd.DataFrame(annotation).drop_duplicates("probe_id").set_index("probe_id")
    return detection, ann


def prepare(root, detection_threshold=0.01):
    raw = root / "data" / "raw"
    out = root / "data" / "prepared"
    out.mkdir(parents=True, exist_ok=True)
    primary, meta = read_matrix(raw / "GSE193795_series_matrix.txt.gz")
    external, emeta = read_matrix(raw / "GSE42774_series_matrix.txt.gz")
    # Freeze measurement identity compatibility before either outcome is mapped.
    common = common_probe_mask(primary.index.tolist(), external.index.tolist())
    pd.Series(common, name="probe_id").to_csv(out / "common_probes.tsv", sep="\t", index=False)
    mask_hash = sha256(out / "common_probes.tsv")
    detection, ann = read_family(raw / "GSE193795_family.soft.gz", primary.index, primary.columns)
    external_detection, external_ann = read_family(raw / "GSE42774_family.soft.gz", external.index, external.columns)
    values = primary.to_numpy(copy=True)
    invalid = (detection > detection_threshold) | ~np.isfinite(detection) | (values < 0) | (values > 1)
    values[invalid] = np.nan
    primary = pd.DataFrame(values, index=primary.index, columns=primary.columns)
    if ((external.to_numpy() < 0) | (external.to_numpy() > 1)).any():
        raise ValueError("External matrix is not beta values within [0,1]")
    external_values = external.to_numpy(copy=True)
    external_invalid = (external_detection > detection_threshold) | ~np.isfinite(external_detection)
    external_values[external_invalid] = np.nan
    external = pd.DataFrame(external_values, index=external.index, columns=external.columns)
    meta["phenotype"] = phenotype(meta, "GSE193795")
    emeta["phenotype"] = phenotype(emeta, "GSE42774")
    if meta.phenotype.value_counts().to_dict() != {"HTN": 44, "preHT": 44, "control": 44}:
        raise ValueError(f"Unexpected primary cohort: {meta.phenotype.value_counts().to_dict()}")
    if emeta.phenotype.value_counts().to_dict() != {"HTN": 8, "control": 8}:
        raise ValueError("Unexpected external cohort size")
    # Chips are reported only if present in deposited metadata. Never invent age/sex.
    def chip(row):
        matches = re.findall(r"(\d{10,})_(R\d+C\d+)", " ".join(row.fillna("").astype(str)))
        return matches[0][0] if matches else "unavailable"
    meta["chip_id"] = meta.apply(chip, axis=1)
    meta.to_csv(out / "primary_metadata.tsv", sep="\t")
    emeta.to_csv(out / "external_metadata.tsv", sep="\t")
    ann.to_csv(out / "probe_annotation.tsv", sep="\t")
    external_ann.to_csv(out / "external_probe_annotation.tsv", sep="\t")
    primary_ids = primary.columns.to_numpy(dtype=str)
    labels = meta.phenotype.to_numpy(dtype=str)
    np.savez_compressed(out / "matrices.npz", primary=primary.to_numpy().T, external=external.loc[common].to_numpy().T,
                        probe_ids=primary.index.to_numpy(dtype=str), common_probes=common,
                        primary_samples=primary_ids, external_samples=external.columns.to_numpy(dtype=str),
                        primary_phenotypes=labels, external_y=(emeta.phenotype == "HTN").to_numpy(dtype=np.int8))
    sample_qc = pd.DataFrame({"sample_id": primary.columns, "missing_fraction_after_detection_mask": invalid.mean(axis=0),
                              "detection_failure_fraction": ((detection > detection_threshold) | ~np.isfinite(detection)).mean(axis=0)})
    sample_qc.to_csv(out / "sample_qc.tsv", sep="\t", index=False)
    pd.DataFrame({"sample_id": external.columns, "detection_failure_fraction": external_invalid.mean(axis=0)}).to_csv(out / "external_sample_qc.tsv", sep="\t", index=False)
    pd.crosstab(meta.phenotype, meta.chip_id).to_csv(out / "phenotype_by_chip.tsv", sep="\t")
    profile_hashes = [hashlib.sha256(np.nan_to_num(row, nan=-1).tobytes()).hexdigest() for row in primary.to_numpy().T]
    if len(set(profile_hashes)) != len(profile_hashes):
        raise ValueError("Duplicate primary beta profiles: do not run sample-wise CV")
    report = {"primary_samples": len(primary.columns), "external_samples": len(external.columns), "primary_probes": len(primary),
              "common_probe_count": len(common), "common_mask_sha256": mask_hash, "mask_uses_outcomes": False,
              "detection_threshold": detection_threshold, "external_detection_p": "Observed per-sample Detection Pval from official GSE42774 family SOFT; fixed per-measurement mask only",
              "primary_covariates": covariate_report(meta), "external_covariates": covariate_report(emeta),
              "biological_unit": "one deposited GSM; duplicate profiles checked; undisclosed relatedness cannot be excluded",
              "sample_exclusion": "No sample excluded; detection-failed measurements become missing; probe QC learned inside CV",
              "external_interpretation": "Exploratory transfer across platform, age range, ancestry and sex; n16 age-matched male study",
              "no_joint_normalization": True, "versions": versions()}
    write_json(out / "qc_manifest.json", report)
    print(json.dumps({k: report[k] for k in ["primary_samples", "external_samples", "primary_probes", "common_probe_count"]}), flush=True)
    return report
