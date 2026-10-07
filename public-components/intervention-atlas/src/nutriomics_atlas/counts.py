from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import TextIO

from .geo import parse_soft_samples
DATA_DIR=Path('data')


FILES = (
    "GSE127530_combinedCounts.txt.gz",
    "GSE127530_fixed_combinedCounts.txt.gz",
)
SAMPLE_LABEL = re.compile(r"^(S\d+)-D(\d+)-(Fast|3hr|6hr)_S(\d+)$", re.IGNORECASE)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_manifest(data_dir: Path) -> dict[str, dict]:
    path = data_dir / "raw" / "geo" / "download_manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    return {row["name"]: row for row in manifest.get("files", []) if row.get("name")}


def _open_gzip(path: Path) -> TextIO:
    return gzip.open(path, "rt", encoding="utf-8", newline="")


def _parse_header(line: str, path: Path) -> tuple[str, ...]:
    fields = tuple(line.rstrip("\r\n").split("\t"))
    if len(fields) != 45 or len(set(fields)) != len(fields):
        raise ValueError(f"Expected 45 unique sample labels in {path.name}; got {len(fields)} fields")
    labels = []
    for label in fields:
        match = SAMPLE_LABEL.fullmatch(label)
        if not match:
            raise ValueError(f"Unrecognized sample label {label!r} in {path.name}")
        labels.append(
            {
                "sample_label": label,
                "subject_id": match.group(1).upper(),
                "study_day": int(match.group(2)),
                "timepoint": match.group(3).lower(),
                "library_suffix": int(match.group(4)),
            }
        )
    combinations = {(row["subject_id"], row["study_day"], row["timepoint"]) for row in labels}
    subjects = {row["subject_id"] for row in labels}
    days = {row["study_day"] for row in labels}
    times = {row["timepoint"] for row in labels}
    if (len(subjects), len(days), len(times), len(combinations)) != (5, 3, 3, 45):
        raise ValueError(
            "The parsed sample labels do not describe the expected 5 people × 3 days × 3 timepoints"
        )
    return fields


def _validate_rows(handle: TextIO, path: Path) -> tuple[int, str, dict[str, int]]:
    digest = hashlib.sha256()
    row_count = 0
    seen: dict[str, int] = {}
    for line_number, line in enumerate(handle, start=2):
        digest.update(line.encode("utf-8"))
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != 46:
            raise ValueError(f"Row {line_number} has {len(fields)} fields (expected 46) in {path.name}")
        feature = fields[0]
        if not feature:
            raise ValueError(f"Empty gene symbol at row {line_number} in {path.name}")
        seen[feature] = seen.get(feature, 0) + 1
        try:
            values = [int(value) for value in fields[1:]]
        except ValueError as error:
            raise ValueError(f"Non-integer count at row {line_number} in {path.name}") from error
        if any(value < 0 for value in values):
            raise ValueError(f"Negative count at row {line_number} in {path.name}")
        row_count += 1
    if row_count < 10_000:
        raise ValueError(f"Only {row_count} gene rows found in {path.name}")
    duplicates = {feature: count for feature, count in seen.items() if count > 1}
    return row_count, digest.hexdigest(), duplicates


def _indexed_rows(path: Path) -> dict[tuple[str, int], list[str]]:
    with _open_gzip(path) as handle:
        next(handle)
        occurrences: dict[str, int] = {}
        rows = {}
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            symbol = fields[0]
            occurrence = occurrences.get(symbol, 0)
            occurrences[symbol] = occurrence + 1
            rows[(symbol, occurrence)] = fields[1:]
    return rows


def _compare_payloads(data_dir: Path, headers: tuple[str, ...]) -> dict:
    folder = data_dir / "raw" / "geo" / "GSE127530"
    original = _indexed_rows(folder / FILES[0])
    fixed = _indexed_rows(folder / FILES[1])
    original_keys = set(original)
    fixed_keys = set(fixed)
    original_only = sorted(original_keys - fixed_keys)
    fixed_only = sorted(fixed_keys - original_keys)
    shared = sorted(original_keys & fixed_keys)
    changed_rows = 0
    changed_cells = 0
    examples = []
    for key in shared:
        old = original[key]
        new = fixed[key]
        changed = [index for index, (a, b) in enumerate(zip(old, new, strict=True)) if a != b]
        if not changed:
            continue
        changed_rows += 1
        changed_cells += len(changed)
        if len(examples) < 20:
            examples.append(
                {
                    "gene_symbol": key[0],
                    "duplicate_symbol_occurrence": key[1] + 1,
                    "changed_cell_count": len(changed),
                    "changed_samples": [headers[index] for index in changed[:10]],
                    "original_values": [old[index] for index in changed[:10]],
                    "fixed_values": [new[index] for index in changed[:10]],
                }
            )
    return {
        "identical_feature_value_rows": not original_only and not fixed_only and changed_rows == 0,
        "original_only_row_count": len(original_only),
        "fixed_only_row_count": len(fixed_only),
        "changed_feature_row_count": changed_rows,
        "changed_count_cell_count": changed_cells,
        "original_only_examples": [symbol for symbol, _ in original_only[:20]],
        "fixed_only_examples": [symbol for symbol, _ in fixed_only[:20]],
        "original_only_rows": [
            {"gene_symbol": symbol, "duplicate_symbol_occurrence": occurrence + 1}
            for symbol, occurrence in original_only
        ],
        "changed_row_examples": examples,
    }


def audit(data_dir: Path = DATA_DIR, output_dir: Path | None = None) -> dict:
    data_dir = data_dir.expanduser().resolve()
    folder = data_dir / "raw" / "geo" / "GSE127530"
    manifest = _read_manifest(data_dir)
    records = []
    headers: list[tuple[str, ...]] = []
    uncompressed_hashes = []
    for name in FILES:
        path = folder / name
        if not path.is_file() or name not in manifest:
            raise FileNotFoundError(f"Missing count file or download manifest record: {path}")
        expected = manifest[name].get("sha256")
        actual = sha256_file(path)
        if not expected or expected != actual:
            raise ValueError(f"Download manifest checksum mismatch for {path}")
        with _open_gzip(path) as handle:
            header = _parse_header(handle.readline(), path)
            rows, data_digest, duplicate_symbols = _validate_rows(handle, path)
        headers.append(header)
        uncompressed_hashes.append(data_digest)
        records.append(
            {
                "name": name,
                "compressed_sha256": actual,
                "uncompressed_data_sha256": data_digest,
                "gene_count": rows,
                "duplicate_symbol_rows": sum(count - 1 for count in duplicate_symbols.values()),
                "duplicate_symbols": duplicate_symbols,
                "sample_count": len(header),
                "bytes": path.stat().st_size,
                "source_url": manifest[name].get("url"),
            }
        )
    if headers[0] != headers[1]:
        raise ValueError("The original and fixed count files have different feature or sample order; manual review required")
    payload_comparison = _compare_payloads(data_dir, headers[1])
    identical_payloads = (
        uncompressed_hashes[0] == uncompressed_hashes[1]
        and payload_comparison["identical_feature_value_rows"]
    )
    expected_cleanup = Counter(
        {
            "N_ambiguous": 1,
            "N_multimapping": 1,
            "N_noFeature": 1,
            "N_unmapped": 1,
            "TTTY17B": 2,
        }
    )
    observed_removed = Counter(payload_comparison["original_only_examples"])
    verified_fixed_cleanup = (
        payload_comparison["fixed_only_row_count"] == 0
        and payload_comparison["changed_feature_row_count"] == 0
        and observed_removed == expected_cleanup
    )
    approved_release = identical_payloads or verified_fixed_cleanup
    family_soft_path = data_dir / "raw" / "geo" / "GSE127530" / "GSE127530_family.soft.gz"
    family_samples = parse_soft_samples(family_soft_path) if family_soft_path.is_file() else {}
    title_index: dict[str, list[tuple[str, str]]] = {}
    for sample_accession, sample in family_samples.items():
        title = str(sample.get("title") or "").strip()
        if title:
            title_index.setdefault(title.casefold(), []).append((sample_accession, title))
    sample_crosswalk = []
    for label in headers[1]:
        base_label = re.sub(r"_S\d+$", "", label, flags=re.IGNORECASE)
        matches = title_index.get(base_label.casefold(), [])
        sample_crosswalk.append(
            {
                "sample_label": label,
                "subject_id": SAMPLE_LABEL.fullmatch(label).group(1).upper(),
                "study_day": int(SAMPLE_LABEL.fullmatch(label).group(2)),
                "timepoint": SAMPLE_LABEL.fullmatch(label).group(3).lower(),
                "library_suffix": int(SAMPLE_LABEL.fullmatch(label).group(4)),
                "geo_sample_accession": matches[0][0] if len(matches) == 1 else None,
                "geo_sample_title": matches[0][1] if len(matches) == 1 else None,
                "geo_sample_match_count": len(matches),
            }
        )
    mapped_accessions = [row["geo_sample_accession"] for row in sample_crosswalk if row["geo_sample_accession"]]
    crosswalk_status = (
        "complete"
        if len(sample_crosswalk) == 45 and len(mapped_accessions) == 45 and len(set(mapped_accessions)) == 45
        else "pending_family_soft_download"
        if not family_soft_path.is_file()
        else "unmatched_or_ambiguous"
    )

    report = {
        "accession": "GSE127530",
        "audit_version": "gse127530-count-audit-v2",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "status": (
            "identical_payloads"
            if identical_payloads
            else "verified_fixed_cleanup"
            if verified_fixed_cleanup
            else "payloads_differ_manual_review"
        ),
        "selected_analysis_file": FILES[1] if approved_release else None,
        "release_decision": (
            "Use the fixed file: all retained gene counts are unchanged; the original-only rows are four non-gene count summaries and two duplicate TTTY17B rows."
            if verified_fixed_cleanup
            else "Use the fixed file because both payloads are identical."
            if identical_payloads
            else "Analysis is blocked pending inspection of the original-versus-fixed differences."
        ),
        "payload_comparison": payload_comparison,
        "sample_design": {
            "people": 5,
            "study_days_per_person": 3,
            "timepoints_per_day": 3,
            "timepoints": ["fast", "3hr", "6hr"],
            "samples": sample_crosswalk,
            "geo_sample_crosswalk_status": crosswalk_status,
            "geo_family_soft_sample_records": len(family_samples),
            "independent_person_count": 5,
        },
        "files": records,
        "interpretation_limit": (
            "Column labels identify person, visit day, and timepoint. A food-component-specific causal effect cannot be inferred from this mixed high-fat challenge."
        ),
    }
    if output_dir is not None:
        output_dir = output_dir.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        destination = output_dir / "gse127530_counts_audit.json"
        temporary = destination.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, destination)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify GEO GSE127530 count files and sample labels")
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    print(json.dumps(audit(args.data_dir, args.output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

