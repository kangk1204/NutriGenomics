"""Export receipts only: synthetic tables, no model fitting or private data.

The exact export AST is compiled so this bounded CI test needs pandas/numpy but
does not import evaluation's model/scientific dependencies. Native file I/O
helpers are loaded unchanged; dependency-version metadata alone is stubbed.
"""
import ast
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


COMPONENT = Path(__file__).resolve().parents[1]
SOURCE = COMPONENT / "src/nutriomics_methylation/evaluation.py"
IO_SOURCE = COMPONENT / "src/nutriomics_methylation/io.py"


@pytest.fixture
def export_with_io():
    spec = importlib.util.spec_from_file_location("export_test_native_io", IO_SOURCE)
    io = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(io)
    node = next(node for node in ast.parse(SOURCE.read_text(encoding="utf-8")).body
                if isinstance(node, ast.FunctionDef) and node.name == "export")
    scope = {"pd": pd, "sha256": io.sha256, "write_json": io.write_json,
             "versions": lambda: {"synthetic_export_test": True}}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(SOURCE), "exec"), scope)
    return scope["export"], io


def write_synthetic_inputs(root, io, auroc):
    for name in ("results", "data/prepared", "src"):
        (root / name).mkdir(parents=True, exist_ok=True)
    rows = [{"analysis": "synthetic_only", "n": 8, "metric": metric,
             "estimate": estimate, "ci95_lower": 0.1, "ci95_upper": 0.9}
            for metric, estimate in (("auroc", auroc), ("average_precision", 0.7), ("brier", 0.15))]
    pd.DataFrame(rows).to_csv(root / "results/performance.tsv", sep="\t", index=False)
    for name in ("results/metrics_all.json", "results/training_config.json",
                 "data/source_manifest.json", "data/prepared/qc_manifest.json"):
        io.write_json(root / name, {"synthetic_fixture": True})
    (root / "src/synthetic.py").write_text("# synthetic fixture only\n", encoding="utf-8")


def assert_current_manifest(root, io, expected=None):
    target = root / "docs/validation"
    manifest = json.loads((target / "artifact_manifest.json").read_text(encoding="utf-8"))
    entries = {record["path"]: record for record in manifest["files"]}
    if expected is None:
        expected = {path.relative_to(root).as_posix() for path in target.iterdir()
                    if path.is_file() and path.name != "artifact_manifest.json"}
    assert set(entries) == expected
    assert "docs/validation/SUMMARY.md" in entries
    for name, record in entries.items():
        assert record["sha256"] == io.sha256(root / name), name
        assert record["bytes"] == (root / name).stat().st_size, name
    assert manifest["not_clinical"] is True
    assert manifest["no_causal_dietary_claim"] is True
    assert manifest["code_files"] == [{"path": "src/synthetic.py", "sha256": io.sha256(root / "src/synthetic.py")}]


def test_first_export_records_current_summary(tmp_path, export_with_io):
    export, io = export_with_io
    write_synthetic_inputs(tmp_path, io, auroc=0.75)
    assert export(tmp_path) == tmp_path / "docs/validation"
    assert_current_manifest(tmp_path, io)


def test_repeated_export_rehashes_changed_summary(tmp_path, export_with_io):
    export, io = export_with_io
    write_synthetic_inputs(tmp_path, io, auroc=0.75)
    export(tmp_path)
    old_summary_hash = io.sha256(tmp_path / "docs/validation/SUMMARY.md")
    write_synthetic_inputs(tmp_path, io, auroc=0.61)
    export(tmp_path)
    assert io.sha256(tmp_path / "docs/validation/SUMMARY.md") != old_summary_hash
    assert_current_manifest(tmp_path, io)



def write_current_optional_inputs(root):
    directory = root / "results/synthetic_ridge/inner"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "marker_stability.tsv").write_text("marker_id\tfrequency\nsynthetic_only\t0.5\n", encoding="utf-8")
    (directory / "current_predictions.tsv").write_text("synthetic_id\tprobability\nfixture_only\t0.5\n", encoding="utf-8")
    (directory / "roc_calibration.svg").write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>synthetic only</text></svg>\n', encoding="utf-8")
    return directory


def expected_current_paths():
    names = ("SUMMARY.md", "performance.tsv", "metrics_all.json", "training_config.json",
             "source_manifest.json", "qc_manifest.json", "synthetic_ridge_inner_marker_stability.tsv",
             "synthetic_ridge_inner_current_predictions.tsv", "synthetic_ridge_inner_roc_calibration.svg")
    return {"docs/validation/" + name for name in names}


def preserved_file_state(path):
    return path.read_bytes(), path.stat().st_size, path.stat().st_mtime_ns


def test_export_omits_seeded_obsolete_file_without_touching_it(tmp_path, export_with_io):
    export, io = export_with_io
    write_synthetic_inputs(tmp_path, io, auroc=0.75)
    write_current_optional_inputs(tmp_path)
    obsolete = tmp_path / "docs/validation/OBSOLETE_KEEP.txt"
    obsolete.parent.mkdir(parents=True, exist_ok=True)
    obsolete.write_bytes(b"synthetic obsolete artifact; preserve exactly\n")
    before = preserved_file_state(obsolete)
    export(tmp_path)
    assert_current_manifest(tmp_path, io, expected_current_paths())
    assert preserved_file_state(obsolete) == before
    manifest = json.loads((obsolete.parent / "artifact_manifest.json").read_text(encoding="utf-8"))
    assert manifest["inventory_scope"] == "files_written_by_current_export"


def test_repeat_export_omits_removed_optional_input_without_touching_old_copy(tmp_path, export_with_io):
    export, io = export_with_io
    write_synthetic_inputs(tmp_path, io, auroc=0.75)
    directory = write_current_optional_inputs(tmp_path)
    retired = directory / "retired_predictions.tsv"
    retired.write_bytes(b"synthetic_id\tprobability\nretired_fixture_only\t0.3\n")
    export(tmp_path)
    previous = tmp_path / "docs/validation/synthetic_ridge_inner_retired_predictions.tsv"
    assert_current_manifest(tmp_path, io, expected_current_paths() | {previous.relative_to(tmp_path).as_posix()})
    before = preserved_file_state(previous)
    old_summary_hash = io.sha256(tmp_path / "docs/validation/SUMMARY.md")
    retired.unlink()
    write_synthetic_inputs(tmp_path, io, auroc=0.61)
    export(tmp_path)
    assert not retired.exists()
    assert preserved_file_state(previous) == before
    assert io.sha256(tmp_path / "docs/validation/SUMMARY.md") != old_summary_hash
    assert_current_manifest(tmp_path, io, expected_current_paths())
