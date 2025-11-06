# -*- coding: utf-8 -*-
"""
setup_data_pipeline.py
----------------------
End-to-end helper that
  1) downloads the raw FDC / CTD / GWAS sources (data/raw/*)
  2) runs preprocess.py to build data/processed/*.parquet

Usage:
    python setup_data_pipeline.py --all
    python setup_data_pipeline.py --datasets ctd gwas --preprocess-only   # skip downloads
    python setup_data_pipeline.py --download-only                         # just fetch raw files

The script simply wraps existing utilities (`download_data.py`, `preprocess.py`)
so you do not have to remember the sequence of commands.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from download_data import download_fdc, download_ctd, download_gwas

ROOT_DIR = Path(__file__).resolve().parent
RAW_DIR = ROOT_DIR / "data" / "raw"
PROC_DIR = ROOT_DIR / "data" / "processed"


def run_preprocess(limit_rows: int | None = None, fuzzy_threshold: int | None = None) -> None:
    """Invoke preprocess.py via a subprocess so CLI behaviour stays consistent."""
    cmd = [sys.executable, str(ROOT_DIR / "preprocess.py")]
    if limit_rows is not None:
        cmd += ["--limit-rows", str(limit_rows)]
    if fuzzy_threshold is not None:
        cmd += ["--fuzzy-threshold", str(fuzzy_threshold)]

    print(f"[+] Running preprocess: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=str(ROOT_DIR))


def ensure_raw_directories() -> None:
    (RAW_DIR / "fdc").mkdir(parents=True, exist_ok=True)
    (RAW_DIR / "ctd").mkdir(parents=True, exist_ok=True)
    (RAW_DIR / "gwas").mkdir(parents=True, exist_ok=True)


def download_selected(datasets: Iterable[str]) -> None:
    ensure_raw_directories()
    ds = {d.lower() for d in datasets}

    if "fdc" in ds or "all" in ds:
        print("[*] Downloading FDC full dataset …")
        download_fdc(RAW_DIR / "fdc")
        print("[✓] FDC download completed.")

    if "ctd" in ds or "all" in ds:
        print("[*] Downloading CTD reports …")
        download_ctd(RAW_DIR / "ctd")
        print("[✓] CTD download completed.")

    if "gwas" in ds or "all" in ds:
        print("[*] Downloading GWAS catalog associations …")
        download_gwas(RAW_DIR / "gwas")
        print("[✓] GWAS download completed.")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Download raw datasets and run preprocessing in one step.")
    ap.add_argument(
        "--datasets",
        nargs="+",
        default=["fdc", "ctd", "gwas"],
        help="Datasets to download (fdc ctd gwas). Use 'all' for every source."
    )
    ap.add_argument(
        "--download-only",
        action="store_true",
        help="Only download raw files; do not run preprocess.py"
    )
    ap.add_argument(
        "--preprocess-only",
        action="store_true",
        help="Skip downloads and only execute preprocess.py"
    )
    ap.add_argument(
        "--limit-rows",
        type=int,
        default=None,
        help="Optional row cap used when running preprocess.py (useful for smoke tests)."
    )
    ap.add_argument(
        "--fuzzy-threshold",
        type=int,
        default=None,
        help="Forwarded to preprocess.py --fuzzy-threshold (defaults to value inside preprocess)."
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()

    if not args.preprocess_only:
        download_selected(args.datasets)
    else:
        print("[!] Skipping downloads (preprocess-only mode).")

    if not args.download_only:
        PROC_DIR.mkdir(parents=True, exist_ok=True)
        run_preprocess(limit_rows=args.limit_rows, fuzzy_threshold=args.fuzzy_threshold)
        print(f"[✓] Preprocessing finished. Outputs in: {PROC_DIR}")
    else:
        print("[!] Skipping preprocessing (download-only mode).")


if __name__ == "__main__":
    main()
