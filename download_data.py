# -*- coding: utf-8 -*-
import argparse
import sys
import gzip
import shutil
import zipfile
from pathlib import Path
from typing import Optional
import requests
from tqdm import tqdm

CTD_REPORT_BASE = "https://ctdbase.org/reports"  # reports 경로는 논문/도구 문서에 예시 존재
# 예시 파일명: CTD_chem_gene_ixns.tsv.gz, CTD_genes_diseases.tsv.gz, CTD_chemicals.tsv.gz, CTD_diseases.tsv.gz

# 2025-04 FDC Full Download (페이지의 "Full Download of All Data Types 04/2025" 링크를 확인)
FDC_FULL_URL = "https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_csv_2025-04-24.zip"

# GWAS Catalog: latest 심볼릭 폴더를 HTTPS로 접근 (FTP 미러를 HTTPS로)
GWAS_BASE_HTTPS = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases"
GWAS_LATEST_DIR = f"{GWAS_BASE_HTTPS}/latest"
GWAS_ONTOLOGY_ZIP = f"{GWAS_LATEST_DIR}/gwas-catalog-associations_ontology-annotated-full.zip"
GWAS_ONTOLOGY_TSV = "gwas-catalog-associations_ontology-annotated.tsv"

def stream_download(url: str, dest: Path, chunk_size: int = 1024 * 1024) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        total = int(r.headers.get("Content-Length", 0))
        with tqdm(total=total, unit='B', unit_scale=True, desc=f"GET {url.split('/')[-1]}") as pbar:
            with open(dest, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))

def gunzip_file(src_gz: Path, dst_path: Optional[Path] = None) -> Path:
    if dst_path is None:
        dst_path = src_gz.with_suffix('')  # drop .gz
    with gzip.open(src_gz, 'rb') as f_in, open(dst_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return dst_path

def unzip_file(zip_path: Path, extract_dir: Path) -> None:
    extract_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, 'r') as zf:
        zf.extractall(extract_dir)

def download_fdc(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    zip_path = out_dir / Path(FDC_FULL_URL).name
    if not zip_path.exists():
        stream_download(FDC_FULL_URL, zip_path)
    # 압축 해제
    unzip_dir = out_dir / "FoodData_Central_csv_2025-04-24"
    if not unzip_dir.exists():
        unzip_file(zip_path, unzip_dir)

def download_ctd(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    # 보고서 경로 기준 파일들
    files = [
        "CTD_chem_gene_ixns.tsv.gz",
        "CTD_genes_diseases.tsv.gz",
        "CTD_chemicals.tsv.gz",
        "CTD_diseases.tsv.gz"
    ]
    for fname in files:
        url = f"{CTD_REPORT_BASE}/{fname}"
        gz_path = out_dir / fname
        if not gz_path.exists():
            try:
                stream_download(url, gz_path)
            except Exception as e:
                print(f"[WARN] Failed to download {url}: {e}", file=sys.stderr)
                print("       수동으로 다운로드해 data/raw/ctd/에 배치하세요.")
                continue
        # 해제본은 전처리 단계에서 on-the-fly로 읽어도 되므로 여기서는 해제 생략 가능

def _extract_gwas_zip(zip_path: Path, dest_tsv: Path) -> None:
    with zipfile.ZipFile(zip_path, 'r') as zf:
        candidates = [name for name in zf.namelist() if name.endswith(".tsv")]
        if not candidates:
            raise ValueError(f"{zip_path} 안에 TSV 파일이 없습니다: {zf.namelist()}")
        # 선호: 정확한 파일명, 없으면 첫번째 TSV
        target = next((name for name in candidates if name.endswith(GWAS_ONTOLOGY_TSV)), candidates[0])
        dest_tsv.parent.mkdir(parents=True, exist_ok=True)
        extracted_path = Path(zf.extract(target, dest_tsv.parent))
        shutil.move(str(extracted_path), dest_tsv)
        # 남은 빈 디렉토리 정리
        extra_parent = extracted_path.parent
        while extra_parent != dest_tsv.parent and extra_parent.exists():
            try:
                extra_parent.rmdir()
            except OSError:
                break
            extra_parent = extra_parent.parent

def download_gwas(out_dir: Path, url: str = GWAS_ONTOLOGY_ZIP) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / GWAS_ONTOLOGY_TSV
    if dest.exists():
        return

    def warn_failure(err_msg: str) -> None:
        print(f"[WARN] GWAS latest HTTPS 접근 실패: {err_msg}", file=sys.stderr)
        print("      대안) ftp 링크 또는 파일을 수동 다운로드하여 data/raw/gwas/에 두세요.", file=sys.stderr)
        print("      예) ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/{release_date}/gwas-catalog-associations_ontology-annotated.tsv", file=sys.stderr)

    try:
        if url.endswith(".zip"):
            zip_path = out_dir / Path(url).name
            if not zip_path.exists():
                stream_download(url, zip_path)
            _extract_gwas_zip(zip_path, dest)
        else:
            stream_download(url, dest)
    except Exception as e:
        fallback_url = f"{GWAS_LATEST_DIR}/{GWAS_ONTOLOGY_TSV}"
        if url != fallback_url:
            try:
                stream_download(fallback_url, dest)
                return
            except Exception as e2:
                warn_failure(f"{e} (zip) / {e2} (direct TSV)")
                return
        warn_failure(str(e))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="project root")
    ap.add_argument("--datasets", nargs="+", default=["fdc", "ctd", "gwas"],
                    help="which datasets to download")
    args = ap.parse_args()

    root = Path(args.root)
    raw_dir = root / "data" / "raw"

    if "fdc" in args.datasets:
        print("[*] Downloading FDC ...")
        download_fdc(raw_dir / "fdc")
        print("[✓] FDC done.")

    if "ctd" in args.datasets:
        print("[*] Downloading CTD ...")
        download_ctd(raw_dir / "ctd")
        print("[✓] CTD done.")

    if "gwas" in args.datasets:
        print("[*] Downloading GWAS Catalog (associations) ...")
        download_gwas(raw_dir / "gwas")
        print("[✓] GWAS done.")

if __name__ == "__main__":
    main()
