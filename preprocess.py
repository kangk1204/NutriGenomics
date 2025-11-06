# -*- coding: utf-8 -*-
"""
preprocess.py (FINAL hardened)
- API 없이 FDC/CTD/GWAS 원시 파일을 읽어 data/processed/* 로 전처리 저장
- CTD 헤더(#) 처리/스키마 변동에 강함. 디버깅 출력 풍부.

산출물:
  data/processed/
    ├─ fdc_foods.parquet
    ├─ fdc_food_nutrient.parquet
    ├─ ctd_chem_synonyms.parquet
    ├─ ctd_chem_gene_human.parquet
    ├─ ctd_gene_disease.parquet
    ├─ ctd_diseases_metabolic.parquet
    ├─ fdc_to_ctd_chemical_map.parquet
    ├─ debug_fdc_to_ctd_mapping.csv
    └─ gwas_assoc_subset.parquet  (옵션)
"""

from __future__ import annotations
import argparse
import sys
import re
from pathlib import Path
from typing import Optional, Tuple, Dict, List
import pandas as pd
import numpy as np
from pandas.errors import ParserError
import gzip

# 퍼지 매칭
try:
    from rapidfuzz import process, fuzz
    HAS_RAPIDFUZZ = True
except Exception:
    HAS_RAPIDFUZZ = False

# -----------------------------
# 경로/환경
# -----------------------------
RAW_FDC_DIR = Path("data/raw/fdc")
RAW_CTD_DIR = Path("data/raw/ctd")
RAW_GWAS_DIR = Path("data/raw/gwas")
PROC_DIR = Path("data/processed")
OVERRIDES_PATH = Path("config/synonyms_override.csv")  # 선택(수동 매핑 보완)

MANUAL_NUTRIENT_CTDS = {
    "calcium, ca": ("MESH:D002115", "calcium"),
    "iron, fe": ("MESH:D007501", "iron"),
    "magnesium, mg": ("MESH:D008274", "magnesium"),
    "sodium, na": ("MESH:D012964", "sodium"),
    "potassium, k": ("MESH:D011188", "potassium"),
    "phosphorus, p": ("MESH:D010738", "phosphorus"),
    "zinc, zn": ("MESH:D015043", "zinc"),
    "copper, cu": ("MESH:D003315", "copper"),
    "manganese, mn": ("MESH:D008357", "manganese"),
    "selenium, se": ("MESH:D012640", "selenium"),
    "boron, b": ("MESH:D001911", "boron"),
    "chlorine, cl": ("MESH:D002701", "chlorine"),
    "chromium, cr": ("MESH:D002806", "chromium"),
    "cobalt, co": ("MESH:D003072", "cobalt"),
    "carbohydrate, by difference": ("MESH:D002245", "carbohydrates"),
    "carbohydrate, by summation": ("MESH:D002245", "carbohydrates"),
    "carbohydrate, other": ("MESH:D002245", "carbohydrates"),
    "ash": ("MESH:D000459", "ash"),
    "choline, total": ("MESH:D002813", "choline"),
    "choline, free": ("MESH:D002813", "choline"),
    "choline, from glycerophosphocholine": ("MESH:D002813", "glycerophosphocholine"),
    "choline, from phosphocholine": ("MESH:D002813", "phosphocholine"),
    "choline, from phosphotidyl choline": ("MESH:D002813", "phosphatidylcholine"),
    "choline, from sphingomyelin": ("MESH:D013162", "sphingomyelin")
}

# -----------------------------
# 로깅/유틸
# -----------------------------
def log(msg: str, level: str = "INFO"):
    print(f"[{level}] {msg}")

def df_info(df: pd.DataFrame, name: str, show_head: int = 3):
    mem_mb = df.memory_usage(deep=True).sum() / (1024 ** 2)
    log(f"{name}: shape={df.shape}, mem={mem_mb:,.2f} MB")
    if show_head > 0:
        with pd.option_context("display.max_columns", 120, "display.width", 160):
            print(df.head(show_head))

def ensure_dirs():
    PROC_DIR.mkdir(parents=True, exist_ok=True)

def normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    # (FDC용) 컬럼명 표준화: lower + 공백/하이픈 -> 언더스코어, 흔한 alias 보정
    def to_snake(c: str) -> str:
        return c.replace(" ", "_").replace("-", "_").lower()
    df = df.copy()
    df.columns = [to_snake(c) for c in df.columns]
    alias = {"unitname": "unit_name", "nutrientname": "name"}
    df.rename(columns={k: v for k, v in alias.items() if k in df.columns}, inplace=True)
    return df

def safe_to_numeric(ser: pd.Series) -> pd.Series:
    return pd.to_numeric(ser, errors="coerce")

# -----------------------------
# FDC 로딩/표준화
# -----------------------------
def find_fdc_paths(fdc_root: Path) -> Tuple[Path, Path, Path, Optional[Path]]:
    cands = list(fdc_root.glob("**/food.csv"))
    if not cands:
        raise FileNotFoundError("food.csv not found under data/raw/fdc. Ensure FDC zip is extracted.")
    food_csv = max(cands, key=lambda p: p.stat().st_size)
    base_dir = food_csv.parent
    nutrient_csv = base_dir / "nutrient.csv"
    food_nutrient_csv = base_dir / "food_nutrient.csv"
    food_category_csv = base_dir / "food_category.csv"
    log(f"Using FDC folder: {base_dir}")
    for f in (food_csv, nutrient_csv, food_nutrient_csv, food_category_csv):
        log(f" - {f.name}: {'exists' if f.exists() else 'MISSING'}")
    if not nutrient_csv.exists() or not food_nutrient_csv.exists():
        raise FileNotFoundError("nutrient.csv or food_nutrient.csv missing next to food.csv")
    return food_csv, nutrient_csv, food_nutrient_csv, (food_category_csv if food_category_csv.exists() else None)

def load_fdc_core(
    fdc_root: Path,
    limit_rows: Optional[int] = None,
    slim: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Optional[pd.DataFrame]]:
    food_csv, nutrient_csv, food_nutrient_csv, food_category_csv = find_fdc_paths(fdc_root)

    # food.csv
    head = pd.read_csv(food_csv, nrows=0, encoding="utf-8")
    log(f"food.csv columns({len(head.columns)}): {head.columns.tolist()}")

    food_usecols = [c for c in ["fdc_id", "data_type", "description", "food_category_id", "publication_date"] if c in head.columns]
    food = pd.read_csv(
        food_csv,
        usecols=food_usecols if food_usecols else None,
        dtype={c: "string" for c in ["data_type", "description"] if c in food_usecols},
        encoding="utf-8",
        low_memory=False,
        nrows=limit_rows
    )
    food = normalize_cols(food)
    if "food_category_id" in food.columns:
        nonnum = food["food_category_id"].astype(str).str.contains(r"[A-Za-z]", na=False).sum()
        if nonnum > 0:
            sample_bad = food.loc[
                food["food_category_id"].astype(str).str.contains(r"[A-Za-z]", na=False),
                "food_category_id"
            ].astype(str).head(5).tolist()
            log(f"Non-numeric values in food_category_id: count={nonnum}, samples={sample_bad}", "WARN")
            log("Coercing food_category_id to numeric (non-numeric -> <NA>)", "WARN")
        food["food_category_id"] = safe_to_numeric(food["food_category_id"]).astype("Int64")
    df_info(food, "FDC food")

    # nutrient.csv
    n_head = pd.read_csv(nutrient_csv, nrows=0, encoding="utf-8")
    nutrient_usecols = [c for c in ["id", "name", "unit_name", "rank"] if c in n_head.columns]
    nutrient = pd.read_csv(
        nutrient_csv,
        usecols=nutrient_usecols if nutrient_usecols else None,
        encoding="utf-8",
        low_memory=False,
        nrows=limit_rows
    )
    nutrient = normalize_cols(nutrient)
    must = ["id", "name"]
    missing = [c for c in must if c not in nutrient.columns]
    if missing:
        raise ValueError(f"nutrient.csv missing required columns: {missing}")
    if "unit_name" not in nutrient.columns:
        for c in nutrient.columns:
            if c.lower() in ("unitname", "unit_name"):
                nutrient.rename(columns={c: "unit_name"}, inplace=True)
                break
        if "unit_name" not in nutrient.columns:
            log("nutrient.csv has no 'unit_name'; mg conversion will be partial.", "WARN")
            nutrient["unit_name"] = np.nan
    df_info(nutrient, "FDC nutrient")

    # food_nutrient.csv (대용량) -> 필수 컬럼만 로딩
    fn_head = pd.read_csv(food_nutrient_csv, nrows=0, encoding="utf-8")
    fn_usecols = [c for c in ["fdc_id", "nutrient_id", "amount"] if c in fn_head.columns]
    if "amount" not in fn_usecols:
        for c in fn_head.columns:
            if c.lower().startswith("amount"):
                fn_usecols.append(c)
                break
    food_nutrient = pd.read_csv(
        food_nutrient_csv,
        usecols=fn_usecols if fn_usecols else None,
        encoding="utf-8",
        low_memory=False,
        nrows=limit_rows
    )
    food_nutrient = normalize_cols(food_nutrient)
    if "amount" not in food_nutrient.columns:
        cand_amt = [c for c in food_nutrient.columns if c.startswith("amount")]
        if cand_amt:
            food_nutrient.rename(columns={cand_amt[0]: "amount"}, inplace=True)
            log(f"Renamed {cand_amt[0]} -> amount", "WARN")
        else:
            raise ValueError("food_nutrient.csv has no 'amount' column")
    food_nutrient["amount"] = safe_to_numeric(food_nutrient["amount"])
    for col in ("fdc_id", "nutrient_id"):
        if col in food_nutrient.columns:
            food_nutrient[col] = safe_to_numeric(food_nutrient[col]).astype("Int64")
        else:
            raise ValueError(f"food_nutrient.csv missing required column: {col}")
    df_info(food_nutrient, "FDC food_nutrient (slim)" if slim else "FDC food_nutrient")

    if "unit_name" in nutrient.columns:
        uni = nutrient["unit_name"].astype(str).fillna("<NA>").value_counts().head(10)
        log(f"nutrient.unit_name top10: {uni.to_dict()}")

    # food_category.csv(옵션)
    food_cat = None
    if food_category_csv:
        try:
            food_cat = pd.read_csv(food_category_csv, encoding="utf-8", low_memory=False, nrows=limit_rows)
            food_cat = normalize_cols(food_cat)
            df_info(food_cat, "FDC food_category")
        except Exception as e:
            log(f"Failed to read food_category.csv: {e}", "WARN")

    return food, nutrient, food_nutrient, food_cat

def convert_units_to_mg(amount: float, unit: Optional[str]) -> float:
    if pd.isna(amount):
        return np.nan
    unit = (unit or "").strip().upper()
    if unit == "MG":
        return amount
    if unit == "UG":
        return amount / 1000.0
    if unit == "G":
        return amount * 1000.0
    return np.nan  # KCAL, IU, kJ 등은 환산 불가

def build_food_nutrient_table(
    food: pd.DataFrame, nutrient: pd.DataFrame, food_nutrient: pd.DataFrame, slim: bool = True
) -> pd.DataFrame:
    nut = nutrient.rename(columns={"id": "nutrient_id"})
    fn = food_nutrient.merge(nut[["nutrient_id", "name", "unit_name"]], on="nutrient_id", how="left")
    fn = fn.merge(food[["fdc_id", "description", "data_type"]], on="fdc_id", how="left")
    fn["amount_mg_per_100g"] = fn.apply(lambda r: convert_units_to_mg(r["amount"], r["unit_name"]), axis=1)
    if slim:
        keep = ["fdc_id", "description", "data_type", "nutrient_id", "name", "unit_name", "amount", "amount_mg_per_100g"]
        fn = fn[keep]
    df_info(fn, "FDC food_nutrient (joined)", show_head=5)
    return fn

# -----------------------------
# CTD 안전 로더 (주석/헤더 복원 + 폴백 스키마)
# -----------------------------
CTD_DEFAULT_SCHEMAS: Dict[str, List[str]] = {
    # 파일명에 포함되는 키 -> 기본 컬럼
    "CTD_chemicals": [
        "ChemicalName","ChemicalID","CasRN","Definition","ParentIDs",
        "TreeNumbers","ParentTreeNumbers","Synonyms"
    ],
    "CTD_diseases": [
        "DiseaseName","DiseaseID","AltDiseaseIDs","Definition","ParentIDs",
        "TreeNumbers","ParentTreeNumbers","Synonyms","SlimMappings"
    ],
    "CTD_genes_diseases": [
        "GeneSymbol","GeneID","DiseaseName","DiseaseID","DirectEvidence",
        "InferenceGeneSymbol","InferenceScore","OmimIDs","PubMedIDs"
    ],
    "CTD_chem_gene_ixns": [
        "ChemicalName","ChemicalID","GeneSymbol","GeneID","Organism",
        "OrganismID","Interaction","InteractionActions","PublicationIDs"
    ],
}

CTD_HEADER_TOKENS: Dict[str, List[str]] = {
    "CTD_chemicals": ["ChemicalName","ChemicalID","Synonyms"],
    "CTD_diseases":  ["DiseaseName","DiseaseID","TreeNumbers"],
    "CTD_genes_diseases": ["GeneSymbol","DiseaseName","DiseaseID"],
    "CTD_chem_gene_ixns": ["ChemicalID","GeneSymbol","InteractionActions"],
}

def _ctd_kind_by_name(path: Path) -> Optional[str]:
    name = path.name
    for key in CTD_DEFAULT_SCHEMAS.keys():
        if key in name:
            return key
    return None

def _peek_header_line(path: Path, sep: str, max_lines: int = 400) -> Tuple[Optional[int], Optional[List[str]], List[str]]:
    """
    파일 앞부분을 스캔하여 '# '로 시작하면서 핵심 토큰을 포함한 헤더 라인을 찾는다.
    반환: (헤더라인 인덱스, 컬럼명 리스트, 미리보기 라인들)
    """
    preview = []
    kind = _ctd_kind_by_name(path)
    tokens_need = CTD_HEADER_TOKENS.get(kind, [])

    open_fn = (lambda p: gzip.open(p, "rt", encoding="utf-8", errors="replace")) if str(path).lower().endswith(".gz") \
              else (lambda p: open(p, "rt", encoding="utf-8", errors="replace"))

    with open_fn(path) as fh:
        for idx, line in enumerate(fh):
            if idx < 50 and line.strip():
                preview.append(line.rstrip("\n"))
            if idx > max_lines:
                break
            if line.startswith("#"):
                cand = line.lstrip("#").strip()
                cols = cand.split(sep)
                # 헤더 후보 판단: 토큰이 모두 포함되는지
                if tokens_need and all(any(tok == c or tok in c for c in cols) for tok in tokens_need):
                    return idx, [c.strip() for c in cols], preview
                # 토큰 리스트가 없더라도 ChemicalID/GeneSymbol 등 대표 키워드가 있으면 사용
                if any(k in cand for k in ["ChemicalID","GeneSymbol","DiseaseID","Synonyms","Interaction"]):
                    return idx, [c.strip() for c in cols], preview
    return None, None, preview

def smart_read_ctd(path: Path) -> pd.DataFrame:
    """
    CTD .tsv(.gz)/.csv(.gz) 안전 로더:
      1) 파일 앞부분에서 '# 헤더' 라인을 탐지해 컬럼 복원
      2) 실패 시 파일 종류별 기본 스키마로 names 지정
      3) 데이터 부분부터 읽어오기(header=None, names=..., skiprows=헤더줄까지)
    """
    name_lower = path.name.lower()
    sep = "\t" if (name_lower.endswith(".tsv") or name_lower.endswith(".tsv.gz")) else ","

    header_idx, header_cols, preview = _peek_header_line(path, sep=sep, max_lines=400)
    if preview:
        log(f"Preview of {path.name} (first non-empty lines up to 50):", "INFO")
        for ln in preview[:10]:
            print("  >>", ln[:200])

    if header_cols is not None:
        log(f"Detected CTD header in {path.name} at line {header_idx}: {header_cols}")
        skiprows = header_idx + 1  # 헤더 줄까지 스킵
        try:
            df = pd.read_csv(
                path,
                sep=sep,
                compression="gzip" if name_lower.endswith(".gz") else None,
                header=None,
                names=header_cols,
                skiprows=skiprows,
                dtype="string",
                low_memory=False,
                engine="c",
            )
            # 끝부분에 또 다른 주석이 있어도 데이터 줄만 있음(일반적으로).
        except ParserError as e:
            log(f"C-engine ParserError for {path.name}: {e} -> fallback python engine w/ skip", "WARN")
            df = pd.read_csv(
                path,
                sep=sep,
                compression="gzip" if name_lower.endswith(".gz") else None,
                header=None,
                names=header_cols,
                skiprows=skiprows,
                dtype="string",
                low_memory=False,
                engine="python",
                on_bad_lines="skip",
            )
    else:
        # 헤더 발견 실패 → 파일 종류별 기본 스키마로 강제
        kind = _ctd_kind_by_name(path)
        default_cols = CTD_DEFAULT_SCHEMAS.get(kind)
        if not default_cols:
            log(f"Cannot detect header and no default schema known for {path.name}. Attempting raw read ...", "WARN")
            # 최후: comment='#'로 단순 읽기(헤더 미지정) → 첫 행을 헤더로 추정
            try:
                df = pd.read_csv(
                    path,
                    sep=sep,
                    compression="gzip" if name_lower.endswith(".gz") else None,
                    dtype="string",
                    low_memory=False,
                    engine="c",
                    comment="#",
                )
            except ParserError:
                df = pd.read_csv(
                    path,
                    sep=sep,
                    compression="gzip" if name_lower.endswith(".gz") else None,
                    dtype="string",
                    low_memory=False,
                    engine="python",
                    on_bad_lines="skip",
                    comment="#",
                )
        else:
            log(f"Header not found in {path.name}. Using default schema for {kind}", "WARN")
            try:
                df = pd.read_csv(
                    path,
                    sep=sep,
                    compression="gzip" if name_lower.endswith(".gz") else None,
                    header=None,
                    names=default_cols,
                    dtype="string",
                    low_memory=False,
                    engine="c",
                    comment="#",
                )
            except ParserError:
                df = pd.read_csv(
                    path,
                    sep=sep,
                    compression="gzip" if name_lower.endswith(".gz") else None,
                    header=None,
                    names=default_cols,
                    dtype="string",
                    low_memory=False,
                    engine="python",
                    on_bad_lines="skip",
                    comment="#",
                )

    # BOM 제거 및 공백열 제거
    df.columns = [c.lstrip("\ufeff") for c in df.columns]
    df = df.dropna(axis=1, how="all")
    log(f"CTD loaded columns ({path.name}): {list(df.columns)}")
    return df

# -----------------------------
# CTD 파이프라인
# -----------------------------
def load_ctd_vocab() -> Tuple[pd.DataFrame, pd.DataFrame]:
    chem_gz = RAW_CTD_DIR / "CTD_chemicals.tsv.gz"
    dis_gz  = RAW_CTD_DIR / "CTD_diseases.tsv.gz"
    if not chem_gz.exists() or not dis_gz.exists():
        raise FileNotFoundError("CTD vocab missing. Place CTD_chemicals.tsv.gz and CTD_diseases.tsv.gz in data/raw/ctd/")
    chem = smart_read_ctd(chem_gz)
    dis  = smart_read_ctd(dis_gz)
    df_info(chem, "CTD chemicals")
    df_info(dis,  "CTD diseases")
    return chem, dis

def build_ctd_chemical_synonyms(chem: pd.DataFrame) -> pd.DataFrame:
    need = ["ChemicalID", "ChemicalName"]
    miss = [c for c in need if c not in chem.columns]
    if miss:
        raise ValueError(f"CTD_chemicals.tsv.gz missing columns: {miss}")

    synonym_cols = [
        c for c in ["Synonyms", "CTDCuratedSynonyms", "MESHSynonyms"] if c in chem.columns
    ]
    syn_sources = chem[["ChemicalID", "ChemicalName"] + synonym_cols].copy()

    rows = []
    for _, r in syn_sources.iterrows():
        cid = str(r["ChemicalID"]).strip()
        cname = str(r["ChemicalName"]).strip().lower()
        if not cid or cid.lower() == "nan" or not cname or cname == "nan":
            continue
        rows.append((cid, cname, "canonical"))
        for col in synonym_cols:
            values = str(r[col]) if pd.notna(r[col]) else ""
            if not values or values.lower() == "nan":
                continue
            for s in values.split("|"):
                s = s.strip().lower()
                if s and s != cname:
                    rows.append((cid, s, f"synonym:{col.lower()}"))

    syn_df = pd.DataFrame(rows, columns=["ChemicalID", "key", "kind"]).drop_duplicates()
    df_info(syn_df, "CTD chemical synonyms", show_head=5)
    return syn_df

def _detect_human_mask(df: pd.DataFrame) -> pd.Series:
    cols = [c for c in df.columns if c.lower() in ("organism", "organismid", "organism_id", "taxid", "tax_id")]
    if not cols:
        return pd.Series([True] * len(df), index=df.index)  # 필터 불가 시 전부 True
    mask = pd.Series([False] * len(df), index=df.index)
    for c in cols:
        col = df[c].astype("string")
        mask = mask | col.str.contains("homo sapiens", case=False, na=False) | col.str.contains("9606", na=False)
    return mask

def load_ctd_chem_gene(limit_rows: Optional[int] = None) -> pd.DataFrame:
    candidates = [RAW_CTD_DIR / "CTD_chem_gene_ixns.tsv.gz", RAW_CTD_DIR / "CTD_chem_gene_ixns.csv.gz"]
    files = [p for p in candidates if p.exists()]
    if not files:
        files = list(RAW_CTD_DIR.glob("CTD_chem_gene_ixns.*"))
        if not files:
            raise FileNotFoundError("CTD_chem_gene_ixns.* not found in data/raw/ctd/")
    fpath = files[0]
    log(f"Reading {fpath.name}")
    df = smart_read_ctd(fpath)
    if limit_rows:
        df = df.head(limit_rows).copy()
    df_info(df, "CTD chem_gene (raw)", show_head=5)

    # References 컬럼 보정
    if "References" not in df.columns:
        alt = [c for c in df.columns if c.lower() in ("references", "reference", "pubmedids", "pubmed_ids", "publicationids")]
        if alt:
            df.rename(columns={alt[0]: "References"}, inplace=True)

    mask = _detect_human_mask(df)
    if mask.sum() == 0:
        log("No organism columns or no human rows detected; skipping organism filter.", "WARN")
    else:
        df = df[mask].copy()

    # ref_count 계산
    if "References" in df.columns:
        df["ref_count"] = df["References"].fillna("").apply(lambda x: 0 if not x else len(str(x).split("|")))
    else:
        df["ref_count"] = 0

    keep = [c for c in ["ChemicalID", "ChemicalName", "GeneSymbol", "GeneID", "InteractionActions", "ref_count"] if c in df.columns]
    df = df[keep].dropna(subset=[c for c in ["ChemicalID", "GeneSymbol"] if c in keep]).drop_duplicates()
    df_info(df, "CTD chem_gene (human, cleaned)", show_head=5)
    return df

def load_ctd_gene_disease(limit_rows: Optional[int] = None) -> pd.DataFrame:
    f_gd = RAW_CTD_DIR / "CTD_genes_diseases.tsv.gz"
    if not f_gd.exists():
        raise FileNotFoundError("CTD_genes_diseases.tsv.gz missing in data/raw/ctd/")
    df = smart_read_ctd(f_gd)
    if limit_rows:
        df = df.head(limit_rows).copy()
    keep = [c for c in [
        "GeneSymbol", "DiseaseName", "DiseaseID", "DirectEvidence",
        "InferenceScore", "PubMedIDs"
    ] if c in df.columns]
    df = df[keep].dropna(subset=["GeneSymbol", "DiseaseID"]).drop_duplicates()
    direct_raw = df.get("DirectEvidence", pd.Series(index=df.index, dtype="string")).fillna("")
    df["direct_evidence"] = direct_raw.astype(str)
    df["has_direct_ev"] = df["direct_evidence"].str.strip().ne("")
    if "InferenceScore" in df.columns:
        df["inference_score"] = pd.to_numeric(df["InferenceScore"], errors="coerce").fillna(0.0)
    else:
        df["inference_score"] = 0.0
    if "PubMedIDs" in df.columns:
        pubmed_series = df["PubMedIDs"].fillna("").astype(str)
        df["pubmed_ids"] = pubmed_series
        df["pubmed_count"] = pubmed_series.apply(
            lambda s: sum(1 for x in s.split("|") if x.strip())
        )
    else:
        df["pubmed_ids"] = ""
        df["pubmed_count"] = 0
    df = df.drop(columns=[c for c in ["InferenceScore", "DirectEvidence", "PubMedIDs"] if c in df.columns])
    df_info(df, "CTD gene_disease (cleaned)", show_head=5)
    return df

def filter_metabolic_diseases(diseases_df: pd.DataFrame) -> pd.DataFrame:
    need = ["DiseaseID", "DiseaseName", "TreeNumbers"]
    miss = [c for c in need if c not in diseases_df.columns]
    if miss:
        raise ValueError(f"CTD_diseases.tsv.gz missing columns: {miss}")
    d = diseases_df[need].copy()
    d["TreeNumbers"] = d["TreeNumbers"].fillna("")
    mask = d["TreeNumbers"].str.contains(r"\bC18", na=False) | d["TreeNumbers"].str.contains(r"\bC19", na=False)
    out = d[mask].drop_duplicates()
    log(f"Metabolic/Endocrine diseases filtered: {out.shape[0]} rows")
    return out

# -----------------------------
# GWAS (옵션)
# -----------------------------
def load_gwas_assoc(limit_rows: Optional[int] = None) -> Optional[pd.DataFrame]:
    assoc = RAW_GWAS_DIR / "gwas-catalog-associations_ontology-annotated.tsv"
    if not assoc.exists():
        log("GWAS associations TSV not found (optional) -> skipping.", "WARN")
        return None
    df = pd.read_csv(assoc, sep="\t", dtype="string", low_memory=False, nrows=limit_rows)
    gene_cols = [c for c in df.columns if "MAPPED_GENE" in c.upper() or "REPORTED GENE" in c.upper()]
    trait_cols = [c for c in df.columns if "DISEASE/TRAIT" in c.upper() or "MAPPED_TRAIT" in c.upper() or "EFO_TERM" in c.upper()]
    keep = sorted(set(gene_cols + trait_cols))
    if not keep:
        log("GWAS association file has no expected gene/trait columns -> skipping.", "WARN")
        return None
    df = df[keep].dropna(how="all")
    df_info(df, "GWAS associations (subset)", show_head=5)
    return df

# -----------------------------
# FDC nutrient ↔ CTD chemical 매핑
# -----------------------------
def normalize_key(s: str) -> str:
    return re.sub(r"\s+", " ", str(s).strip().lower())

def strip_parentheses(text: str) -> str:
    return re.sub(r"\([^)]*\)", "", text)

def clean_descriptor(text: str) -> str:
    t = text.replace("total", "").replace("by difference", "")
    t = re.sub(r"from .*", "", t)
    return t

def _generate_candidate_keys(nutrient_name: str) -> List[str]:
    raw = str(nutrient_name)
    variants = []
    base = normalize_key(raw)
    variants.append(base)
    no_paren = normalize_key(strip_parentheses(raw))
    variants.append(no_paren)
    simple = normalize_key(re.sub(r"[^a-z0-9]+", " ", raw))
    variants.append(simple)
    cleaned = normalize_key(clean_descriptor(raw))
    variants.append(cleaned)
    for sep in [",", "/", "-", "(" , ")"]:
        if sep in raw:
            parts = [normalize_key(p) for p in raw.split(sep) if p.strip()]
            variants.extend(parts)
    final = []
    for v in variants:
        v = v.strip()
        if v:
            final.append(v)
    return list(dict.fromkeys(final))

def map_fdc_nutrients_to_ctd(
    fn_joined: pd.DataFrame,
    chem_syn: pd.DataFrame,
    overrides_path: Optional[Path] = None,
    fuzzy_threshold: int = 85
) -> pd.DataFrame:
    chem_syn = chem_syn.copy()
    chem_syn["key"] = chem_syn["key"].astype(str).str.strip().str.lower()
    keys = chem_syn["key"].unique().tolist()

    # 수동 오버라이드
    overrides: Dict[str, Tuple[str, str]] = {}
    if overrides_path and overrides_path.exists():
        try:
            ov = pd.read_csv(overrides_path)
            for _, r in ov.iterrows():
                nrm = normalize_key(r["nutrient_name"])
                overrides[nrm] = (str(r.get("chemical_id") or ""), normalize_key(r["ctd_key"]))
            log(f"Loaded overrides: {len(overrides)} rows from {overrides_path}")
        except Exception as e:
            log(f"Failed to read overrides: {e}", "WARN")
    for nrm, (chem_id, key_hint) in MANUAL_NUTRIENT_CTDS.items():
        overrides.setdefault(nrm, (chem_id, key_hint))

    uniq_nutrients = (
        fn_joined[["name", "unit_name"]]
        .drop_duplicates()
        .sort_values("name", kind="stable")
        .reset_index(drop=True)
    )

    rows = []
    missing = []
    for _, r in uniq_nutrients.iterrows():
        nname = str(r["name"])
        norm = normalize_key(nname)
        variants = _generate_candidate_keys(nname)
        if not variants:
            variants = [norm]
        chem_id = None
        matched_key = None
        how = "none"
        fuzzy_score = np.nan

        # 1) overrides
        if norm in overrides:
            cid, k = overrides[norm]
            if cid:
                chem_id = cid; matched_key = k; how = "override_cid"
            elif k:
                cand = chem_syn[chem_syn["key"] == k]
                if len(cand):
                    chem_id = cand["ChemicalID"].iloc[0]; matched_key = k; how = "override_key"
        if chem_id is None:
            for var in variants:
                if var in overrides:
                    cid, k = overrides[var]
                    if cid:
                        chem_id = cid; matched_key = k; how = "override_variant"
                        break
                    elif k:
                        cand = chem_syn[chem_syn["key"] == k]
                        if len(cand):
                            chem_id = cand["ChemicalID"].iloc[0]; matched_key = k; how = "override_variant"
                            break

        # 2) exact
        if chem_id is None:
            for var in variants:
                cand = chem_syn[chem_syn["key"] == var]
                if len(cand):
                    chem_id = cand["ChemicalID"].iloc[0]; matched_key = var; how = "exact"
                    break

        # 3) fuzzy
        if chem_id is None and HAS_RAPIDFUZZ:
            search_name = variants[0]
            result = process.extractOne(search_name, keys, scorer=fuzz.WRatio, score_cutoff=fuzzy_threshold)
            if result:
                matched_key, fuzzy_score = result[0], float(result[1])
                chem_id = chem_syn[chem_syn["key"] == matched_key]["ChemicalID"].iloc[0]
                how = "fuzzy"

        if chem_id is None:
            missing.append(nname)

        rows.append({
            "fdc_nutrient_name": nname,
            "unit_name": r.get("unit_name"),
            "ChemicalID": chem_id,
            "matched_key": matched_key,
            "match_method": how,
            "fuzzy_score": fuzzy_score
        })

    map_df = pd.DataFrame(rows)
    n_total = len(map_df); n_mapped = map_df["ChemicalID"].notna().sum()
    log(f"FDC->CTD mapping: mapped {n_mapped}/{n_total} ({(n_mapped/n_total*100 if n_total else 0):.1f}%)")
    if missing:
        log(f"Unmapped nutrients (first 20): {missing[:20]}", "WARN")
    map_df.sort_values(["match_method", "fuzzy_score"], ascending=[True, False]).to_csv(
        PROC_DIR / "debug_fdc_to_ctd_mapping.csv", index=False, encoding="utf-8"
    )
    return map_df

# -----------------------------
# 진단 도우미
# -----------------------------
def diagnose_fdc(limit_rows: Optional[int] = None):
    try:
        food, nutrient, food_nutrient, food_cat = load_fdc_core(RAW_FDC_DIR, limit_rows=limit_rows, slim=True)
        if "food_category_id" in food.columns:
            vc = food["food_category_id"].astype("Int64").value_counts(dropna=True).head(10)
            log(f"food_category_id top10: {vc.to_dict()}")
        na_amt = food_nutrient["amount"].isna().sum()
        log(f"food_nutrient.amount NaN: {na_amt}/{len(food_nutrient)}")
        if "unit_name" in nutrient.columns:
            units = nutrient["unit_name"].dropna().unique().tolist()[:20]
            log(f"nutrient.unit_name unique(sample): {units}")
    except Exception as e:
        log(f"FDC diagnose failed: {e}", "ERROR")

def diagnose_ctd(limit_rows: Optional[int] = None):
    try:
        chem, dis = load_ctd_vocab()
        # 샘플 동의어 빌드
        syn = build_ctd_chemical_synonyms(chem.head(5000))
        df_info(syn, "CTD chemical synonyms (sample)", show_head=5)
        cgi = load_ctd_chem_gene(limit_rows=limit_rows)
        gd  = load_ctd_gene_disease(limit_rows=limit_rows)
        dmeta = filter_metabolic_diseases(dis)
        df_info(dmeta, "CTD diseases (metabolic subset)", show_head=5)
    except Exception as e:
        log(f"CTD diagnose failed: {e}", "ERROR")

def diagnose_gwas(limit_rows: Optional[int] = None):
    try:
        gwas = load_gwas_assoc(limit_rows=limit_rows)
        if gwas is not None:
            df_info(gwas, "GWAS (diag)", show_head=5)
    except Exception as e:
        log(f"GWAS diagnose failed: {e}", "ERROR")

# -----------------------------
# 메인
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diagnose", choices=["fdc", "ctd", "gwas", "all"], help="진단 모드 실행")
    ap.add_argument("--limit-rows", type=int, default=None, help="디버그용 상위 N행만 읽기")
    ap.add_argument("--slim", type=str, default="true", help="조인 결과 슬림 컬럼만 저장(true/false)")
    ap.add_argument("--fuzzy-threshold", type=int, default=92, help="퍼지 매칭 점수 하한(0~100)")
    args = ap.parse_args()

    slim = str(args.slim).strip().lower() in ("1", "true", "yes", "y")

    log(f"Python: {sys.version.split()[0]}")
    log(f"pandas: {pd.__version__}")
    log(f"Working dir: {Path('.').resolve()}")
    ensure_dirs()

    if args.diagnose:
        if args.diagnose in ("fdc", "all"):
            log("=== DIAG: FDC ===")
            diagnose_fdc(limit_rows=args.limit_rows)
        if args.diagnose in ("ctd", "all"):
            log("=== DIAG: CTD ===")
            diagnose_ctd(limit_rows=args.limit_rows)
        if args.diagnose in ("gwas", "all"):
            log("=== DIAG: GWAS ===")
            diagnose_gwas(limit_rows=args.limit_rows)
        return

    # ----- FDC -----
    try:
        log("Loading FDC ...")
        food, nutrient, food_nutrient, food_cat = load_fdc_core(RAW_FDC_DIR, limit_rows=args.limit_rows, slim=slim)
        fn_joined = build_food_nutrient_table(food, nutrient, food_nutrient, slim=slim)
        # 저장
        food_cols = ["fdc_id", "description", "data_type"]
        if "food_category_id" in food.columns:
            food_cols.append("food_category_id")
        if "publication_date" in food.columns:
            food_cols.append("publication_date")
        food[food_cols].to_parquet(PROC_DIR / "fdc_foods.parquet", index=False)
        fn_joined.to_parquet(PROC_DIR / "fdc_food_nutrient.parquet", index=False)
        log("Saved: fdc_foods.parquet, fdc_food_nutrient.parquet")
    except Exception as e:
        log(f"FDC preprocessing failed: {e}", "ERROR")
        raise

    # ----- CTD -----
    try:
        log("Loading CTD vocab ...")
        chem_vocab, dis_vocab = load_ctd_vocab()
        chem_syn = build_ctd_chemical_synonyms(chem_vocab)
        chem_syn.to_parquet(PROC_DIR / "ctd_chem_synonyms.parquet", index=False)
        log("Saved: ctd_chem_synonyms.parquet")

        log("Loading CTD chem_gene ...")
        cgi = load_ctd_chem_gene(limit_rows=args.limit_rows)
        cgi.to_parquet(PROC_DIR / "ctd_chem_gene_human.parquet", index=False)
        log("Saved: ctd_chem_gene_human.parquet")

        log("Loading CTD gene_disease ...")
        gd = load_ctd_gene_disease(limit_rows=args.limit_rows)
        gd.to_parquet(PROC_DIR / "ctd_gene_disease.parquet", index=False)
        log("Saved: ctd_gene_disease.parquet")

        log("Filtering metabolic/endocrine diseases (MeSH C18/C19) ...")
        dmeta = filter_metabolic_diseases(dis_vocab)
        dmeta.to_parquet(PROC_DIR / "ctd_diseases_metabolic.parquet", index=False)
        log("Saved: ctd_diseases_metabolic.parquet")
    except Exception as e:
        log(f"CTD preprocessing failed: {e}", "ERROR")
        raise

    # ----- 매핑 -----
    try:
        log("Mapping FDC nutrients to CTD chemicals ...")
        map_df = map_fdc_nutrients_to_ctd(
            fn_joined, chem_syn, overrides_path=OVERRIDES_PATH, fuzzy_threshold=int(args.fuzzy_threshold)
        )
        map_df.to_parquet(PROC_DIR / "fdc_to_ctd_chemical_map.parquet", index=False)
        log("Saved: fdc_to_ctd_chemical_map.parquet, debug_fdc_to_ctd_mapping.csv")
    except Exception as e:
        log(f"Mapping failed: {e}", "ERROR")
        raise

    # ----- GWAS (옵션) -----
    try:
        log("Loading GWAS associations (optional) ...")
        gwas = load_gwas_assoc(limit_rows=args.limit_rows)
        if gwas is not None:
            gwas.to_parquet(PROC_DIR / "gwas_assoc_subset.parquet", index=False)
            log("Saved: gwas_assoc_subset.parquet")
        else:
            log("GWAS not available. Skipped.")
    except Exception as e:
        log(f"GWAS preprocessing failed: {e}", "ERROR")

    log("All preprocessing steps completed successfully. ✅")

if __name__ == "__main__":
    main()
