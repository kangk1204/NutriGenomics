# -*- coding: utf-8 -*-
"""
KFRI NutriGenomics Food-Health Database Analyzer
데이터 구조 분석 결과를 반영한 실무형 버전

수정사항:
- 검증 기준을 데이터 현실에 맞게 조정
- 화학물질 ID 정규화 로직 개선
- 경로 품질 점수 계산 방식 수정
- 실용적인 과학적 기준 적용
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo
from pathlib import Path
import networkx as nx
from collections import defaultdict, Counter
import html
import json
import pickle
import time
import shutil
import re
import warnings
from typing import Dict, List, Tuple, Optional, Set
from tqdm import tqdm
import math

METABOLIC_KEYWORDS = [
    "metabolic", "obesity", "diabetes", "glucose", "lipid", "insulin",
    "triglyceride", "hdl", "ldl", "cholesterol", "nafld", "nash",
    "fatty", "steatosis", "bmi", "hyperglycemia", "dyslipidemia"
]

DEFAULT_DV = {
    "Vitamin C, total ascorbic acid": ("MG", 90.0),
    "Vitamin D (D2 + D3)": ("UG", 20.0),
    "Vitamin E (alpha-tocopherol)": ("MG", 15.0),
    "Vitamin K (phylloquinone)": ("UG", 120.0),
    "Vitamin A, RAE": ("UG", 900.0),
    "Folate, total": ("UG", 400.0),
    "Choline, total": ("MG", 550.0),
    "Niacin": ("MG", 16.0),
    "Riboflavin": ("MG", 1.3),
    "Thiamin": ("MG", 1.2),
    "Vitamin B-6": ("MG", 1.7),
    "Vitamin B-12": ("UG", 2.4),
    "Calcium, Ca": ("MG", 1300.0),
    "Iron, Fe": ("MG", 18.0),
    "Magnesium, Mg": ("MG", 420.0),
    "Phosphorus, P": ("MG", 1250.0),
    "Potassium, K": ("MG", 4700.0),
    "Sodium, Na": ("MG", 2300.0),
    "Zinc, Zn": ("MG", 11.0),
    "Selenium, Se": ("UG", 55.0),
    "Fiber, total dietary": ("G", 28.0),
    "Protein": ("G", 50.0),
    "Total lipid (fat)": ("G", 78.0)
}

PRO_RISK_GENES: Set[str] = {
    "TNF", "IL6", "IL1B", "CRP", "CCL2", "CCL5", "CXCL8", "NFKB1", "RELA",
    "JUN", "FOS", "MAPK8", "MAPK14", "CASP3", "CASP8", "BAX", "DDIT3"
}

ANTI_RISK_GENES: Set[str] = {
    "CAT", "SOD1", "SOD2", "GPX1", "GPX3", "GSR", "HMOX1", "NQO1", "IL10",
    "PPARG", "PPARGC1A", "ADIPOQ", "INSR", "SLC2A4", "AKT1", "PRKAA1", "PRKAA2"
}

DIRECT_EVIDENCE_PRIORS = {
    "therapeutic": 1.0,
    "dosage": 0.9,
    "marker/mechanism": 0.8,
    "activity": 0.7,
    "expression": 0.6
}

warnings.filterwarnings('ignore')


def nan_iqr(values) -> float:
    """Compute interquartile range ignoring NaNs."""
    if values is None:
        return 0.0
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return 0.0
    q75 = np.nanpercentile(arr, 75)
    q25 = np.nanpercentile(arr, 25)
    return float(q75 - q25)

MANUAL_NUTRIENT_OVERRIDES = {
    "Vitamin C, total ascorbic acid": ("MESH:D014786", 0.92),
    "Vitamin A, RAE": ("MESH:D014801", 0.9),
    "Vitamin D (D2 + D3)": ("MESH:D014807", 0.9),
    "Vitamin K (phylloquinone)": ("MESH:D014810", 0.9),
    "Vitamin E (alpha-tocopherol)": ("MESH:D014805", 0.88),
    "Choline, total": ("MESH:D002813", 0.85),
    "Protein": ("MESH:D011506", 0.75),
    "Total lipid (fat)": ("MESH:D010999", 0.75),
    "Fiber, total dietary": ("MESH:D005345", 0.8),
    "Folate, total": ("MESH:D005492", 0.85),
    "Iron, Fe": ("MESH:D007501", 0.9),
    "Magnesium, Mg": ("MESH:D008274", 0.88),
    "Potassium, K": ("MESH:D011188", 0.88),
    "Calcium, Ca": ("MESH:D002115", 0.92),
    "Zinc, Zn": ("MESH:D015043", 0.9)
}

EXCLUDED_FOOD_KEYWORDS = [
    "water enhancer",
    "drink enhancer",
    "liquid water",
    "flavor drops",
    "energy shot",
    "pre-workout",
    "supplement mix",
    "powder drink mix",
    "dietary supplement"
]

class KFRINutriGenomicsAnalyzer:
    def __init__(self, data_dir: str = "data/processed", cache_dir: str = "cache_nutriomics"):
        self.data_dir = Path(data_dir)
        self.cache_dir = Path(cache_dir)
        
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        self.data = {}
        self.validated_connections = {}
        self.evidence_scores = {}
        self.norm_stats: Dict[str, Dict[str, float]] = {}
        self.gwas_df: Optional[pd.DataFrame] = None
        self.gwas_support_map: Dict[str, float] = {}
        self.metabolic_keywords = METABOLIC_KEYWORDS

        # 데이터 현실에 맞춘 실용적 기준
        self.quality_thresholds = {
            'min_literature_support': 1,        # 최소 문헌 지원 수
            'min_inference_score': 2.5,         # 최소 추론 점수
            'min_mapping_confidence': 0.5,      # 매핑 신뢰도 최소치
            'min_direct_evidence_ratio': 0.0,
            'min_dose_threshold_mg': 0.5,       # 0.5 mg/100g 이상만 기본 통과
            'chemical_id_confidence': 0.3,
            'pathway_significance': 0.1,
            'min_pathway_quality': 0.6
        }
        self.min_dv_ratio = 0.01  # 최소 1% 일일권장량 기여
        self.min_micro_mg = 0.001  # 1µg
        self.max_nutrients_per_food = 5
        self.max_genes_per_nutrient = 3
        self.max_diseases_per_gene = 3
        self.max_pathways = 6000
        self.dv_lookup_mg: Dict[str, float] = {}
        for nutrient, (unit, value) in DEFAULT_DV.items():
            if unit == 'UG':
                mg_value = value / 1000.0
            elif unit == 'G':
                mg_value = value * 1000.0
            else:
                mg_value = value
            self.dv_lookup_mg[nutrient.lower()] = mg_value

        self.score_weights = {
            'dose': 0.20,
            'mapping': 0.15,
            'literature': 0.15,
            'inference': 0.20,
            'support': 0.20,
            'direction': 0.10
        }
        total_weight = sum(self.score_weights.values())
        if not math.isclose(total_weight, 1.0, rel_tol=1e-6):
            self.score_weights = {k: v / total_weight for k, v in self.score_weights.items()}

        self.min_direction_alignment = 0.15
        self.gene_role_map = self._build_gene_role_map()
        self.min_pubmed_support = 2
        self.food_exclusion_regex = re.compile("|".join(re.escape(k) for k in EXCLUDED_FOOD_KEYWORDS), re.IGNORECASE)
        
        # 간소화된 생물학적 타당성 평가
        self.nutrient_categories = {
            'vitamins': ['vitamin', 'tocopherol', 'retinol', 'thiamin', 'riboflavin', 'niacin', 'ascorbic', 'biotin'],
            'minerals': ['calcium', 'iron', 'zinc', 'magnesium', 'phosphorus', 'potassium', 'sodium'],
            'fatty_acids': ['fatty', 'oleic', 'linoleic', 'palmitic', 'epa', 'dha', 'mufa', 'pufa'],
            'amino_acids': ['alanine', 'arginine', 'histidine', 'glutamic', 'glycine'],
            'metabolites': ['cholesterol', 'caffeine', 'glucose', 'choline', 'glutathione']
        }
        nutrient_terms = sorted({kw for kws in self.nutrient_categories.values() for kw in kws})
        self._nutrient_keyword_regex = re.compile(r"(?:%s)" % "|".join(re.escape(k) for k in nutrient_terms), re.IGNORECASE)

    def load_data_efficiently(self):
        """효율적 데이터 로딩"""
        print("Loading data for KFRI NutriGenomics analysis...")
        
        files_config = {
            'foods': {'file': 'fdc_foods.parquet'},
            'food_nutrients': {
                'file': 'fdc_food_nutrient.parquet',
                'sample': 400000,
                'columns': ['description', 'name', 'amount_mg_per_100g'],
                'stats_column': 'amount_mg_per_100g'
            },
            'chem_synonyms': {'file': 'ctd_chem_synonyms.parquet'},
            'chem_genes': {
                'file': 'ctd_chem_gene_human.parquet',
                'sample': 250000,
                'columns': ['ChemicalID', 'GeneSymbol', 'ref_count', 'InteractionActions'],
                'stats_column': 'ref_count'
            },
            'gene_diseases': {
                'file': 'ctd_gene_disease.parquet',
                'sample': 350000,
                'columns': [
                    'GeneSymbol', 'DiseaseID', 'DiseaseName',
                    'inference_score', 'has_direct_ev', 'direct_evidence', 'pubmed_count'
                ],
                'stats_column': 'inference_score'
            },
            'metabolic_diseases': {'file': 'ctd_diseases_metabolic.parquet'},
            'nutrient_mapping': {'file': 'fdc_to_ctd_chemical_map.parquet'}
        }
        self.dataset_paths: Dict[str, Path] = {}
        self.files_config = files_config
        
        for key, config in files_config.items():
            filepath = self.data_dir / config['file']
            if filepath.exists():
                self.dataset_paths[key] = filepath
                df = self._read_with_optional_columns(filepath, config.get('columns'))

                sample_size = config.get('sample')
                if sample_size and len(df) > sample_size:
                    df = df.sample(sample_size, random_state=42)

                self.data[key] = df
                print(f"Loaded {key}: {len(df):,} rows")
            else:
                print(f"[WARN] {filepath} not found – skipping {key}")

        self._prepare_normalization_stats()
        self._maybe_load_gwas()

    @staticmethod
    def _read_with_optional_columns(filepath: Path, columns: Optional[List[str]]):
        if filepath.suffix == '.parquet':
            try:
                return pd.read_parquet(filepath, columns=columns)
            except (KeyError, ValueError):
                return pd.read_parquet(filepath)
        usecols = columns if columns else None
        return pd.read_csv(filepath, usecols=usecols)

    def _get_full_series(self, dataset_key: str, column: str) -> Optional[pd.Series]:
        config = self.files_config.get(dataset_key, {}) if hasattr(self, "files_config") else {}
        sampled = bool(config.get('sample'))
        df = self.data.get(dataset_key)
        if df is not None and column in df.columns and not sampled:
            return df[column]
        path = self.dataset_paths.get(dataset_key)
        if not path:
            return None
        try:
            if path.suffix == '.parquet':
                series = pd.read_parquet(path, columns=[column])[column]
            else:
                series = pd.read_csv(path, usecols=[column])[column]
            return series
        except Exception as exc:
            print(f"[WARN] Failed to load column '{column}' for {dataset_key}: {exc}")
            return None

    def _prepare_normalization_stats(self):
        """점수 계산에 사용할 분포 통계 준비"""
        stats_map: Dict[str, Dict[str, float]] = {}

        amount_series = self._get_full_series('food_nutrients', 'amount_mg_per_100g')
        if amount_series is not None:
            amount = pd.to_numeric(amount_series, errors='coerce').clip(lower=0)
            log_amount = np.log1p(amount.dropna())
            if not log_amount.empty:
                median = float(np.nanmedian(log_amount))
                iqr = nan_iqr(log_amount)
                scale = max(iqr / 1.35, 1e-3)
                stats_map['nutrient_amount'] = {'median': median, 'scale': scale}

        ref_series = self._get_full_series('chem_genes', 'ref_count')
        if ref_series is not None:
            ref_count = pd.to_numeric(ref_series, errors='coerce').clip(lower=0)
            log_ref = np.log1p(ref_count.dropna())
            if not log_ref.empty:
                median = float(np.nanmedian(log_ref))
                iqr = nan_iqr(log_ref)
                scale = max(iqr / 1.35, 1e-3)
                stats_map['reference_count'] = {'median': median, 'scale': scale}
                q75_raw = float(np.nanpercentile(ref_count.dropna(), 75))
                self.quality_thresholds['min_literature_support'] = max(
                    self.quality_thresholds['min_literature_support'],
                    int(max(1, round(q75_raw)))
                )

        inf_series = self._get_full_series('gene_diseases', 'inference_score')
        if inf_series is not None:
            inf = pd.to_numeric(inf_series, errors='coerce').clip(lower=0)
            inf = inf.dropna()
            if not inf.empty:
                median = float(np.nanmedian(inf))
                iqr = nan_iqr(inf)
                scale = max(iqr / 1.35, 1e-3)
                stats_map['inference_score'] = {'median': median, 'scale': scale}
                q80 = float(np.nanpercentile(inf, 80))
                self.quality_thresholds['min_inference_score'] = max(
                    self.quality_thresholds['min_inference_score'],
                    q80
                )

        self.norm_stats = stats_map

    def _sigmoid_scale(self, value: float, key: str, use_log: bool=False) -> float:
        stats_conf = self.norm_stats.get(key)
        if stats_conf is None or value is None or (isinstance(value, float) and np.isnan(value)):
            return 0.0
        x = float(value)
        if use_log:
            x = np.log1p(max(x, 0.0))
        mid = stats_conf['median']
        scale = stats_conf['scale']
        if scale <= 0:
            return 0.0
        return float(1.0 / (1.0 + np.exp(-(x - mid) / scale)))

    def _build_gene_role_map(self) -> Dict[str, str]:
        role = {}
        for g in PRO_RISK_GENES:
            role[g.upper()] = 'pro'
        for g in ANTI_RISK_GENES:
            role[g.upper()] = 'anti'
        override_path = Path('config/gene_roles.csv')
        if override_path.exists():
            try:
                override_df = pd.read_csv(override_path)
                if {'gene', 'role'}.issubset({c.lower() for c in override_df.columns}):
                    gene_col = next(c for c in override_df.columns if c.lower() == 'gene')
                    role_col = next(c for c in override_df.columns if c.lower() == 'role')
                    for _, row in override_df.iterrows():
                        gene = str(row[gene_col]).strip().upper()
                        r = str(row[role_col]).strip().lower()
                        if gene and r in ('pro', 'anti'):
                            role[gene] = r
            except Exception as exc:
                print(f"Failed to read gene role overrides: {exc}")
        return role

    def _maybe_load_gwas(self) -> None:
        path = self.data_dir / 'gwas_assoc_subset.parquet'
        if not path.exists():
            print("GWAS subset parquet not found; GWAS support disabled")
            self.gwas_df = None
            return
        try:
            self.gwas_df = pd.read_parquet(path)
            print(f"GWAS evidence loaded: {self.gwas_df.shape}")
        except Exception as exc:
            print(f"Failed to load GWAS parquet: {exc}")
            self.gwas_df = None

    def _build_gwas_support_map(self, genes: Set[str]) -> Dict[str, float]:
        if self.gwas_df is None or self.gwas_df.empty:
            return {}
        genes = {str(g).upper() for g in genes if isinstance(g, str) and g.strip()}
        if not genes:
            return {}
        gene_cols = [c for c in self.gwas_df.columns if "GENE" in c.upper()]
        trait_cols = [c for c in self.gwas_df.columns if ("TRAIT" in c.upper()) or ("EFO" in c.upper())]
        if not gene_cols or not trait_cols:
            return {}
        trait_text = self.gwas_df[trait_cols].astype(str).agg(" ".join, axis=1).str.lower()
        kw_pattern = "|".join(re.escape(k.lower()) for k in self.metabolic_keywords)
        trait_mask = trait_text.str.contains(kw_pattern, na=False)
        subset = self.gwas_df[trait_mask].copy()
        if subset.empty:
            return {}
        gene_text = subset[gene_cols].astype(str).agg(" ".join, axis=1)
        support_map: Dict[str, float] = {}
        for gene in genes:
            pattern = re.compile(rf"\b{re.escape(gene)}\b", flags=re.IGNORECASE)
            hits = gene_text.str.contains(pattern, na=False)
            if hits.any():
                support_map[gene] = 1.3
        return support_map

    def _ensure_gwas_support_map(self, genes: Set[str]) -> None:
        if not genes:
            return
        upper_genes = {str(g).upper() for g in genes if isinstance(g, str) and g.strip()}
        missing = upper_genes - set(self.gwas_support_map.keys())
        if not missing:
            return
        new_map = self._build_gwas_support_map(missing)
        self.gwas_support_map.update(new_map)

    def _dose_ratio_component(self, nutrient: str, amount_mg: float) -> float:
        if nutrient is None or amount_mg is None or amount_mg <= 0:
            return 0.0
        dv_mg = self.dv_lookup_mg.get(str(nutrient).lower())
        if not dv_mg or dv_mg <= 0:
            return 0.0
        ratio = amount_mg / dv_mg
        return float(max(0.0, min(1.0, 1.0 - math.exp(-ratio))))

    @staticmethod
    def _convert_amount_to_mg(amount: float, unit: str) -> float:
        if amount is None or unit is None:
            return 0.0
        unit = str(unit).upper()
        if unit == 'MG':
            return float(amount)
        if unit == 'UG':
            return float(amount) / 1000.0
        if unit == 'G':
            return float(amount) * 1000.0
        return 0.0

    def _parse_direction_label(self, action_text: Optional[str]) -> str:
        if not action_text or not isinstance(action_text, str):
            return 'unknown'
        text = action_text.lower()
        up_keywords = ['increase', 'activate', 'upregulat', 'induce']
        down_keywords = ['decrease', 'downregulat', 'inhibit', 'suppress', 'reduce']
        up = any(k in text for k in up_keywords)
        down = any(k in text for k in down_keywords)
        if up and down:
            return 'mixed'
        if up:
            return 'up'
        if down:
            return 'down'
        return 'unknown'

    def _compute_direction_alignment(self, gene_symbol: str, direction: str) -> float:
        if not gene_symbol:
            return 0.5
        role = self.gene_role_map.get(str(gene_symbol).upper())
        if role is None:
            return 0.5
        direction = direction or 'unknown'
        if direction == 'unknown':
            return 0.5
        if role == 'pro':
            if direction == 'down':
                return 1.0
            if direction == 'up':
                return 0.0
            if direction == 'mixed':
                return 0.5
        elif role == 'anti':
            if direction == 'up':
                return 1.0
            if direction == 'down':
                return 0.0
            if direction == 'mixed':
                return 0.5
        return 0.5
    
    def normalize_chemical_id(self, chem_id) -> str:
        """개선된 화학물질 ID 정규화"""
        if pd.isna(chem_id):
            return None
        
        chem_id = str(chem_id).strip()
        
        # 알려진 접두사 제거 (분석 결과 기반)
        prefixes = ['MESH:', 'mesh:', 'CTD:', 'ctd:', 'PUBCHEM:', 'pubchem:']
        for prefix in prefixes:
            if chem_id.startswith(prefix):
                chem_id = chem_id[len(prefix):]
        
        return chem_id.upper()

    def _normalize_chem_id_series(self, series: pd.Series) -> pd.Series:
        if series.empty:
            return series.astype("string")
        normalized = series.astype("string").str.strip()
        normalized = normalized.str.replace(r'^(?:MESH:|CTD:|PUBCHEM:)', '', regex=True, flags=re.IGNORECASE)
        normalized = normalized.str.upper()
        return normalized.where(normalized != "", None)
    
    def build_chemical_mapping(self) -> Dict[str, str]:
        """정규화된 화학물질 매핑"""
        print("Building chemical mapping...")
        
        if not all(k in self.data for k in ['nutrient_mapping', 'chem_genes']):
            return {}
        
        mapping = self.data['nutrient_mapping'].copy()
        chem_genes = self.data['chem_genes'].copy()
        
        mapping_ids = mapping['ChemicalID'] if 'ChemicalID' in mapping.columns else pd.Series(index=mapping.index, dtype="string")
        chem_ids = chem_genes['ChemicalID'] if 'ChemicalID' in chem_genes.columns else pd.Series(index=chem_genes.index, dtype="string")
        mapping['normalized_chem_id'] = self._normalize_chem_id_series(mapping_ids)
        chem_genes['normalized_chem_id'] = self._normalize_chem_id_series(chem_ids)
        
        successful_mapping = mapping[mapping['normalized_chem_id'].notna()].copy()
        print(f"Found {len(successful_mapping)} mapped nutrients in mapping table")
        
        chem_gene_ids = set(chem_genes['normalized_chem_id'].dropna())
        successful_mapping = successful_mapping[successful_mapping['normalized_chem_id'].isin(chem_gene_ids)].copy()
        
        if successful_mapping.empty:
            print("No nutrient-chemical pairs passed normalized ID validation")
            return {}
        
        method_series = successful_mapping['match_method'].fillna('unknown').astype(str).str.lower() if 'match_method' in successful_mapping.columns else pd.Series('unknown', index=successful_mapping.index)
        fuzzy_scores = pd.to_numeric(successful_mapping['fuzzy_score'], errors='coerce').fillna(0.0) if 'fuzzy_score' in successful_mapping.columns else pd.Series(0.0, index=successful_mapping.index)
        nutrient_is_known = successful_mapping['fdc_nutrient_name'].astype(str).apply(self._is_known_nutrient).astype(float)

        method_lookup = {
            'exact': 0.95,
            'fuzzy': 0.8,
            'synonym': 0.7,
            'manual': 0.9,
            'bridge': 0.65,
            'unknown': 0.55
        }

        raw = (
            0.55 * method_series.map(lambda m: method_lookup.get(m, 0.55)).astype(float) +
            0.25 * np.clip(fuzzy_scores / 100.0, 0.0, 1.0) +
            0.20 * nutrient_is_known
        )
        successful_mapping['confidence'] = (
            1.0 / (1.0 + np.exp(-(raw - 0.6) / 0.15))
        ).clip(0.0, 1.0)
        
        deduped = (
            successful_mapping
            .sort_values('confidence', ascending=False)
            .drop_duplicates(subset=['fdc_nutrient_name'])
            .loc[:, ['fdc_nutrient_name', 'normalized_chem_id', 'confidence', 'ChemicalID']]
            .rename(columns={
                'fdc_nutrient_name': 'nutrient',
                'normalized_chem_id': 'chemical_id',
                'ChemicalID': 'original_id'
            })
        )
        
        validated_mapping = deduped.set_index('nutrient').to_dict(orient='index')
        validated_mapping = self._inject_manual_mappings(validated_mapping)
        print(f"Successfully validated {len(validated_mapping)} chemical mappings")
        return validated_mapping
    
    def _inject_manual_mappings(self, mapping: Dict[str, Dict[str, object]]) -> Dict[str, Dict[str, object]]:
        existing_lower = {k.lower(): k for k in mapping.keys()}
        for nutrient, (chem_id, confidence) in MANUAL_NUTRIENT_OVERRIDES.items():
            key_lower = nutrient.lower()
            if key_lower in existing_lower:
                continue
            mapping[nutrient] = {
                'chemical_id': chem_id.replace("MESH:", ""),
                'original_id': chem_id,
                'confidence': confidence
            }
        return mapping
    
    def _is_known_nutrient(self, nutrient_name: str) -> bool:
        """알려진 영양소인지 확인"""
        if not nutrient_name:
            return False
        return bool(self._nutrient_keyword_regex.search(str(nutrient_name)))

    def _is_excluded_food(self, food_name: str) -> bool:
        if not food_name:
            return False
        return bool(self.food_exclusion_regex.search(food_name))
    
    def build_connections_efficiently(self):
        """효율적인 연결 구축"""
        print("Building connections with validated logic...")
        
        # 화학물질 매핑 검증
        validated_mappings = self.build_chemical_mapping()
        
        self.validated_connections = {
            'food_to_nutrient': defaultdict(list),
            'nutrient_to_gene': defaultdict(list),
            'gene_to_disease': defaultdict(list),
            'complete_pathways': defaultdict(list)
        }
        
        # Step 1: Food → Nutrient (현실적 기준)
        print("Building Food → Nutrient connections...")
        if 'food_nutrients' in self.data:
            fn = self.data['food_nutrients']
            
            # 최소 임계값만 적용
            amount_series = pd.to_numeric(fn.get('amount_mg_per_100g', pd.Series(0, index=fn.index)), errors='coerce')
            nutrient_names = fn.get('name', pd.Series('Unknown', index=fn.index)).astype(str)
            dv_values = nutrient_names.str.lower().map(self.dv_lookup_mg).replace({0: np.nan})
            ratio_series = amount_series.divide(dv_values).replace([np.inf, -np.inf], np.nan)
            dose_mask = amount_series >= self.quality_thresholds['min_dose_threshold_mg']
            ratio_mask = ratio_series >= self.min_dv_ratio
            micro_mask = amount_series >= self.min_micro_mg
            selection_mask = amount_series.notna() & ((dose_mask) | (ratio_mask.fillna(False) & micro_mask))
            quality_fn = fn[selection_mask].copy()
            
            if not quality_fn.empty:
                foods = quality_fn.get('description', pd.Series('Unknown', index=quality_fn.index))
                nutrients = quality_fn.get('name', pd.Series('Unknown', index=quality_fn.index))
                prepared = pd.DataFrame({
                    'food': foods.astype(str).str.slice(0, 40),
                    'nutrient': nutrients.astype(str),
                    'amount': amount_series.loc[quality_fn.index].fillna(0.0),
                    'dv_ratio': ratio_series.loc[quality_fn.index].fillna(0.0)
                })
                
                for food, grp in prepared.groupby('food'):
                    if self._is_excluded_food(food):
                        continue
                    records = [
                        {'nutrient': nutrient, 'amount': float(amount), 'dv_ratio': max(0.0, float(dv_ratio))}
                        for nutrient, amount, dv_ratio in zip(grp['nutrient'], grp['amount'], grp['dv_ratio'])
                    ]
                    if records:
                        records.sort(key=lambda r: (r['dv_ratio'], r['amount']), reverse=True)
                        self.validated_connections['food_to_nutrient'][food] = records[:self.max_nutrients_per_food]
        
        print(f"Connected {len(self.validated_connections['food_to_nutrient'])} foods to nutrients")
        
        # Step 2: Nutrient → Gene (검증된 매핑 활용)
        print("Building Nutrient → Gene connections...")
        if validated_mappings and 'chem_genes' in self.data:
            cg = self.data['chem_genes'].copy()
            chem_ids = cg['ChemicalID'] if 'ChemicalID' in cg.columns else pd.Series(index=cg.index, dtype="string")
            cg['normalized_chem_id'] = self._normalize_chem_id_series(chem_ids)
            cg['ref_count_numeric'] = pd.to_numeric(cg.get('ref_count', 0), errors='coerce').fillna(0.0)
            cg = cg[cg['ref_count_numeric'] >= self.quality_thresholds['min_literature_support']]
            
            vm_df = (
                pd.DataFrame.from_dict(validated_mappings, orient='index')
                .reset_index()
                .rename(columns={'index': 'nutrient'})
            )
            
            merged = (
                cg.merge(vm_df, left_on='normalized_chem_id', right_on='chemical_id', how='inner')
                .dropna(subset=['GeneSymbol'])
            )
            
            for nutrient, grp in merged.groupby('nutrient'):
                self._ensure_gwas_support_map(set(grp['GeneSymbol']))
                genes = [
                    {
                        'gene': str(gene),
                        'chemical_id': chem_id,
                        'ref_count': float(ref_count),
                        'confidence': conf,
                        'interaction_actions': str(action),
                        'direction': self._parse_direction_label(action),
                        'gwas_support': float(self.gwas_support_map.get(str(gene).upper(), 1.0))
                    }
                    for gene, chem_id, ref_count, confidence, action in zip(
                        grp['GeneSymbol'],
                        grp['chemical_id'],
                        grp['ref_count_numeric'],
                        grp['confidence'],
                        grp['InteractionActions'] if 'InteractionActions' in grp.columns else pd.Series('', index=grp.index)
                    )
                    if str(gene) and (conf := float(confidence)) >= self.quality_thresholds['min_mapping_confidence']
                ]
                if genes:
                    genes.sort(key=lambda g: (g['confidence'], g['ref_count'], g.get('gwas_support', 1.0)), reverse=True)
                    self.validated_connections['nutrient_to_gene'][nutrient] = genes
            
            connection_count = sum(len(v) for v in self.validated_connections['nutrient_to_gene'].values())
            print(f"Connected {len(self.validated_connections['nutrient_to_gene'])} nutrients to {connection_count} gene relationships")
        else:
            print("Skipping Nutrient → Gene connections (insufficient mapping or chem gene data)")
        
        # Step 3: Gene → Disease (완화된 기준)
        print("Building Gene → Disease connections...")
        if 'gene_diseases' in self.data:
            gd = self.data['gene_diseases'].copy()
            gd['inference_score'] = pd.to_numeric(gd.get('inference_score', 0), errors='coerce').fillna(0.0)
            gd['has_direct_ev'] = gd.get('has_direct_ev', False).fillna(False).astype(bool)
            gd['GeneSymbol'] = gd.get('GeneSymbol', pd.Series('', index=gd.index)).astype(str)
            gd['DiseaseName'] = gd.get('DiseaseName', pd.Series('Unknown', index=gd.index)).astype(str)
            gd['DiseaseID'] = gd.get('DiseaseID', pd.Series('', index=gd.index))
            
            quality_gd = gd[gd['inference_score'] >= self.quality_thresholds['min_inference_score']]
            aggregated = (
                quality_gd
                .groupby(['GeneSymbol', 'DiseaseID', 'DiseaseName'], dropna=False)
                .agg(
                    has_direct_ev=('has_direct_ev', 'max'),
                    direct_evidence=('direct_evidence', lambda x: "|".join(sorted({str(v).strip() for v in x if isinstance(v, str) and v.strip()}))),
                    inference_score=('inference_score', 'max'),
                    evidence_count=('GeneSymbol', 'size'),
                    pubmed_count=('pubmed_count', 'sum')
                )
                .reset_index()
            )

            metabolic = self.data.get('metabolic_diseases')
            if metabolic is not None and 'DiseaseID' in metabolic.columns:
                meta_ids = set(metabolic['DiseaseID'].dropna().astype(str))
                if meta_ids:
                    aggregated['DiseaseID'] = aggregated['DiseaseID'].astype(str)
                    aggregated = aggregated[aggregated['DiseaseID'].isin(meta_ids)]
                    if aggregated.empty:
                        print("No gene-disease rows remain after metabolic filter")
                        self.validated_connections['gene_to_disease'] = defaultdict(list)
                        return
                    aggregated['DiseaseID'] = aggregated['DiseaseID'].astype(str)
            
            filtered = aggregated[
                (aggregated['inference_score'] >= self.quality_thresholds['min_inference_score']) &
                ((aggregated['evidence_count'] >= 2) | (aggregated['pubmed_count'] >= self.min_pubmed_support))
            ]
            for gene, grp in filtered.groupby('GeneSymbol'):
                if not gene:
                    continue
                diseases = []
                for _, row in grp.iterrows():
                    disease_name = str(row['DiseaseName'])
                    disease_id = row['DiseaseID']
                    codes = [c.strip().lower() for c in str(row.get('direct_evidence', '')).split('|') if c.strip()]
                    diseases.append({
                        'disease': disease_name[:40],
                        'disease_id': disease_id,
                        'has_direct_evidence': bool(row['has_direct_ev']),
                        'inference_score': float(row['inference_score']),
                        'evidence_count': int(row['evidence_count']),
                        'pubmed_count': int(row.get('pubmed_count', 0) or 0),
                        'direct_evidence_codes': codes
                    })
                if diseases:
                    self.validated_connections['gene_to_disease'][gene] = diseases
            
            gene_connections = sum(len(v) for v in self.validated_connections['gene_to_disease'].values())
            print(f"Connected {len(self.validated_connections['gene_to_disease'])} genes to diseases ({gene_connections} relationships)")
        else:
            print("Skipping Gene → Disease connections (dataset not available)")
        
        # Step 4: 완전한 경로 구축 (관대한 기준)
        self._build_complete_pathways()
    
    def _build_complete_pathways(self):
        """완전한 경로 구축"""
        print("Building complete pathways...")
        
        pathway_count = 0
        pathway_qualities = []
        
        foods_sorted = sorted(
            self.validated_connections['food_to_nutrient'].items(),
            key=lambda item: (-len(item[1]), item[0])
        )
        max_total = self.max_pathways if self.max_pathways else float('inf')
        stop_building = False
        
        for food, nutrient_records in foods_sorted:
            if stop_building:
                break
            food_pathways = []
            
            # 식품의 영양소들
            for nut_info in nutrient_records[:self.max_nutrients_per_food]:
                nutrient = nut_info['nutrient']
                
                if nutrient in self.validated_connections['nutrient_to_gene']:
                    # 영양소의 유전자들
                    for gene_info in self.validated_connections['nutrient_to_gene'][nutrient][:self.max_genes_per_nutrient]:
                        gene = gene_info['gene']
                        
                        if gene in self.validated_connections['gene_to_disease']:
                            direction_alignment = self._compute_direction_alignment(gene, gene_info.get('direction'))
                            if direction_alignment < self.min_direction_alignment:
                                continue
                            enriched_gene_info = dict(gene_info)
                            enriched_gene_info['direction_alignment'] = direction_alignment
                            # 유전자의 질병들
                            for disease_info in self.validated_connections['gene_to_disease'][gene][:self.max_diseases_per_gene]:
                                
                                quality_score, components = self._calculate_pathway_quality_score(nut_info, enriched_gene_info, disease_info)
                                
                                if quality_score >= self.quality_thresholds['min_pathway_quality']:
                                    pathway = {
                                        'food': food,
                                        'nutrient': nutrient,
                                        'nutrient_amount': nut_info['amount'],
                                        'gene': gene,
                                        'chemical_id': gene_info['chemical_id'],
                                        'ref_count': gene_info['ref_count'],
                                        'confidence': gene_info['confidence'],
                                        'disease': disease_info['disease'],
                                        'disease_id': disease_info['disease_id'],
                                        'has_direct_evidence': disease_info['has_direct_evidence'],
                                        'inference_score': disease_info['inference_score'],
                                        'evidence_count': disease_info.get('evidence_count', 1),
                                        'pubmed_count': disease_info.get('pubmed_count', 0),
                                        'direct_evidence_codes': "|".join(disease_info.get('direct_evidence_codes', [])),
                                        'pathway_quality_score': quality_score,
                                        'evidence_level': self._classify_evidence_level(
                                            quality_score,
                                            disease_info.get('has_direct_evidence'),
                                            disease_info.get('evidence_count', 1),
                                            components['inference'],
                                            disease_info.get('direct_evidence_codes', []),
                                            disease_info.get('pubmed_count', 0)
                                        ),
                                        'score_dose': components['dose'],
                                        'score_mapping': components['mapping'],
                                        'score_literature': components['literature'],
                                        'score_inference': components['inference'],
                                        'score_support': components['support'],
                                        'score_direction': components['direction'],
                                        'score_pubmed': components['pubmed'],
                                        'direction_alignment': direction_alignment,
                                        'gene_direction': enriched_gene_info.get('direction', 'unknown'),
                                        'gwas_factor': components.get('gwas_factor', 1.0)
                                    }
                                    
                                    food_pathways.append(pathway)
                                    pathway_qualities.append(quality_score)
                                    pathway_count += 1
                                    
                                    if pathway_count >= max_total:
                                        stop_building = True
                                        break
                            if stop_building:
                                break
                    if stop_building:
                        break
            
            if food_pathways:
                self.validated_connections['complete_pathways'][food] = food_pathways
            if stop_building:
                break
        
        # 통계 출력
        if pathway_qualities:
            avg_quality = np.mean(pathway_qualities)
            high_quality = sum(1 for q in pathway_qualities if q >= 0.75)
            moderate_quality = sum(
                1 for q in pathway_qualities
                if self.quality_thresholds['min_pathway_quality'] <= q < 0.75
            )
            
            print(f"Built {pathway_count} complete pathways")
            print(f"Average quality: {avg_quality:.3f}")
            print(f"High quality (≥0.75): {high_quality}")
            print(f"Moderate quality (≥{self.quality_thresholds['min_pathway_quality']:.2f} & <0.75): {moderate_quality}")
        else:
            print("No pathways built - check quality thresholds")
    
    def _calculate_pathway_quality_score(self, nut_info, gene_info, disease_info) -> Tuple[float, Dict[str, float]]:
        """데이터 기반 정규화를 적용한 품질 점수 계산"""
        amount_component = self._sigmoid_scale(nut_info.get('amount', 0.0), 'nutrient_amount', use_log=True)
        ratio_component = self._dose_ratio_component(nut_info.get('nutrient'), nut_info.get('amount'))
        dv_ratio_component = np.clip(nut_info.get('dv_ratio', 0.0), 0.0, 1.0)
        dose_component = 0.4 * amount_component + 0.4 * ratio_component + 0.2 * dv_ratio_component

        mapping_confidence = float(np.clip(gene_info.get('confidence', 0.0), 0.0, 1.0))
        literature_component = self._sigmoid_scale(gene_info.get('ref_count', 0.0), 'reference_count', use_log=True)
        inference_component = self._sigmoid_scale(disease_info.get('inference_score', 0.0), 'inference_score')
        gwas_factor = float(max(1.0, gene_info.get('gwas_support', 1.0)))
        inference_component = float(np.clip(inference_component * gwas_factor, 0.0, 1.0))

        evidence_count = float(disease_info.get('evidence_count', 1) or 1)
        pubmed_count = float(disease_info.get('pubmed_count', 0) or 0)
        pubmed_component = np.clip(pubmed_count / 5.0, 0.0, 1.0)
        replication_component = np.clip((evidence_count - 1.0) / 3.0 + 0.5 * pubmed_component, 0.0, 1.0)

        direct_codes = disease_info.get('direct_evidence_codes', []) or []
        if disease_info.get('has_direct_evidence') and direct_codes:
            prior = max(DIRECT_EVIDENCE_PRIORS.get(code, 0.7) for code in direct_codes)
            support_component = min(1.0, prior + 0.3 * replication_component)
        elif disease_info.get('has_direct_evidence'):
            support_component = min(1.0, 0.65 + 0.35 * replication_component)
        else:
            support_component = 0.1 + 0.45 * replication_component
            support_component = min(0.6, support_component)

        direction_alignment = float(gene_info.get('direction_alignment', 0.5))
        direction_component = np.clip(direction_alignment, 0.0, 1.0)

        components = {
            'dose': dose_component,
            'mapping': mapping_confidence,
            'literature': literature_component,
            'inference': inference_component,
            'support': support_component,
            'direction': direction_component,
            'pubmed': pubmed_component
        }
        components['gwas_factor'] = np.clip(gwas_factor, 1.0, 1.5)

        score = sum(self.score_weights[k] * components.get(k, 0.0) for k in self.score_weights)

        return float(np.clip(score, 0.0, 1.0)), components
    
    def _classify_evidence_level(self,
                                 score: float,
                                 has_direct: bool,
                                 evidence_count: int,
                                 inference_component: float,
                                 direct_codes: List[str],
                                 pubmed_count: int) -> str:
        """증거 수준 분류"""
        evidence_count = int(max(0, evidence_count))
        pubmed_count = int(max(0, pubmed_count))
        direct_codes = [c.lower() for c in (direct_codes or [])]

        if has_direct:
            has_therapeutic = any("therapeutic" in c for c in direct_codes)
            has_marker = any("marker" in c or "mechanism" in c for c in direct_codes)
            if has_therapeutic and (pubmed_count >= 2 or evidence_count >= 3) and score >= 0.65:
                return "Strong Evidence"
            if (has_therapeutic or has_marker) and score >= 0.55:
                return "Moderate Evidence"
            return "Preliminary Evidence"

        if evidence_count >= 3 and pubmed_count >= 3 and score >= 0.6 and inference_component >= 0.65:
            return "Moderate Evidence"
        if evidence_count >= 2 and score >= 0.4:
            return "Preliminary Evidence"
        if score >= 0.3:
            return "Preliminary Evidence"
        return "Weak Evidence"
    
    def export_results(self, output_dir: str = "nutrigenomics_csv_exports"):
        """분석 결과 CSV 내보내기"""
        print("Exporting results to CSV...")
        
        export_path = Path(output_dir)
        export_path.mkdir(exist_ok=True)
        
        # 완전한 경로들
        all_pathways = []
        for food, pathways in self.validated_connections['complete_pathways'].items():
            all_pathways.extend(pathways)
        
        if all_pathways:
            pathways_df = pd.DataFrame(all_pathways)
            pathways_df.to_csv(export_path / "pathways.csv", index=False)
            print(f"Exported {len(all_pathways)} pathways")
            
            # 품질별 분류
            high_quality = [p for p in all_pathways if p['pathway_quality_score'] >= 0.75]
            moderate_quality = [
                p for p in all_pathways
                if self.quality_thresholds['min_pathway_quality'] <= p['pathway_quality_score'] < 0.75
            ]
            
            if high_quality:
                pd.DataFrame(high_quality).to_csv(export_path / "high_quality_pathways.csv", index=False)
                print(f"Exported {len(high_quality)} high-quality pathways")
            
            if moderate_quality:
                pd.DataFrame(moderate_quality).to_csv(export_path / "moderate_quality_pathways.csv", index=False)
        
        # 연결 통계
        quality_scores = [p['pathway_quality_score'] for p in all_pathways] if all_pathways else []
        evidence_counts = [p.get('evidence_count', 1) for p in all_pathways] if all_pathways else []
        pubmed_counts = [p.get('pubmed_count', 0) for p in all_pathways] if all_pathways else []
        direct_evidence_total = sum(1 for p in all_pathways if p.get('has_direct_evidence'))

        unique_foods = {p['food'] for p in all_pathways}
        unique_nutrients = {p['nutrient'] for p in all_pathways}
        unique_genes = {p['gene'] for p in all_pathways}

        stats = {
            'total_pathways': len(all_pathways),
            'unique_foods': len(unique_foods),
            'unique_nutrients': len(unique_nutrients),
            'unique_genes': len(unique_genes),
            'avg_quality_score': float(np.mean(quality_scores)) if quality_scores else 0.0,
            'median_quality_score': float(np.median(quality_scores)) if quality_scores else 0.0,
            'p90_quality_score': float(np.percentile(quality_scores, 90)) if quality_scores else 0.0,
            'high_quality_count': len([p for p in all_pathways if p['pathway_quality_score'] >= 0.75]),
            'moderate_quality_count': len([
                p for p in all_pathways
                if self.quality_thresholds['min_pathway_quality'] <= p['pathway_quality_score'] < 0.75
            ]),
            'direct_evidence_count': direct_evidence_total,
            'mean_evidence_support': float(np.mean(evidence_counts)) if evidence_counts else 0.0,
            'max_evidence_support': int(np.max(evidence_counts)) if evidence_counts else 0,
            'mean_pubmed_support': float(np.mean(pubmed_counts)) if pubmed_counts else 0.0,
            'max_pubmed_support': int(np.max(pubmed_counts)) if pubmed_counts else 0
        }

        print("[SUMMARY] Pathway inventory")
        print(f" - total pathways: {stats['total_pathways']} (high>=0.75: {stats['high_quality_count']}, moderate≥{self.quality_thresholds['min_pathway_quality']:.2f}: {stats['moderate_quality_count']})")
        print(f" - unique foods: {stats['unique_foods']} | nutrients: {stats['unique_nutrients']} | genes: {stats['unique_genes']}")
        print(f" - mean quality: {stats['avg_quality_score']:.3f} (median {stats['median_quality_score']:.3f}, p90 {stats['p90_quality_score']:.3f})")
        print(f" - mean CTD associations per pathway: {stats['mean_evidence_support']:.2f} (max {stats['max_evidence_support']})")
        print(f" - mean PubMed refs per pathway: {stats['mean_pubmed_support']:.2f} (max {stats['max_pubmed_support']})")
        print(f" - pathways with direct evidence: {stats['direct_evidence_count']}")

        stats_df = pd.DataFrame([stats])
        stats_df.to_csv(export_path / "summary_stats.csv", index=False)
        print("Exported summary statistics")
        
        return export_path

    def create_network_visualization(self):
        """네트워크 시각화"""
        print("Creating network visualization...")
        
        all_pathways = []
        for pathways in self.validated_connections['complete_pathways'].values():
            all_pathways.extend(pathways)
        
        if not all_pathways:
            return go.Figure().add_annotation(
                text="No pathways generated - check data and thresholds",
                xref="paper", yref="paper", x=0.5, y=0.5,
                font=dict(size=16, color="red")
            )
        
        # 상위 품질 경로들로 네트워크 구성
        top_pathways = sorted(all_pathways, key=lambda x: x['pathway_quality_score'], reverse=True)[:50]
        
        G = nx.Graph()
        layer_map = {'food': 0, 'nutrient': 1, 'gene': 2, 'disease': 3}
        rng = np.random.default_rng(42)
        
        # 노드 및 엣지 추가
        for pathway in top_pathways:
            food = pathway['food']
            nutrient = pathway['nutrient']
            gene = pathway['gene']
            disease = pathway['disease']
            quality = pathway['pathway_quality_score']
            
            # 노드 추가 (품질에 따른 크기)
            base_size = 8 + (quality * 15)
            G.add_node(food, type='food', size=base_size + 8, quality=quality)
            G.add_node(nutrient, type='nutrient', size=base_size + 4, quality=quality)
            G.add_node(gene, type='gene', size=base_size, quality=quality)
            G.add_node(disease, type='disease', size=max(base_size - 2, 6), quality=quality)
            
            # 엣지 추가
            G.add_edge(food, nutrient, weight=quality)
            G.add_edge(nutrient, gene, weight=quality)
            G.add_edge(gene, disease, weight=quality)
        
        if len(G) == 0:
            return go.Figure().add_annotation(
                text="No pathways generated - check data and thresholds",
                xref="paper", yref="paper", x=0.5, y=0.5,
                font=dict(size=16, color="red")
            )

        nx.set_node_attributes(G, {node: layer_map.get(attrs.get('type'), 0) for node, attrs in G.nodes(data=True)}, 'subset')
        pos = nx.multipartite_layout(G, subset_key='subset', align='horizontal', scale=2.0)
        jitter = dict(zip(G.nodes(), rng.normal(scale=0.02, size=(len(G), 2))))
        pos = {node: (coords[0] + jitter[node][0], coords[1] + jitter[node][1]) for node, coords in pos.items()}
        
        edge_trace = []
        for node_u, node_v, data in G.edges(data=True):
            x0, y0 = pos[node_u]
            x1, y1 = pos[node_v]
            weight = max(data.get('weight', 0.1), 0.05)
            edge_trace.append(go.Scatter(
                x=[x0, x1, None], y=[y0, y1, None],
                mode='lines',
                line=dict(
                    width=1.2 + weight * 1.6,
                    color=f'rgba(52, 73, 94,{0.35 + min(weight, 0.9) * 0.45})'
                ),
                hoverinfo='none',
                showlegend=False
            ))
        
        colors = {'food': '#E74C3C', 'nutrient': '#1ABC9C', 'gene': '#2980B9', 'disease': '#F39C12'}
        node_traces = []
        for node_type in ['food', 'nutrient', 'gene', 'disease']:
            nodes_of_type = [node for node, attrs in G.nodes(data=True) if attrs.get('type') == node_type]
            if not nodes_of_type:
                continue
            
            x_coords = [pos[node][0] for node in nodes_of_type]
            y_coords = [pos[node][1] for node in nodes_of_type]
            sizes = [G.nodes[node].get('size', 10) for node in nodes_of_type]
            qualities = [G.nodes[node].get('quality', 0.0) for node in nodes_of_type]
            texts = [node if len(node) <= 16 else f"{node[:13]}..." for node in nodes_of_type]
            
            node_traces.append(go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='markers+text',
                marker=dict(
                    size=sizes,
                    color=colors[node_type],
                    opacity=0.9,
                    line=dict(width=1.2, color='#FFFFFF'),
                    sizemode='diameter'
                ),
                text=texts,
                textposition="middle center",
                textfont=dict(size=8, color='#FFFFFF'),
                name=node_type.title(),
                hovertemplate=(
                    f"<b>{node_type.title()}</b><br>"
                    + "Name: %{customdata[0]}<br>"
                    + "Quality Score: %{customdata[1]:.3f}<br>"
                    + "<extra></extra>"
                ),
                customdata=list(zip(nodes_of_type, qualities))
            ))
        
        fig = go.Figure(data=edge_trace + node_traces)
        fig.update_layout(
            title=f"<b>KFRI NutriGenomics Network</b> &nbsp;|&nbsp; Top {len(top_pathways)} Pathways",
            showlegend=True,
            hovermode='closest',
            height=700,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='#FBFCFC',
            paper_bgcolor='#FBFCFC',
            legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='center', x=0.5),
            margin=dict(l=30, r=30, t=80, b=30)
        )
        
        return fig

    @staticmethod
    def _render_table(headers: List[str], rows: List[List[str]]) -> str:
        header_html = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
        if rows:
            body_html = "".join(
                "<tr>" + "".join(f"<td>{html.escape(str(cell))}</td>" for cell in row) + "</tr>"
                for row in rows
            )
        else:
            body_html = f"<tr><td colspan='{len(headers)}' class='empty-cell'>No data available</td></tr>"
        return (
            "<table class='summary-table'>"
            "<thead><tr>" + header_html + "</tr></thead>"
            "<tbody>" + body_html + "</tbody>"
            "</table>"
        )
    
    def generate_report(self, output_file: str = "food_health_report.html"):
        """종합 리포트 생성"""
        print("Generating comprehensive report...")
        
        start_time = time.time()
        
        # 데이터 처리
        self.load_data_efficiently()
        self.build_connections_efficiently()
        
        # CSV 내보내기
        csv_path = self.export_results()
        
        # 시각화
        network_fig = self.create_network_visualization()
        
        # HTML 변환
        network_html = pyo.plot(network_fig, output_type='div', include_plotlyjs=False)
        
        # 통계 계산
        all_pathways = []
        for pathways in self.validated_connections['complete_pathways'].values():
            all_pathways.extend(pathways)
        
        total_pathways = len(all_pathways)
        total_foods = len(self.validated_connections['complete_pathways'])
        total_nutrients = len(self.validated_connections['nutrient_to_gene'])
        total_genes = len(self.validated_connections['gene_to_disease'])
        
        pathways_df = pd.DataFrame(all_pathways) if all_pathways else pd.DataFrame()
        if not pathways_df.empty:
            avg_quality = float(pathways_df['pathway_quality_score'].mean())
            median_quality = float(pathways_df['pathway_quality_score'].median())
            p90_quality = float(pathways_df['pathway_quality_score'].quantile(0.9))
            high_quality = int((pathways_df['pathway_quality_score'] >= 0.75).sum())
            with_direct_evidence = int(pathways_df['has_direct_evidence'].sum())
            evidence_level_counts = pathways_df['evidence_level'].value_counts().to_dict()
            avg_evidence_support = float(pathways_df['evidence_count'].mean())
            max_evidence_support = int(pathways_df['evidence_count'].max())
            avg_pubmed_support = float(pathways_df['pubmed_count'].mean())
            max_pubmed_support = int(pathways_df['pubmed_count'].max())
        else:
            avg_quality = 0.0
            median_quality = 0.0
            p90_quality = 0.0
            high_quality = 0
            with_direct_evidence = 0
            evidence_level_counts = {}
            avg_evidence_support = 0.0
            max_evidence_support = 0
            avg_pubmed_support = 0.0
            max_pubmed_support = 0
        
        processing_time = time.time() - start_time
        csv_dir_display = html.escape(str(csv_path))
        
        quality_points = [
            ("Average quality", f"{avg_quality:.3f}"),
            ("Median quality", f"{median_quality:.3f}"),
            ("90th percentile", f"{p90_quality:.3f}"),
            ("High-quality pathways", f"{high_quality:,}"),
            ("Direct evidence pathways", f"{with_direct_evidence:,}")
        ]
        quality_list_html = "".join(
            f"<li><span>{html.escape(label)}</span><strong>{html.escape(value)}</strong></li>"
            for label, value in quality_points
        )
        evidence_items = [
            f"<li><span>{html.escape(level)}</span><strong>{count:,}</strong></li>"
            for level, count in sorted(evidence_level_counts.items(), key=lambda x: (-x[1], x[0]))
        ]
        if avg_evidence_support > 0:
            evidence_items.append(
                f"<li><span>Avg CTD associations</span><strong>{avg_evidence_support:.1f}</strong></li>"
            )
            evidence_items.append(
                f"<li><span>Max CTD associations</span><strong>{max_evidence_support:,}</strong></li>"
            )
        if avg_pubmed_support > 0:
            evidence_items.append(
                f"<li><span>Avg PubMed references</span><strong>{avg_pubmed_support:.1f}</strong></li>"
            )
            evidence_items.append(
                f"<li><span>Max PubMed references</span><strong>{max_pubmed_support:,}</strong></li>"
            )
        evidence_list_html = "".join(evidence_items) or "<li><span>Evidence levels</span><strong>No pathways</strong></li>"
        
        if not pathways_df.empty:
            food_summary_df = (
                pathways_df.groupby('food')
                .agg(pathway_count=('food', 'size'), avg_quality=('pathway_quality_score', 'mean'))
                .sort_values(['pathway_count', 'avg_quality'], ascending=[False, False])
                .head(5)
                .reset_index()
            )
            top_food_rows = [
                [str(idx + 1), row['food'], f"{int(row['pathway_count']):,}", f"{row['avg_quality']:.3f}"]
                for idx, row in food_summary_df.iterrows()
            ]
            disease_summary_df = (
                pathways_df.groupby('disease')
                .agg(pathway_count=('disease', 'size'), avg_quality=('pathway_quality_score', 'mean'))
                .sort_values(['pathway_count', 'avg_quality'], ascending=[False, False])
                .head(5)
                .reset_index()
            )
            top_disease_rows = [
                [str(idx + 1), row['disease'], f"{int(row['pathway_count']):,}", f"{row['avg_quality']:.3f}"]
                for idx, row in disease_summary_df.iterrows()
            ]
            top_pathways_rows = [
                [
                    rec.get('food', ''),
                    rec.get('nutrient', ''),
                    rec.get('gene', ''),
                    rec.get('disease', ''),
                    f"{rec.get('pathway_quality_score', 0):.3f}",
                    rec.get('evidence_level', ''),
                    int(rec.get('evidence_count', 0))
                ]
                for rec in pathways_df.sort_values('pathway_quality_score', ascending=False).head(8).to_dict('records')
            ]
        else:
            top_food_rows = []
            top_disease_rows = []
            top_pathways_rows = []
        
        top_foods_table = self._render_table(["#", "Food", "Pathways", "Avg Quality"], top_food_rows)
        top_diseases_table = self._render_table(["#", "Disease", "Pathways", "Avg Quality"], top_disease_rows)
        top_pathways_table = self._render_table(
            ["Food", "Nutrient", "Gene", "Disease", "Quality", "Evidence Level", "Evidence Count"],
            top_pathways_rows
        )
        
        if with_direct_evidence:
            direct_evidence_text = f"{with_direct_evidence:,} pathways with direct evidence"
        else:
            direct_evidence_text = "Direct evidence not detected"
        if avg_evidence_support > 0:
            direct_evidence_text += f" · mean CTD associations {avg_evidence_support:.1f}"
        if avg_pubmed_support > 0:
            direct_evidence_text += f" · mean PubMed refs {avg_pubmed_support:.1f}"
        processing_time_text = f"{processing_time:.1f} seconds"
        
        html_template = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <meta charset="utf-8">
            <title>KFRI NutriGenomics Food-Health Database Analysis</title>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <style>
                :root {{
                    --bg-gradient: linear-gradient(135deg, #0F2027 0%, #203A43 45%, #2C5364 100%);
                    --card-bg: rgba(255, 255, 255, 0.12);
                    --accent: #1ABC9C;
                    --accent-2: #F39C12;
                    --text-light: #ECF0F1;
                    --text-muted: rgba(236, 240, 241, 0.76);
                    --table-border: rgba(255, 255, 255, 0.15);
                }}
                * {{ box-sizing: border-box; }}
                body {{
                    margin: 0;
                    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                    background: var(--bg-gradient);
                    color: var(--text-light);
                }}
                .page {{
                    max-width: 1300px;
                    margin: 0 auto;
                    padding: 48px 36px 64px;
                }}
                .hero {{
                    text-align: center;
                    padding: 48px 32px;
                    border-radius: 28px;
                    background: rgba(255, 255, 255, 0.06);
                    box-shadow: 0 30px 60px rgba(0, 0, 0, 0.35);
                    backdrop-filter: blur(16px);
                }}
                .hero h1 {{
                    margin: 0;
                    font-size: 2.8rem;
                    letter-spacing: 0.03em;
                }}
                .hero p {{
                    margin: 12px 0 0;
                    font-size: 1.1rem;
                    color: var(--text-muted);
                }}
                .hero .badge {{
                    display: inline-block;
                    margin-top: 24px;
                    padding: 10px 18px;
                    border-radius: 999px;
                    background: rgba(26, 188, 156, 0.18);
                    color: #E8F8F5;
                    font-size: 0.95rem;
                    letter-spacing: 0.04em;
                }}
                .metric-grid {{
                    margin-top: 48px;
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                    gap: 24px;
                }}
                .metric-card {{
                    background: var(--card-bg);
                    padding: 24px;
                    border-radius: 20px;
                    border: 1px solid rgba(255, 255, 255, 0.08);
                    box-shadow: 0 18px 40px rgba(0, 0, 0, 0.25);
                }}
                .metric-label {{
                    font-size: 0.9rem;
                    text-transform: uppercase;
                    letter-spacing: 0.08em;
                    color: var(--text-muted);
                    margin-bottom: 12px;
                }}
                .metric-value {{
                    font-size: 2.3rem;
                    font-weight: 600;
                    color: var(--text-light);
                }}
                .content-section {{
                    margin-top: 56px;
                }}
                .section-title {{
                    font-size: 2.1rem;
                    margin-bottom: 24px;
                    color: #FDFEFE;
                    letter-spacing: 0.02em;
                }}
                .insight-grid {{
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(260px, 1fr));
                    gap: 24px;
                    margin-bottom: 36px;
                }}
                .insight-card {{
                    background: var(--card-bg);
                    padding: 24px;
                    border-radius: 20px;
                    border: 1px solid rgba(255, 255, 255, 0.08);
                    box-shadow: 0 18px 40px rgba(0, 0, 0, 0.25);
                }}
                .insight-card h3 {{
                    margin: 0 0 16px;
                    font-size: 1.2rem;
                }}
                .insight-card ul {{
                    padding: 0;
                    margin: 0;
                    list-style: none;
                }}
                .insight-card li {{
                    display: flex;
                    justify-content: space-between;
                    align-items: baseline;
                    padding: 10px 0;
                    border-bottom: 1px solid rgba(255, 255, 255, 0.08);
                    font-size: 0.98rem;
                }}
                .insight-card li:last-child {{
                    border-bottom: none;
                }}
                .insight-card span {{
                    color: var(--text-muted);
                }}
                .chart-section {{
                    background: rgba(255, 255, 255, 0.08);
                    border-radius: 24px;
                    padding: 32px;
                    box-shadow: 0 24px 50px rgba(0, 0, 0, 0.35);
                    margin-bottom: 36px;
                }}
                .summary-grid {{
                    display: grid;
                    grid-template-columns: repeat(auto-fit, minmax(340px, 1fr));
                    gap: 24px;
                }}
                .summary-table {{
                    width: 100%;
                    border-collapse: collapse;
                    background: var(--card-bg);
                    border-radius: 18px;
                    overflow: hidden;
                    border: 1px solid var(--table-border);
                }}
                .summary-table th, .summary-table td {{
                    padding: 14px 16px;
                    border-bottom: 1px solid var(--table-border);
                    text-align: left;
                    font-size: 0.95rem;
                }}
                .summary-table th {{
                    background: rgba(255, 255, 255, 0.06);
                    text-transform: uppercase;
                    letter-spacing: 0.05em;
                    font-size: 0.85rem;
                }}
                .summary-table tbody tr:nth-child(odd) {{
                    background: rgba(255, 255, 255, 0.03);
                }}
                .summary-table td:last-child {{
                    font-weight: 600;
                }}
                .empty-cell {{
                    text-align: center;
                    color: var(--text-muted);
                    padding: 20px;
                }}
                .callout {{
                    margin-top: 36px;
                    padding: 18px 24px;
                    border-radius: 18px;
                    background: rgba(243, 156, 18, 0.15);
                    color: #FAD7A0;
                    border: 1px solid rgba(243, 156, 18, 0.4);
                    display: flex;
                    align-items: center;
                    gap: 16px;
                }}
                .callout strong {{
                    color: #FDEBD0;
                    font-weight: 600;
                }}
                .footer {{
                    margin-top: 64px;
                    text-align: center;
                    color: var(--text-muted);
                    font-size: 0.9rem;
                }}
                @media (max-width: 720px) {{
                    .page {{
                        padding: 32px 18px;
                    }}
                    .hero {{
                        padding: 36px 24px;
                    }}
                    .metric-grid {{
                        grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
                    }}
                }}
            </style>
        </head>
        <body>
            <div class="page">
                <header class="hero">
                    <p>KFRI NutriGenomics · Knowledge Graph</p>
                    <h1>Food → Nutrient → Gene → Disease Pathway Summary</h1>
                    <div class="badge">Validation thresholds tuned to empirical distributions</div>
                    
                    <div class="metric-grid">
                        <div class="metric-card">
                            <div class="metric-label">Total pathways</div>
                            <div class="metric-value">{total_pathways:,}</div>
                        </div>
                        <div class="metric-card">
                            <div class="metric-label">Foods connected</div>
                            <div class="metric-value">{total_foods:,}</div>
                        </div>
                        <div class="metric-card">
                            <div class="metric-label">Nutrients mapped</div>
                            <div class="metric-value">{total_nutrients:,}</div>
                        </div>
                        <div class="metric-card">
                            <div class="metric-label">Genes connected</div>
                            <div class="metric-value">{total_genes:,}</div>
                        </div>
                        <div class="metric-card">
                            <div class="metric-label">Avg quality score</div>
                            <div class="metric-value">{avg_quality:.3f}</div>
                        </div>
                        <div class="metric-card">
                            <div class="metric-label">Compute time</div>
                            <div class="metric-value">{processing_time_text}</div>
                        </div>
                    </div>
                </header>
                
                <section class="content-section">
                    <h2 class="section-title">Insight Highlights</h2>
                    <div class="insight-grid">
                        <div class="insight-card">
                            <h3>Quality Distribution</h3>
                            <ul>{quality_list_html}</ul>
                        </div>
                        <div class="insight-card">
                            <h3>Evidence Levels</h3>
                            <ul>{evidence_list_html}</ul>
                        </div>
                        <div class="insight-card">
                            <h3>Direct Evidence</h3>
                            <p>{direct_evidence_text}</p>
                            <p>CSV exports available in <strong>{csv_dir_display}</strong></p>
                        </div>
                    </div>
                </section>
                
                <section class="content-section">
                    <h2 class="section-title">KFRI NutriGenomics Network Visualization</h2>
                    <div class="chart-section">
                        {network_html}
                    </div>
                </section>
                
                <section class="content-section">
                    <h2 class="section-title">Top Contributors</h2>
                    <div class="summary-grid">
                        <div>{top_foods_table}</div>
                        <div>{top_diseases_table}</div>
                    </div>
                    <div class="summary-grid" style="margin-top:24px;">
                        <div style="grid-column:1 / -1;">
                            {top_pathways_table}
                        </div>
                    </div>
                </section>
                
                <div class="callout">
                    <strong>Evidence integration:</strong>
                    <span>Scoring now blends nutrient exposure, mapping confidence, literature strength, and inference evidence through data-driven normalization for each pathway.</span>
                </div>
                
                <footer class="footer">
                    Generated in {processing_time_text} · Plotly figure elements: {len(network_fig.data)} · High quality pathways: {high_quality:,}
                </footer>
            </div>
        </body>
        </html>
        """
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_template)
        
        print(f"\nReport generated: {output_file}")
        print(f"Results: {total_pathways} pathways, {high_quality} high quality")
        print(f"CSV data exported to: {csv_path}")
        
        return output_file

def main():
    print("KFRI NutriGenomics Food-Health Analyzer")
    print("Applying data-driven validation standards...")
    
    analyzer = KFRINutriGenomicsAnalyzer()
    output_file = analyzer.generate_report()
    
    print(f"\nOpen '{output_file}' to view the analysis!")

if __name__ == "__main__":
    main()
