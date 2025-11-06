#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Compound->Gene->Disease Pathway Analyzer (v5.8, 2025-09-15)

Key updates (v5.8)
- Keeps the PubMed response disk cache (v5.7) and enriches entity/context tracking for species, study design, and dose.
- Splits the Gene/Disease layers into Drug->Gene / Gene->Disease / Drug->Disease tables with rule annotations.
- Surfaces evidence snippets and species tags in reports/CSVs while branching non-human models by context specificity.
- Reduces coffee/caffeine ambiguity and tags non-HGNC gene symbols with their species.
- Preserves prior v5.6~v5.7 safeguards (aggressive organ/specimen drops, stable up/down labels, etc.).
"""

from __future__ import annotations
import argparse, csv, datetime as dt, json, os, re, time, statistics as stats, sys, hashlib
from xml.sax.saxutils import escape
from dataclasses import dataclass, asdict, field
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple, Any
import xml.etree.ElementTree as ET
import requests
from tqdm.auto import tqdm

# ---------- Optional plotting/report ----------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import networkx as nx
    HAVE_GRAPH = True
except Exception:
    HAVE_GRAPH = False

try:
    from reportlab.lib.pagesizes import A4, LETTER
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    HAVE_REPORTLAB = True
except Exception:
    HAVE_REPORTLAB = False

# ---------------- Constants ----------------
NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
HGNC_TSV_GCS  = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_TSV_FTP  = "https://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt"

DEFAULT_DELAY_NO_KEY = 0.35
DEFAULT_DELAY_WITH_KEY = 0.12

SENT_SPLIT = re.compile(r"(?<=[\.\?\!])\s+(?=[A-Z(])", re.UNICODE)
DEFAULT_CACHE_DIR = os.path.join("cache", "pubmed")

SPECIES_MAP = {
    "humans": ("Homo sapiens", "Human"),
    "human": ("Homo sapiens", "Human"),
    "male": ("Homo sapiens", "Human"),
    "female": ("Homo sapiens", "Human"),
    "mice": ("Mus musculus", "Mouse"),
    "rats": ("Rattus norvegicus", "Rat"),
    "dogs": ("Canis lupus familiaris", "Dog"),
    "cats": ("Felis catus", "Cat"),
    "rabbits": ("Oryctolagus cuniculus", "Rabbit"),
    "caenorhabditis elegans": ("Caenorhabditis elegans", "C. elegans"),
    "c. elegans": ("Caenorhabditis elegans", "C. elegans"),
    "zebrafish": ("Danio rerio", "Zebrafish"),
    "cell line": ("Cell culture", "Cellular"),
    "cell lines": ("Cell culture", "Cellular"),
    "organ culture techniques": ("Organ culture", "Organ culture")
}
SPECIES_PRIORITY = ["Human","Mouse","Rat","C. elegans","Zebrafish","Dog","Cat","Rabbit","Cellular","Organ culture"]

DOSE_PATTERN = re.compile(r"\b\d+(?:\.\d+)?\s*(?:mg|μg|ug|g|kg|µg|ng|ml|mM|µM|uM|%)(?:\s*/\s*(?:kg|day|d|h))?", re.IGNORECASE)
TIME_PATTERN = re.compile(r"\b\d+(?:\.\d+)?\s*(?:h|hr|hrs|hour|hours|day|days|week|weeks|month|months)\b", re.IGNORECASE)

# ---------- caching ----------
class PubmedCache:
    """Lightweight JSON file cache for PubMed E-utilities responses."""
    def __init__(self, base_dir: str = DEFAULT_CACHE_DIR, ttl_days: int = 7,
                 enabled: bool = True, force_refresh: bool = False):
        self.enabled = enabled
        self.force_refresh = force_refresh
        self.base_dir = base_dir
        self.ttl = None if ttl_days <= 0 else dt.timedelta(days=ttl_days)
        self.metrics = Counter()
        if self.enabled:
            os.makedirs(self.base_dir, exist_ok=True)

    def _path(self, kind: str, key: str) -> str:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", key)
        return os.path.join(self.base_dir, f"{kind}_{safe}.json")

    def _is_fresh(self, path: str) -> bool:
        if not os.path.exists(path):
            return False
        if self.force_refresh:
            return False
        if not self.ttl:
            return True
        mtime = dt.datetime.fromtimestamp(os.path.getmtime(path))
        return (dt.datetime.now() - mtime) <= self.ttl

    def _read_json(self, path: str) -> Optional[Dict[str, Any]]:
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return None

    def _write_json(self, path: str, payload: Dict[str, Any]) -> None:
        tmp = f"{path}.tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False)
        os.replace(tmp, path)

    def snapshot(self) -> Dict[str, int]:
        return dict(self.metrics)

    def fetch_esearch(self, key: str) -> Optional[List[str]]:
        if not self.enabled:
            self.metrics["esearch_miss"] += 1
            return None
        path = self._path("esearch", key)
        if self._is_fresh(path):
            data = self._read_json(path)
            if data and isinstance(data.get("pmids"), list):
                self.metrics["esearch_hit"] += 1
                return data["pmids"]
        self.metrics["esearch_miss"] += 1
        return None

    def store_esearch(self, key: str, pmids: List[str]) -> None:
        if not self.enabled:
            return
        path = self._path("esearch", key)
        self._write_json(path, {"pmids": pmids, "saved": dt.datetime.now().isoformat()})

    def fetch_detail(self, pmid: str) -> Optional[Dict[str, Any]]:
        if not self.enabled:
            self.metrics["efetch_miss"] += 1
            return None
        path = self._path("pmid", pmid)
        if self._is_fresh(path):
            data = self._read_json(path)
            if data:
                self.metrics["efetch_hit"] += 1
                return data
        self.metrics["efetch_miss"] += 1
        return None

    def store_detail(self, pmid: str, payload: Dict[str, Any]) -> None:
        if not self.enabled:
            return
        path = self._path("pmid", pmid)
        self._write_json(path, payload)

# ---------- interaction cues ----------
UP_PATS   = [re.compile(p, re.IGNORECASE) for p in [
    r"\bup[\-\s]?regulat\w+", r"\bincreas\w+", r"\benhanc\w+", r"\belevat\w+",
    r"\bactivat\w+", r"\bover[-\s]?express\w+", r"\binduc\w+", r"\baugmen\w+"
]]
DOWN_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bdown[\-\s]?regulat\w+", r"\bdecreas\w+", r"\breduc\w+", r"\bsuppress\w+",
    r"\binhibit\w+", r"\bknock[-\s]?down\w*", r"\bsilenc\w+", r"\bloss\b"
]]
BIND_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bbind\w+\b", r"\binteract\w+\b", r"\btarget\w+\b", r"\baffinity\b",
    r"\bligand\b", r"\breceptor\b", r"\bcomplex\b"
]]
MODU_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bmodulat\w+\b", r"\bregulat\w+\b", r"\baffect\w+\b", r"\binfluenc\w+\b", r"\balter\w+\b"
]]
DIRECT_HINTS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bdirectly\b", r"\bspecifically\b", r"\bselectively\b", r"\btargets?\b",
    r"\bupregulates?\b", r"\bdownregulates?\b", r"\binduces?\b", r"\bsuppresses?\b",
    r"\binhibits?\b", r"\bactivates?\b"
]]

# 유전다형성/유전자 문맥
GENETIC_ASSOC_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bpolymorphism(s)?\b", r"\bsnp(s)?\b", r"\bgenotype(s)?\b", r"\ballele(s)?\b",
    r"\bvariant(s)?\b", r"\bmutation(s)?\b"
]]
GENE_CTX_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bgene(s)?\b", r"\bmRNA\b", r"\bprotein(s)?\b", r"\bexpression\b", r"\bactivity\b",
    r"\benzyme(s)?\b", r"\breceptor(s)?\b", r"\btranscription\b", r"\bpromoter\b",
    r"\bknock[-\s]?down\b", r"\bsilenc\w+", r"\bmethylation\b", r"\bpathway(s)?\b",
    r"\bpolymorphism(s)?\b", r"\bgenotype(s)?\b", r"\bSNP(s)?\b", r"\bvariant(s)?\b"
]]
# 음료/역학 문맥(SSB 오탐 방지)
NEG_BEVERAGE_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bsugar[-\s]?sweetened\b", r"\bbeverage(s)?\b", r"\bdrinks?\b", r"\bdrinkers\b"
]]

# ---------- disease polarity ----------
DIS_POS = [re.compile(p, re.IGNORECASE) for p in [
    r"\bpromot\w+", r"\benhanc\w+", r"\bincreas\w+", r"\baccelerat\w+",
    r"\bprogression\b", r"\baggravat\w+", r"\bexacerbat\w+", r"\bmetastasi\w+",
    r"\bproliferat\w+", r"\bresistan\w+", r"\brisk\s*factor\b"
]]
DIS_NEG = [re.compile(p, re.IGNORECASE) for p in [
    r"\bsuppress\w+", r"\binhibit\w+", r"\battenuat\w+", r"\breduc\w+", r"\bprotect\w+",
    r"\bprevent\w+", r"\bblock\w+", r"\bchemo|radio\s*sensitiz\w+", r"\bgood\b.*\bprognos\w+",
    r"\bimprov\w+\s+\bsurvival\b", r"\btumor\s*suppressor\b"
]]
DIS_NEU = [re.compile(p, re.IGNORECASE) for p in [
    r"\bno\b\s+(?:significant\s+)?(?:association|effect|change)\b", r"\bnot\s+associate\w*\b",
    r"\bunchang\w+\b", r"\bnonsignificant\b", r"\bns\b"
]]

# ---------- disease terms + filters ----------
DISEASE_TERMS = [
    "disease","disorder","syndrome","cancer","tumor","tumour","carcinoma","sarcoma",
    "neoplasm","metastasis","infection","sepsis","inflammation","arthritis","asthma",
    "colitis","ulcerative colitis","ibd","hepatitis","cirrhosis","fibrosis","diabetes",
    "obesity","hypertension","stroke","ischemia","infarction","nephropathy","neuropathy",
    "retinopathy","dermatitis","psoriasis","alzheimer","parkinson","als","autism",
    "schizophrenia","depression","pneumonia","tuberculosis","hiv","aids"
]
DISEASE_PATTERNS = [(t, re.compile(rf"\b{re.escape(t)}\b", re.IGNORECASE)) for t in DISEASE_TERMS]

# 범주/통계/방법론 + 장기/검체(organ/specimen) 배제 (하드 드롭용, 모두 소문자로 보관)
NON_DISEASE_EXCLUDE = set([
    "animals","humans","materials","water","male","female","cell line","cells",
    "biomarkers","risk factors","cohort studies","case-control studies",
    "odds ratio","hazard ratio","confidence intervals",
    # 장기/검체
    "stomach","liver","kidney","brain","heart","lung","lungs","colon","prostate",
    "breast","pancreas","skin","muscle","adipose","tissue","tissues","serum",
    "plasma","urine","blood","milk","saliva","feces","faeces","hair","nails",
    "ocular","retina","cornea","bone"
])

# ---------- compound disambiguation ----------
EXCL_ENZYME_PATS = [re.compile(p, re.IGNORECASE) for p in [
    r"\blecithin\s*[:-]?\s*cholesterol\s+acyltransferase\b",   # LCAT
    r"\blecithin[-\s]?cholesterol\s+acyltransferase\b",
    r"\bLecithin\s*:\s*cholesterol\s+acyltransferase\b",
    r"\blecithin\s*[:-]?\s*retinol\s+acyltransferase\b",       # LRAT
    r"\blecithin[-\s]?retinol\s+acyltransferase\b",
]]
COFFEE_SYNS = [
    "coffee","coffee intake","coffee consumption","coffee drinkers","drinking coffee",
    "caffeinated coffee","decaffeinated coffee","decaf coffee",
    "espresso","instant coffee","filtered coffee","unfiltered coffee","brewed coffee"
]
COFFEE_EXPOSURE_RE = re.compile(r"\bcoffee\b|\bcups?\s*/?\s*day\b|\bcups?\s+per\s+day\b|\bcoffee\s+drinkers?\b", re.IGNORECASE)

# ---------- dataclasses ----------
@dataclass
class DiseaseMention:
    name: str
    mesh_id: Optional[str]=None
    source: str="text"
    normalized: str=""

@dataclass
class CompoundGeneInteraction:
    compound: str
    gene: str
    gene_display: Optional[str]
    interaction_type: str     # upregulation/downregulation/binding/modulation/genetic_association
    evidence_strength: str    # strong/moderate/weak
    context: str
    pmid: str
    confidence_score: float
    species: str = "unspecified"
    species_tags: List[str] = field(default_factory=list)
    study_context: Dict[str, Any] = field(default_factory=dict)
    evidence_tags: List[str] = field(default_factory=list)

@dataclass
class GeneDiseasePath:
    gene: str
    disease: str
    disease_id: Optional[str]
    disease_source: str
    relationship: str         # risk_factor/protective/associated
    effect_direction: str     # positive/negative/neutral (gene->disease)
    evidence_count: int
    pmids: List[str]
    species: str = "unspecified"
    rule: str = ""

@dataclass
class CompoundPathway:
    compound: str
    gene: str
    gene_display: Optional[str]
    main_interaction: str
    main_strength: str
    gene_interactions: List[CompoundGeneInteraction]
    disease_associations: List[GeneDiseasePath]
    overall_effects: Dict[str,str]     # disease -> predicted_effect
    confidence_scores: Dict[str,float] # disease -> conf
    composition_rules: Dict[str,str] = field(default_factory=dict)
    species: str = "unspecified"

@dataclass
class AnalysisResult:
    compound: str
    total_papers: int
    direct_interactions: int
    pathways: List[CompoundPathway]
    summary_stats: Dict[str, Any]
    methodology_notes: List[str]

# ---------- utils ----------
def getenv_str(name: str, default: str="") -> str:
    return (os.environ.get(name, default) or "").strip()

def ncbi_delay() -> float:
    return DEFAULT_DELAY_WITH_KEY if getenv_str("NCBI_API_KEY") else DEFAULT_DELAY_NO_KEY

def robust_request(url: str, method: str="GET", params: Optional[Dict[str,str]]=None,
                   data: Optional[Dict[str,str]]=None, headers: Optional[Dict[str,str]]=None,
                   timeout: int=30, max_retries: int=3) -> requests.Response:
    params = params or {}; data = data or {}
    headers = headers or {"User-Agent":"CompoundAnalyzer/5.8 (+https://example.org)"}
    base, backoff = 1.0, 2.0; last_exc=None
    for attempt in range(max_retries):
        try:
            if method.upper()=="POST":
                r = requests.post(url, data=data, headers=headers, timeout=timeout)
            else:
                r = requests.get(url, params=params, headers=headers, timeout=timeout)
            if r.status_code==200:
                time.sleep(0.30); return r
            if r.status_code==429:
                wait=base*(backoff**attempt)*2; tqdm.write(f"[WARN] Rate limited, wait {wait:.1f}s"); time.sleep(wait); continue
            if r.status_code>=500:
                wait=base*(backoff**attempt); tqdm.write(f"[WARN] Server {r.status_code}, retry {wait:.1f}s"); time.sleep(wait); continue
            r.raise_for_status()
        except requests.RequestException as e:
            last_exc=e; wait=base*(backoff**attempt); tqdm.write(f"[WARN] Request failed (try {attempt+1}): {e}"); time.sleep(wait)
    if last_exc: raise last_exc
    raise RuntimeError(f"Failed request to {url} after {max_retries} attempts")

# ---------- PubMed ----------
MONTH_MAP={'Jan':1,'Feb':2,'Mar':3,'Apr':4,'May':5,'Jun':6,'Jul':7,'Aug':8,'Sep':9,'Oct':10,'Nov':11,'Dec':12}
def month_to_mm(x:str)->str:
    if not x: return "01"
    x=str(x).strip()
    if x.isdigit():
        try: v=max(1,min(12,int(x))); return f"{v:02d}"
        except: return "01"
    return f"{MONTH_MAP.get(x[:3].title(),1):02d}"
def day_to_dd(d:str)->str:
    if not d: return "01"
    d=str(d).strip()
    if d.isdigit():
        try: v=max(1,min(31,int(d))); return f"{v:02d}"
        except: return "01"
    return "01"
def parse_date(article: ET.Element)->str:
    def norm(y,m,d):
        y=(y or "").strip()
        if not y: return ""
        return f"{y}-{month_to_mm(m)}-{day_to_dd(d)}"
    ad = article.find(".//ArticleDate")
    if ad is not None:
        return norm(ad.findtext("Year"), ad.findtext("Month"), ad.findtext("Day"))
    pd = article.find(".//Journal/JournalIssue/PubDate")
    if pd is not None:
        return norm(pd.findtext("Year"), pd.findtext("Month"), pd.findtext("Day"))
    md = article.findtext(".//Journal/JournalIssue/PubDate/MedlineDate") or ""
    m = re.search(r"(\d{4})", md)
    return f"{m.group(1)}-01-01" if m else ""

def esearch(term:str, retmax:int, since:Optional[str], until:Optional[str], sort:str="relevance",
            cache: Optional[PubmedCache]=None)->List[str]:
    email=getenv_str("NCBI_EMAIL","anonymous@example.com")
    apikey=getenv_str("NCBI_API_KEY","")
    params={"db":"pubmed","term":term,"retmode":"json","retmax":str(retmax),"sort":sort,"email":email}
    if apikey: params["api_key"]=apikey
    if since or until:
        params["datetype"]="pdat"
        if since: params["mindate"]=since
        if until: params["maxdate"]=until
    cache_key=None
    if cache:
        payload={"term":term,"retmax":retmax,"since":since,"until":until,"sort":sort}
        cache_key=hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()
        cached=cache.fetch_esearch(cache_key)
        if cached is not None:
            return cached
    js = robust_request(f"{NCBI_BASE}/esearch.fcgi", params=params, timeout=25).json()
    ids=js.get("esearchresult",{}).get("idlist",[]) or []
    if cache and cache_key:
        cache.store_esearch(cache_key, ids)
    return ids

def efetch_details(pmids: List[str], progress: bool=True,
                   cache: Optional[PubmedCache]=None)->Dict[str,Dict]:
    email=getenv_str("NCBI_EMAIL","anonymous@example.com")
    apikey=getenv_str("NCBI_API_KEY","")
    out: Dict[str,Dict] = {}
    remaining: List[str] = []
    if cache:
        for pmid in pmids:
            cached=cache.fetch_detail(pmid)
            if cached:
                out[pmid]=cached
            else:
                remaining.append(pmid)
    else:
        remaining=list(pmids)
    batches=[remaining[i:i+200] for i in range(0,len(remaining),200)]
    pbar=tqdm(total=len(pmids), desc="Fetching details", disable=not progress)
    if cache and out:
        pbar.update(len(out))
    for batch in batches:
        params={"db":"pubmed","id":",".join(batch),"retmode":"xml","email":email}
        if apikey: params["api_key"]=apikey
        root=ET.fromstring(robust_request(f"{NCBI_BASE}/efetch.fcgi", params=params, timeout=45).text)
        for art in root.findall(".//PubmedArticle"):
            pmid = art.findtext(".//PMID") or ""
            article = art.find(".//Article")
            if not article: continue
            tnode = article.find(".//ArticleTitle")
            title = "".join(tnode.itertext()) if tnode is not None else ""
            parts=[]
            for ab in article.findall(".//Abstract/AbstractText"):
                lbl=ab.get("Label"); txt="".join(ab.itertext()).strip()
                parts.append(f"{lbl}: {txt}" if lbl else txt)
            abstract=" ".join(p for p in parts if p)
            journal=(article.findtext(".//Journal/Title") or "").strip()
            pub_date=parse_date(article)
            mesh_terms=[]
            mesh_detail=[]
            for mh in art.findall(".//MeshHeading/DescriptorName"):
                tx=(mh.text or "").strip()
                if tx:
                    mesh_terms.append(tx)
                    mesh_detail.append({
                        "name": tx,
                        "ui": mh.get("UI"),
                        "is_major": (mh.get("MajorTopicYN") == "Y")
                    })
            pub_types=[]
            for pt in art.findall(".//PublicationType"):
                if pt.text:
                    pub_types.append(pt.text.strip())
            keywords=[]
            for kw in art.findall(".//Keyword"):
                if kw.text:
                    keywords.append(kw.text.strip())
            primary_species, species_tags = infer_species(mesh_terms)
            meta={
                "pmid":pmid,
                "title":title,
                "abstract":abstract,
                "journal":journal,
                "pub_date":pub_date,
                "mesh_terms":mesh_terms,
                "mesh_detail":mesh_detail,
                "publication_types":pub_types,
                "keywords":keywords,
                "primary_species":primary_species,
                "species_tags":species_tags
            }
            out[pmid]=meta
            if cache:
                cache.store_detail(pmid, meta)
        pbar.update(len(batch))
    pbar.close()
    return out

def sentence_split(text:str)->List[str]:
    text=(text or "").strip()
    if not text: return []
    sents=re.split(SENT_SPLIT,text)
    return [s.strip() for s in sents if s and len(s.strip())>2]

def _normalize_token(text:str)->str:
    return re.sub(r"[^a-z0-9]+","_", text.lower()).strip("_")

def infer_species(mesh_terms: List[str])->Tuple[str,List[str]]:
    found=[]
    for term in mesh_terms or []:
        low=term.lower()
        if low in SPECIES_MAP:
            display=SPECIES_MAP[low][1]
            found.append(display)
    # dedupe while preserving order
    unique=[]
    for f in found:
        if f not in unique:
            unique.append(f)
    if not unique:
        return ("unspecified", [])
    for priority in SPECIES_PRIORITY:
        if priority in unique:
            return (priority, unique)
    return (unique[0], unique)

def extract_dose_mentions(texts: List[str])->List[str]:
    hits=set()
    for t in texts:
        for m in DOSE_PATTERN.finditer(t):
            hits.add(m.group(0).strip())
    return sorted(hits)

def extract_time_mentions(texts: List[str])->List[str]:
    hits=set()
    for t in texts:
        for m in TIME_PATTERN.finditer(t):
            hits.add(m.group(0).strip())
    return sorted(hits)

def phrase_pattern(text:str)->re.Pattern:
    lower=text.lower().strip()
    if not lower:
        return re.compile(r"$^")
    escaped=re.escape(lower).replace(r"\ ", r"\s+")
    return re.compile(rf"\b{escaped}\b")

# ---------- HGNC ----------
def ensure_dir(p:str): os.makedirs(p,exist_ok=True)
def download_hgnc(cache_dir:str="cache", refresh_days:int=30)->str:
    ensure_dir(cache_dir)
    path=os.path.join(cache_dir,"hgnc_complete_set.txt")
    need=True
    if os.path.exists(path):
        mtime=dt.datetime.fromtimestamp(os.path.getmtime(path))
        need=(dt.datetime.now()-mtime)>dt.timedelta(days=refresh_days)
    if not need: return path
    for url in [HGNC_TSV_GCS, HGNC_TSV_FTP, HGNC_TSV_FTP.replace("https://","http://")]:
        try:
            tqdm.write(f"[INFO] Downloading HGNC from {url}")
            txt=robust_request(url, timeout=90).text
            with open(path,"w",encoding="utf-8") as f: f.write(txt)
            return path
        except Exception as e:
            tqdm.write(f"[WARN] HGNC source failed: {e}")
    raise RuntimeError("Failed to obtain HGNC TSV")

def build_hgnc_alias_index(tsv_path:str):
    alias_index: Dict[str, Dict[str,str]] = {}
    symbol_meta: Dict[str, Dict[str,str]] = {}
    def norm_forms(s:str)->List[str]:
        if not s: return []
        x=s.upper().replace("–","-").replace("—","-").replace("−","-")
        x=re.sub(r"[^A-Z0-9\-]","",x)
        if not x: return []
        return list({x, x.replace("-","")})
    with open(tsv_path,"r",encoding="utf-8") as f:
        reader=csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym=(row.get("symbol") or row.get("Approved symbol") or "").strip()
            if not sym: continue
            symU=sym.upper()
            hgnc_id=(row.get("hgnc_id") or row.get("HGNC ID") or row.get("HGNC ID(supplied by HGNC)") or "").strip()
            ensg=(row.get("ensembl_gene_id") or row.get("Ensembl gene ID") or "").strip()
            symbol_meta[symU]={"symbol":symU,"hgnc_id":hgnc_id,"ensembl_id":ensg}
            syns=[sym]
            for k in ["alias_symbol","Alias symbols","previous symbols","prev_symbol"]:
                val=(row.get(k) or "").strip()
                if val: syns.extend([x.strip() for x in val.split("|") if x.strip()])
            for s in set(syns):
                for k in norm_forms(s):
                    alias_index.setdefault(k, symbol_meta[symU])
    tqdm.write(f"[INFO] HGNC lexicon ready: {len(alias_index):,} keys for {len(symbol_meta):,} symbols.")
    return alias_index, symbol_meta

def match_hgnc(token:str, alias_index:Dict[str,Dict[str,str]]):
    forms=[token.upper(), token.upper().replace("-","")]
    for k in forms:
        if k in alias_index: return alias_index[k]
    return None

def scan_genes_hgnc(text:str, alias_index:Dict[str,Dict[str,str]])->List[Dict[str,str]]:
    if not text: return []
    candidates=set(re.findall(r"\b[A-Za-z][A-Za-z0-9\-]{1,15}\b", text))
    metas, seen=[], set()
    def is_gene_like(t:str)->bool:
        if t.upper()==t and len(t)>=3: return True
        if "-" in t and t.upper()==t: return True
        if any(c.isdigit() for c in t) and t.upper()==t: return True
        return False
    for tok in candidates:
        if not is_gene_like(tok): continue
        meta=match_hgnc(tok, alias_index)
        if not meta: continue
        sym=meta["symbol"]
        if sym not in seen:
            seen.add(sym); metas.append(meta)
    return metas

# ---------- disease extraction ----------
SUFFIX_DISEASE_RE = re.compile(r"\b\w+(?:itis|opathy|oma|omas)\b", re.IGNORECASE)

def extract_diseases(title:str, abstract:str, mesh_terms:List[str], mesh_detail:List[Dict[str,Any]])->List[DiseaseMention]:
    """MeSH + 키워드 기반 질병 후보 추출 (장기/검체/범주어 하드드롭 포함)"""
    text=f"{title} {abstract}"
    non_disease_norm={_normalize_token(ex) for ex in NON_DISEASE_EXCLUDE}
    found: Dict[str,DiseaseMention]={}
    # MeSH 기반
    for m in mesh_detail or []:
        ml=(m.get("name") or "").strip()
        if not ml: continue
        low=ml.lower()
        if low in NON_DISEASE_EXCLUDE:
            continue
        if (
            any(k in low for k in [
                "disease","disorder","syndrome","cancer","tumor","tumour",
                "carcinoma","sarcoma","neoplasm","fibrosis","hepatitis",
                "diabetes","infection","alzheimer","parkinson","psoriasis","colitis","retinopathy","neuropathy","nephropathy"
            ]) or SUFFIX_DISEASE_RE.search(low)
        ):
            norm=_normalize_token(ml)
            found.setdefault(norm, DiseaseMention(name=ml, mesh_id=m.get("ui"), source="mesh", normalized=norm))
    # 키워드 기반
    for term, pat in DISEASE_PATTERNS:
        if pat.search(text):
            if term.lower() not in NON_DISEASE_EXCLUDE:
                norm=_normalize_token(term)
                found.setdefault(norm, DiseaseMention(name=term, mesh_id=None, source="text", normalized=norm))
    # 접미 기반
    for token in re.findall(r"\b[A-Za-z][A-Za-z\-]{3,}\b", text):
        low=token.lower()
        if low in NON_DISEASE_EXCLUDE:
            continue
        if SUFFIX_DISEASE_RE.search(low):
            norm=_normalize_token(token)
            found.setdefault(norm, DiseaseMention(name=token, mesh_id=None, source="suffix", normalized=norm))
    # 최종 하드드롭
    ret=[]
    for dm in found.values():
        if dm.normalized in non_disease_norm:
            continue
        ret.append(dm)
    ret.sort(key=lambda d: d.name.lower())
    return ret

def classify_disease_polarity(snippet:str)->str:
    t=snippet.lower()
    pos=any(p.search(t) for p in DIS_POS)
    neg=any(p.search(t) for p in DIS_NEG)
    neu=any(p.search(t) for p in DIS_NEU)
    if neu and not pos and not neg: return "neutral"
    if pos and not neg: return "positive"
    if neg and not pos: return "negative"
    return "neutral"

# ---------- compound detection ----------
INTERVENTION_HINTS = [re.compile(p, re.IGNORECASE) for p in [
    r"\btreated\s+with\b", r"\badministered\b", r"\bsupplement(?:ed|ation)\b",
    r"\bdiet(?:ary)?\b", r"\bintake\b", r"\bconsumption\b", r"\bconsume[sd]?\b",
    r"\bdrinking\b", r"\bdrinkers\b", r"\bdrinks?\b",
    r"\bfed\b", r"\bfeeding\b", r"\bdose\b", r"\bmg/kg\b", r"\bexposed\s+to\b",
    r"\bco[-\s]?treated\b", r"\bin\s+vitro\s+with\b", r"\bin\s+vivo\s+with\b",
    r"\bhabitual\b", r"\bcups?\s*/\s*day\b", r"\bcups?\s+per\s+day\b", r"\bper\s+day\b"
]]
def build_compound_detector(compound:str, allow_ambiguous_syn:bool=False, include_caffeine_proxy:bool=False):
    c=compound.strip().lower()
    syns={c}
    if c=="coffee":
        syns.update(COFFEE_SYNS)
        if include_caffeine_proxy:
            syns.update(["caffeine"])
    elif c in {"lecithin","phosphatidylcholine"}:
        syns.update(["lecithins","phosphatidylcholine","phosphatidylcholines","pc"])
        if not allow_ambiguous_syn: syns.discard("pc")
    elif c=="caffeine":
        syns.update(["1,3,7-trimethylxanthine"])
    pos_pats=[re.compile(rf"\b{re.escape(s)}\b", re.IGNORECASE) for s in syns]
    inter_pats=INTERVENTION_HINTS
    def is_compound(s:str, require_intervention:bool=True)->bool:
        if not s: return False
        if any(p.search(s) for p in EXCL_ENZYME_PATS): return False  # LCAT/LRAT 등
        if not any(p.search(s) for p in pos_pats): return False
        if c=="caffeine" and COFFEE_EXPOSURE_RE.search(s) and not re.search(r"\bcaffein", s, re.IGNORECASE):
            return False
        if require_intervention and not any(p.search(s) for p in inter_pats): return False
        return True
    return is_compound

# ---------- analyzer ----------
class CompoundGeneAnalyzer:
    def __init__(self, progress: bool=True, require_intervention: bool=True,
                 allow_ambiguous_syn: bool=False, window:int=0, human_only:bool=False,
                 coffee_include_caffeine: bool=False):
        self.progress=progress
        self.alias_index=None
        self.require_intervention=require_intervention
        self.allow_ambiguous_syn=allow_ambiguous_syn
        self.window=max(0,int(window))
        self.human_only=human_only
        self.coffee_include_caffeine=coffee_include_caffeine

    def analyze(self, compound:str, max_papers:int=1000, include_metabolites:bool=False,
                min_evidence_threshold:int=2, since:Optional[str]=None, until:Optional[str]=None,
                sort:str="relevance", hgnc_cache:str="cache",
                pubmed_cache: Optional[PubmedCache]=None, min_confidence: float=0.0)->AnalysisResult:
        tqdm.write(f"[INFO] Starting analysis for: {compound}")
        # HGNC
        try:
            tsv=download_hgnc(cache_dir=hgnc_cache)
            self.alias_index,_=build_hgnc_alias_index(tsv)
        except Exception as e:
            tqdm.write(f"[WARN] HGNC unavailable -> regex-only fallback. {e}")
            self.alias_index=None
        # queries
        for_query=self._build_queries(compound, include_metabolites)
        all_pmids=[]
        per_q=max(50, max_papers // max(1,len(for_query)))
        for qtype, q in for_query.items():
            pmids=esearch(q, retmax=per_q, since=since, until=until, sort=sort, cache=pubmed_cache)
            tqdm.write(f"[INFO] {qtype}: {len(pmids)} PMIDs")
            all_pmids.extend(pmids); time.sleep(ncbi_delay())
        pmids=list(dict.fromkeys(all_pmids))
        if not pmids:
            return AnalysisResult(compound,0,0,[],{},["No PMIDs found."])
        details=efetch_details(pmids, progress=self.progress, cache=pubmed_cache)
        if self.human_only:
            before=len(details)
            details={k:v for k,v in details.items() if "Humans" in (v.get("mesh_terms") or [])}
            tqdm.write(f"[INFO] Human-only filter: {before} -> {len(details)} records")
        interactions=self._extract_cgi(compound, details)
        tqdm.write(f"[INFO] Candidate interactions: {len(interactions)} (unique genes: {len(set(i.gene for i in interactions))})")
        before_filter=len(interactions)
        if min_confidence>0:
            interactions=[i for i in interactions if i.confidence_score>=min_confidence]
            dropped=before_filter-len(interactions)
            if dropped>0:
                tqdm.write(f"[INFO] Confidence filter {min_confidence:.2f}: dropped {dropped} interactions.")
        gene_disease=self._infer_gene_disease(details, interactions)
        # 장기/검체 안전장치(경로 직전 하드 필터)
        non_disease_norm={_normalize_token(ex) for ex in NON_DISEASE_EXCLUDE}
        gene_disease=[gd for gd in gene_disease if _normalize_token(gd.disease) not in non_disease_norm]
        tqdm.write(f"[INFO] Gene-disease pairs before threshold: {len(gene_disease)}")
        pathways=self._build_pathways(compound, interactions, gene_disease, min_evidence_threshold)
        result=AnalysisResult(
            compound=compound,
            total_papers=len(details),
            direct_interactions=sum(1 for i in interactions if i.evidence_strength=="strong"),
            pathways=pathways,
            summary_stats=self._stats(
                interactions,
                gene_disease,
                cache_metrics=pubmed_cache.snapshot() if pubmed_cache else None,
                min_confidence=min_confidence,
                pre_filter_interactions=before_filter
            ),
            methodology_notes=self._notes()
        )
        return result

    def _build_queries(self, compound:str, include_metabolites:bool)->Dict[str,str]:
        c_raw=compound.strip()
        c_low=c_raw.lower()
        def quoted_ta(term:str)->str:
            return f"({term})" if "[" in term else f"\"{term}\"[Title/Abstract]"
        if "coffee" in c_low and (" or " not in c_low and " and " not in c_low):
            syns=COFFEE_SYNS.copy()
            base="(" + " OR ".join([quoted_ta(s) for s in syns]) + ")"
            if self.coffee_include_caffeine:
                base="(" + base + " OR \"caffeine\"[Title/Abstract])"
        else:
            base=quoted_ta(c_raw)
        qs={
            "direct_gene_interaction": f'{base} AND ( "gene expression"[Title/Abstract] OR upregulation OR downregulation OR "gene activity"[Title/Abstract] )',
            "molecular_mechanism":     f'{base} AND ( "molecular mechanism"[Title/Abstract] OR pathway OR "target gene"[Title/Abstract] OR "gene regulation"[Title/Abstract] )',
            "pharmacogenomics":        f'{base} AND ( pharmacogenomics OR "drug target"[Title/Abstract] OR "gene polymorphism"[Title/Abstract] )'
        }
        if include_metabolites:
            qs["metabolite_effects"]=f'({base} OR "metabolite of {c_raw}"[Title/Abstract]) AND gene'
        return qs

    def _extract_cgi(self, compound:str, details:Dict[str,Dict])->List[CompoundGeneInteraction]:
        """유전자-화합물 상호작용 추출(상/하향 부호 안정화)"""
        out: List[CompoundGeneInteraction]=[]
        it=tqdm(details.values(), desc="Extracting interactions", disable=not self.progress)
        include_caff=("coffee" in compound.lower())
        is_compound = build_compound_detector(
            "coffee" if "coffee"==compound.lower() else compound,
            allow_ambiguous_syn=self.allow_ambiguous_syn,
            include_caffeine_proxy=include_caff and self.coffee_include_caffeine
        )
        for meta in it:
            title, abstract = meta["title"], meta["abstract"]
            text=f"{title} {abstract}"
            sents=sentence_split(text)
            sents_l=[s.lower() for s in sents]
            primary_species = meta.get("primary_species","unspecified")
            species_tags = meta.get("species_tags",[])
            pub_types = meta.get("publication_types",[])
            keywords = meta.get("keywords",[])
            mesh_terms = meta.get("mesh_terms",[])
            # genes
            if self.alias_index:
                gene_meta_list=scan_genes_hgnc(text, self.alias_index)
                genes=[m["symbol"] for m in gene_meta_list]
                gene_meta_map={m["symbol"]:m for m in gene_meta_list}
            else:
                genes=sorted(set(re.findall(r"\b[A-Z][A-Z0-9\-]{1,15}\b", text)))
                gene_meta_map={}
            for g in genes:
                gl=g.lower()
                g_re=re.compile(rf"\b{re.escape(gl)}\b")
                votes=Counter(); direct_hits=0; contexts=[]
                evidence_tags=[]
                for j, sl in enumerate(sents_l):
                    if not g_re.search(sl): continue
                    # 유전자 문맥(없으면 스킵; 단 CYP1A2/ADORA2A는 예외적으로 약하게 허용)
                    gene_ctx = any(p.search(sl) for p in GENE_CTX_PATS)
                    if not gene_ctx and g.upper() not in {"CYP1A2","ADORA2A"}:
                        continue
                    # compound within window
                    window_idx=range(max(0, j-self.window), min(len(sents_l), j+self.window+1))
                    comp_k=None
                    for k in window_idx:
                        if is_compound(sents_l[k], require_intervention=self.require_intervention):
                            comp_k=k; break
                    if comp_k is None:
                        continue
                    same_sentence = (comp_k==j)
                    # [v5.6] 동일 문장 가중치 1.5, 이웃 0.5
                    base = 1.5 if same_sentence else 0.5
                    # direct hints
                    direct_here=any(p.search(sl) for p in DIRECT_HINTS)
                    direct_win = any(p.search(sents_l[k]) for k in window_idx for p in DIRECT_HINTS)
                    if direct_here or direct_win:
                        base *= 1.3; direct_hits+=1
                        evidence_tags.append("direct_hint")
                    # [v5.6] 방향 키워드에 우선 가중치
                    if any(p.search(sl) for p in UP_PATS):   votes["upregulation"]   += 1.2*base
                    if any(p.search(sl) for p in DOWN_PATS): votes["downregulation"] += 1.2*base
                    if any(p.search(sl) for p in BIND_PATS): votes["binding"]        += 0.8*base
                    if any(p.search(sl) for p in MODU_PATS): votes["modulation"]     += 0.6*base
                    if any(p.search(sl) for p in GENETIC_ASSOC_PATS): votes["genetic_association"] += 0.9*(0.9 if same_sentence else 0.6)
                    # SSB beverage 보호
                    if g.upper()=="SSB":
                        if any(p.search(sents_l[k]) for k in window_idx for p in NEG_BEVERAGE_PATS):
                            votes.clear(); break
                    # 컨텍스트 저장
                    if len(contexts)<3:
                        ctx=sents[j]
                        if not same_sentence:
                            ctx = ctx + " || " + sents[comp_k]
                        contexts.append(ctx)
                if not votes: continue
                itype, score = max(votes.items(), key=lambda x:x[1])
                # 강도 추정(가중 합 기반)
                if direct_hits>=2 or (itype in {"upregulation","downregulation"} and score>=2.8):
                    strength, conf = "strong", min(1.0, score/4.0)
                elif score>=1.6:
                    strength, conf = "moderate", min(1.0, score/3.0)
                else:
                    strength, conf = "weak", 0.40
                ctx=" | ".join(contexts)
                gene_meta=gene_meta_map.get(g, match_hgnc(g, self.alias_index) if self.alias_index else None)
                gene_symbol = gene_meta["symbol"] if gene_meta else g
                if not gene_meta and primary_species!="unspecified":
                    gene_display=f"{gene_symbol} ({primary_species})"
                else:
                    gene_display=gene_symbol
                dose_mentions=extract_dose_mentions(contexts)
                time_mentions=extract_time_mentions(contexts)
                study_context={
                    "species": species_tags or ([primary_species] if primary_species!="unspecified" else []),
                    "publication_types": pub_types,
                    "mesh_terms": mesh_terms[:15],
                    "keywords": keywords[:12],
                    "gene_raw": g
                }
                if dose_mentions: study_context["dose_mentions"]=dose_mentions
                if time_mentions: study_context["time_mentions"]=time_mentions
                if compound.lower()=="caffeine" and any("coffee" in kw.lower() for kw in keywords):
                    study_context.setdefault("alerts", []).append("coffee_keyword")
                out.append(CompoundGeneInteraction(
                    compound=compound,
                    gene=gene_display,
                    gene_display=gene_display,
                    interaction_type=itype,
                    evidence_strength=strength,
                    context=ctx,
                    pmid=meta["pmid"],
                    confidence_score=round(conf,3),
                    species=primary_species,
                    species_tags=species_tags,
                    study_context=study_context,
                    evidence_tags=sorted(set(evidence_tags))
                ))
        # de-dup
        uniq={}
        for i in out:
            key=(i.pmid, i.compound, i.gene, i.interaction_type)
            if key not in uniq or uniq[key].confidence_score < i.confidence_score:
                uniq[key]=i
        return sorted(uniq.values(), key=lambda x:x.confidence_score, reverse=True)

    def _infer_gene_disease(self, details:Dict[str,Dict], interactions: List[CompoundGeneInteraction])->List[GeneDiseasePath]:
        pmid2meta={v["pmid"]:v for v in details.values()}
        per_gene: Dict[str, Dict[str, Dict[str,Any]]] = defaultdict(
            lambda: defaultdict(lambda: {"count":0,"pmids":set(),"votes":[],"mention":None,"species_counter":Counter(),"pub_types":Counter()})
        )
        for i in interactions:
            meta=pmid2meta.get(i.pmid)
            if not meta: continue
            diseases=extract_diseases(meta["title"], meta["abstract"], meta.get("mesh_terms",[]), meta.get("mesh_detail",[]))
            if not diseases: continue
            sentences=sentence_split(f'{meta["title"]} {meta["abstract"]}')
            sentences_l=[s.lower() for s in sentences]
            gene_pat=phrase_pattern(i.gene_display or i.gene)
            for d in diseases:
                if not d.name: continue
                d_pat=phrase_pattern(d.name)
                sent_hits=0; pol_votes=[]
                for orig, lower in zip(sentences, sentences_l):
                    if gene_pat.search(lower) and d_pat.search(lower):
                        pol_votes.append(classify_disease_polarity(orig)); sent_hits+=1
                if sent_hits>0:
                    key=d.normalized or _normalize_token(d.name)
                    bucket=per_gene[i.gene][key]
                    bucket["count"] += 1
                    bucket["pmids"].add(i.pmid)
                    bucket["votes"].extend(pol_votes or ["neutral"])
                    bucket["mention"]=d
                    bucket["species_counter"].update([i.species or meta.get("primary_species","unspecified") or "unspecified"])
                    bucket["pub_types"].update(meta.get("publication_types",[]) or ["unspecified"])
        out=[]
        for g, dmap in per_gene.items():
            for _, obj in dmap.items():
                mention=obj.get("mention")
                if not mention:
                    continue
                cnt=obj["count"]; pmids=sorted(obj["pmids"])
                pol=Counter(obj["votes"])
                if pol:
                    eff = "positive" if pol["positive"]>pol["negative"] else ("negative" if pol["negative"]>pol["positive"] else "neutral")
                else:
                    eff = "neutral"
                rel = "risk_factor" if eff=="positive" else ("protective" if eff=="negative" else "associated")
                species_counter=obj.get("species_counter", Counter())
                species = species_counter.most_common(1)[0][0] if species_counter else "unspecified"
                out.append(GeneDiseasePath(
                    gene=g,
                    disease=mention.name,
                    disease_id=mention.mesh_id,
                    disease_source=mention.source,
                    relationship=rel,
                    effect_direction=eff,
                    evidence_count=cnt,
                    pmids=pmids[:3],
                    species=species,
                    rule=""
                ))
        return sorted(out, key=lambda x:(x.evidence_count, {"positive":2,"neutral":1,"negative":2}[x.effect_direction]), reverse=True)

    def _pick_main_interaction(self, ints: List[CompoundGeneInteraction])->Tuple[str,str]:
        if not ints: return ("modulation","weak")
        w_map={"strong":3.0,"moderate":2.0,"weak":1.0}
        score=Counter()
        for i in ints:
            score[i.interaction_type]+=w_map.get(i.evidence_strength,1.0)
        main = score.most_common(1)[0][0]
        strengths=["weak","moderate","strong"]
        best = "weak"
        for i in ints:
            if i.interaction_type==main and strengths.index(i.evidence_strength) > strengths.index(best):
                best = i.evidence_strength
        return main, best

    def _build_pathways(self, compound:str, interactions: List[CompoundGeneInteraction],
                        gene_disease: List[GeneDiseasePath], min_evidence:int)->List[CompoundPathway]:
        by_gene_int: Dict[str, List[CompoundGeneInteraction]] = defaultdict(list)
        for i in interactions: by_gene_int[i.gene].append(i)
        by_gene_dis: Dict[str, List[GeneDiseasePath]] = defaultdict(list)
        for gd in gene_disease:
            if gd.evidence_count>=min_evidence:
                by_gene_dis[gd.gene].append(gd)
        pathways=[]
        human_alias={"human","Human","Homo sapiens"}
        for g, ints in by_gene_int.items():
            if g not in by_gene_dis: continue
            dis_list=by_gene_dis[g]
            main, main_strength = self._pick_main_interaction(ints)
            species_counter=Counter(i.species for i in ints if i.species)
            primary_species=species_counter.most_common(1)[0][0] if species_counter else "unspecified"
            interaction_types=set(i.interaction_type for i in ints)
            gene_display=ints[0].gene_display if ints else g
            overall={}; confs={}; comp_rules={}

            def predict_effect(main_type:str, relationship:str, species_label:str, interaction_set:set, gd_obj:GeneDiseasePath)->Tuple[str,str]:
                base_rule=f"compose({main_type},{relationship})"
                if gd_obj.effect_direction=="neutral":
                    pred="uncertain"
                elif main_type=="genetic_association" or relationship=="associated":
                    pred="association_only"
                elif species_label and species_label not in human_alias and species_label!="unspecified":
                    pred="context_specific"
                elif {"upregulation","downregulation"}.issubset(interaction_set):
                    pred="context_specific"
                elif main_type=="upregulation":
                    pred="potentially_harmful" if relationship=="risk_factor" else ("potentially_beneficial" if relationship=="protective" else "uncertain")
                elif main_type=="downregulation":
                    pred="potentially_beneficial" if relationship=="risk_factor" else ("potentially_harmful" if relationship=="protective" else "uncertain")
                else:
                    pred="uncertain"
                if species_label and species_label not in {"unspecified"}:
                    rule=f"{base_rule} [{species_label}] -> {pred}"
                else:
                    rule=f"{base_rule} -> {pred}"
                return pred, rule

            for gd in dis_list:
                if gd.species=="unspecified" and primary_species!="unspecified":
                    gd.species=primary_species
                pred, rule = predict_effect(main, gd.relationship, gd.species or primary_species, interaction_types, gd)
                overall[gd.disease] = pred
                gd.rule = rule
                comp_rules[gd.disease]=rule
                max_conf=max(i.confidence_score for i in ints) if ints else 0.5
                dis_conf=min(1.0, gd.evidence_count/5.0)
                conf_val=round((max_conf+dis_conf)/2,3)
                if overall[gd.disease] in {"context_specific","association_only"}:
                    conf_val=round(conf_val*0.8,3)
                confs[gd.disease]=conf_val
            pathways.append(CompoundPathway(
                compound=compound,
                gene=g,
                gene_display=gene_display,
                main_interaction=main,
                main_strength=main_strength,
                gene_interactions=ints,
                disease_associations=dis_list,
                overall_effects=overall,
                confidence_scores=confs,
                composition_rules=comp_rules,
                species=primary_species
            ))
        return sorted(pathways, key=lambda p:(max(p.confidence_scores.values()) if p.confidence_scores else 0), reverse=True)

    def _stats(self, interactions: List[CompoundGeneInteraction], gene_disease: List[GeneDiseasePath],
               cache_metrics: Optional[Dict[str,int]]=None, min_confidence: float=0.0,
               pre_filter_interactions: Optional[int]=None)->Dict[str,Any]:
        data = {
            "total_gene_interactions": len(interactions),
            "strong_interactions": sum(1 for i in interactions if i.evidence_strength=="strong"),
            "affected_genes": len(set(i.gene for i in interactions)),
            "associated_diseases": len(set(gd.disease for gd in gene_disease)),
            "interaction_types": dict(Counter(i.interaction_type for i in interactions)),
            "avg_confidence": round(stats.mean([i.confidence_score for i in interactions]) if interactions else 0.0, 3)
        }
        if pre_filter_interactions is not None:
            data["pre_filter_interactions"]=pre_filter_interactions
            data["filtered_out_interactions"]=max(0, pre_filter_interactions-len(interactions))
        if min_confidence>0:
            data["min_confidence_threshold"]=round(min_confidence,3)
        if cache_metrics:
            data["cache"]=cache_metrics
        species_counter=Counter(i.species for i in interactions if i.species and i.species!="unspecified")
        if species_counter:
            data["species_annotations"]=dict(species_counter.most_common(6))
        pub_counter=Counter()
        for i in interactions:
            pub_counter.update(pt for pt in (i.study_context.get("publication_types") or []) if pt)
        if pub_counter:
            data["publication_types"]=dict(pub_counter.most_common(6))
        alert_counter=Counter()
        for i in interactions:
            alert_counter.update(i.study_context.get("alerts", []))
        if alert_counter:
            data["alerts"]=dict(alert_counter)
        return data

    def _notes(self)->List[str]:
        return [
            "Normalizes human genes with the HGNC lexicon, including previous symbols and aliases.",
            "Filters out LCAT/LRAT style 'lecithin ... acyltransferase' enzyme contexts.",
            "Treats coffee exposures via consumption verbs (consumption, drinking, cups/day, etc.).",
            "When --coffee-include-caffeine is set, caffeine is added as a proxy analyte.",
            "Sentence window (--window N): 0 keeps same sentences, 1 also keeps neighbor sentences.",
            "Gene->disease polarity: positive=pro-disease, negative=protective, neutral=mixed or unclear.",
            "Primary interactions rely on strength-weighted voting (Strong>Moderate>Weak).",
            "Organ/specimen plus taxonomy/statistical terms are excluded as non-diseases.",
            "[v5.6] Matches '...itis/...oma/...opathy' only as suffixes to avoid organ false positives.",
            "[v5.6] Same-sentence weighting and direction keywords stabilize up/down labels.",
            "[v5.6] Abbreviations without gene context (e.g., SSB) are not treated as genes.",
            "[v5.7] PubMed caching (auto/force refresh/disable, TTL configurable) speeds reruns.",
            "[v5.8] Stores species/design/dose context on edges and outputs Drug->Gene->Disease rules plus evidence snippets."
        ]

# ---------- reporting ----------
def prepare_pdf_fonts(font_regular_path: Optional[str], font_bold_path: Optional[str])->Tuple[Optional[str], Optional[str]]:
    """Register user-provided fonts for PDF rendering and return font aliases."""
    if not HAVE_REPORTLAB:
        if font_regular_path or font_bold_path:
            tqdm.write("[WARN] reportlab unavailable; PDF font options ignored.")
        return (None, None)

    def _register(path: Optional[str], alias_prefix: str)->Optional[str]:
        if not path:
            return None
        if not os.path.isfile(path):
            tqdm.write(f"[WARN] PDF font not found: {path}")
            return None
        alias=f"{alias_prefix}_{hashlib.sha1(path.encode('utf-8')).hexdigest()[:6]}"
        try:
            pdfmetrics.registerFont(TTFont(alias, path))
            return alias
        except Exception as e:
            tqdm.write(f"[WARN] Failed to register PDF font '{path}': {e}")
            return None

    body_alias=_register(font_regular_path, "BodyFont")
    bold_alias=_register(font_bold_path, "BodyFontBold")
    if not body_alias and bold_alias:
        body_alias=bold_alias
    if not bold_alias and body_alias:
        bold_alias=body_alias
    return (body_alias, bold_alias)

class CompoundReportGenerator:
    def __init__(self, result: AnalysisResult, pdf_font: Optional[str]=None, pdf_font_bold: Optional[str]=None):
        self.result=result
        self.pdf_font=pdf_font
        self.pdf_font_bold=pdf_font_bold or pdf_font

    def generate(self, out_dir:str="compound_analysis", page_size:str="A4")->Dict[str,str]:
        os.makedirs(out_dir, exist_ok=True)
        ts=dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        paths={}
        # CSV
        csv_path=os.path.join(out_dir, f"{self.result.compound}_paths_{ts}.csv")
        self._csv(csv_path); paths["csv"]=csv_path
        # JSON
        json_path=os.path.join(out_dir, f"{self.result.compound}_result_{ts}.json")
        with open(json_path,"w",encoding="utf-8") as f: json.dump(asdict(self.result), f, ensure_ascii=False, indent=2)
        paths["json"]=json_path
        # Graph
        if HAVE_GRAPH and self.result.pathways:
            png_path=os.path.join(out_dir, f"{self.result.compound}_network_{ts}.png")
            self._graph(png_path); paths["network"]=png_path
        # PDF
        if HAVE_REPORTLAB:
            pdf_path=os.path.join(out_dir, f"{self.result.compound}_report_{ts}.pdf")
            self._pdf(pdf_path, page_size); paths["pdf"]=pdf_path
        return paths

    def _csv(self, path:str):
        with open(path,"w",newline="",encoding="utf-8") as f:
            w=csv.writer(f)
            w.writerow([
                "Compound",
                "Gene (species)",
                "Drug->Gene",
                "Gene->Disease",
                "Drug->Disease (predicted)",
                "Rule",
                "Confidence",
                "Species",
                "Study context",
                "Evidence snippets",
                "Evidence tags"
            ])
            for p in self.result.pathways:
                gene_label=p.gene_display or p.gene
                pub_counter=Counter()
                dose_hits=set(); time_hits=set(); alerts=set()
                for inter in p.gene_interactions:
                    pub_counter.update(inter.study_context.get("publication_types") or [])
                    dose_hits.update(inter.study_context.get("dose_mentions", []))
                    time_hits.update(inter.study_context.get("time_mentions", []))
                    alerts.update(inter.study_context.get("alerts", []))
                context_bits=[]
                if pub_counter:
                    context_bits.append("pub:" + ", ".join(pt for pt,_ in pub_counter.most_common(3)))
                if dose_hits:
                    context_bits.append("dose:" + ", ".join(sorted(dose_hits)))
                if time_hits:
                    context_bits.append("time:" + ", ".join(sorted(time_hits)))
                if alerts:
                    context_bits.append("alerts:" + ", ".join(sorted(alerts)))
                study_context=" | ".join(context_bits)
                evidence_tags=Counter()
                for inter in p.gene_interactions:
                    evidence_tags.update(inter.evidence_tags)
                evidence_tag_str=", ".join(f"{k}:{v}" for k,v in evidence_tags.most_common()) if evidence_tags else ""
                for gd in p.disease_associations:
                    species=gd.species or p.species
                    evidence=[]
                    for pmid in gd.pmids:
                        snippet=""
                        for inter in p.gene_interactions:
                            if inter.pmid==pmid:
                                snippet=inter.context
                                break
                        if snippet:
                            snippet=snippet if len(snippet)<=160 else snippet[:160]+"..."
                            evidence.append(f"{pmid}: {snippet}")
                        else:
                            evidence.append(pmid)
                    evidence_str=" || ".join(evidence)
                    w.writerow([
                        p.compound,
                        gene_label,
                        f"{p.main_interaction} ({p.main_strength})",
                        f"{gd.relationship} ({gd.effect_direction})",
                        p.overall_effects.get(gd.disease,"uncertain"),
                        p.composition_rules.get(gd.disease, gd.rule),
                        p.confidence_scores.get(gd.disease,0.0),
                        species,
                        study_context,
                        evidence_str,
                        evidence_tag_str
                    ])

    def _graph(self, png:str):
        if not self.result.pathways: return
        G=nx.DiGraph(); c=self.result.compound; G.add_node(c, color="red")
        for p in self.result.pathways:
            gene_node=p.gene_display or p.gene
            G.add_node(gene_node, color="skyblue")
            G.add_edge(c, gene_node, interaction=p.main_interaction, weight=3)
            for gd in p.disease_associations:
                G.add_node(gd.disease, color="lightgreen")
                G.add_edge(gene_node, gd.disease, relation=gd.relationship, weight=max(1,gd.evidence_count))
        plt.figure(figsize=(14,10))
        pos=nx.spring_layout(G, k=2.2, iterations=60)
        colors=[G.nodes[n].get("color","gray") for n in G.nodes()]
        nx.draw(G, pos, node_color=colors, node_size=900, with_labels=True, font_size=8, arrows=True, edge_color="gray")
        plt.title(f"Compound->Gene->Disease network: {c}")
        plt.tight_layout(); plt.savefig(png, dpi=300); plt.close()

    def _pdf(self, path:str, page_size:str="A4"):
        ps=A4 if page_size.upper()=="A4" else LETTER
        doc=SimpleDocTemplate(path, pagesize=ps, rightMargin=1.5*cm,leftMargin=1.5*cm, topMargin=1.5*cm,bottomMargin=1.5*cm)
        sty=getSampleStyleSheet()
        body_font=self.pdf_font or sty["BodyText"].fontName
        bold_font=self.pdf_font_bold or (self.pdf_font or sty["Heading1"].fontName)
        H1=ParagraphStyle("H1", parent=sty["Heading1"], fontSize=18, spaceAfter=6, fontName=bold_font)
        H2=ParagraphStyle("H2", parent=sty["Heading2"], fontSize=14, spaceAfter=6, fontName=bold_font)
        small=ParagraphStyle("Small", parent=sty["BodyText"], fontSize=9, leading=12, fontName=body_font)
        body_style=ParagraphStyle("Body", parent=sty["BodyText"], fontSize=11, leading=14, fontName=body_font)
        table_text=ParagraphStyle("TableText", parent=sty["BodyText"], fontSize=9, leading=11, spaceAfter=0, spaceBefore=0, fontName=body_font)
        table_head=ParagraphStyle("TableHeader", parent=sty["Heading4"], fontSize=10, leading=12, alignment=1, fontName=bold_font)

        def cell(val: Any, header: bool=False):
            style = table_head if header else table_text
            if isinstance(val, float):
                text = f"{val:.2f}"
            else:
                text = str(val)
            text = escape(text).replace("\n", "<br/>")
            return Paragraph(text, style)

        story=[]
        story.append(Paragraph(f"Compound Analysis Report: <b>{self.result.compound}</b>", H1))
        story.append(Paragraph(f"Generated: {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", small))
        s=self.result.summary_stats
        story.append(Paragraph(
            f"<b>Total papers:</b> {self.result.total_papers} | "
            f"<b>Direct interactions (strong):</b> {self.result.direct_interactions} | "
            f"<b>Affected genes:</b> {s.get('affected_genes',0)} | "
            f"<b>Associated diseases:</b> {s.get('associated_diseases',0)} | "
            f"<b>Avg confidence:</b> {s.get('avg_confidence',0):.2f}", body_style))
        story.append(Spacer(1,0.3*cm))
        legend=("Legend - Drug->Gene / Gene->Disease / Drug->Disease columns report directionality and predicted outcomes, "
                "and the Rule column keeps the composition rule plus species tags when needed.")
        story.append(Paragraph(legend, small)); story.append(Spacer(1,0.3*cm))

        if self.result.pathways:
            story.append(Paragraph("Top Pathways", H2))
            headers=["Gene (species)","Drug->Gene","Gene->Disease","Drug->Disease","Rule","Conf.","PMIDs"]
            rows=[[cell(h, header=True) for h in headers]]
            flat=[]
            for p in self.result.pathways:
                for gd in p.disease_associations:
                    flat.append((p,gd))
            flat=sorted(flat, key=lambda x: x[0].confidence_scores.get(x[1].disease,0.0), reverse=True)[:12]
            for p, gd in flat:
                conf=p.confidence_scores.get(gd.disease,0.0)
                pred=p.overall_effects.get(gd.disease,"uncertain")
                gene_label=p.gene_display or p.gene
                species=gd.species or p.species
                row_values=[
                    f"{gene_label}" + (f" [{species}]" if species and species.lower() not in gene_label.lower() and species!="unspecified" else ""),
                    f"{p.main_interaction} ({p.main_strength})",
                    f"{gd.relationship} ({gd.effect_direction})",
                    pred,
                    p.composition_rules.get(gd.disease, gd.rule),
                    f"{conf:.2f}",
                    ", ".join(gd.pmids)
                ]
                rows.append([cell(val) for val in row_values])
            col_widths=[3.0*cm, 2.6*cm, 2.8*cm, 2.8*cm, 4.0*cm, 1.5*cm, 3.5*cm]
            tbl=Table(rows, colWidths=col_widths)
            tbl.setStyle(TableStyle([
                ("BACKGROUND",(0,0),(-1,0),colors.lightgrey),("FONTNAME",(0,0),(-1,0),"Helvetica-Bold"),
                ("GRID",(0,0),(-1,-1),0.25,colors.grey),("FONTSIZE",(0,0),(-1,0),10),("VALIGN",(0,0),(-1,-1),"TOP"),
                ("ALIGN",(5,1),(5,-1),"CENTER")
            ]))
            story.append(tbl); story.append(Spacer(1,0.3*cm))
            story.append(Paragraph("Evidence Snippets", H2))
            for p in self.result.pathways[:6]:
                gene_label=p.gene_display or p.gene
                story.append(Paragraph(f"<b>{escape(gene_label)}</b>", table_text))
                for inter in p.gene_interactions[:3]:
                    snippet=inter.context if len(inter.context)<=200 else inter.context[:200]+"..."
                    story.append(Paragraph(f"{escape(inter.pmid)} - {escape(snippet)}", small))
                story.append(Spacer(1,0.2*cm))
        else:
            story.append(Paragraph("No significant pathways identified.", body_style))
        story.append(PageBreak()); story.append(Paragraph("Methodology", H2))
        for note in self.result.methodology_notes:
            story.append(Paragraph(f"- {note}", body_style))
        doc.build(story)

# ---------- CLI ----------
def main():
    ap=argparse.ArgumentParser(description="Compound->Gene->Disease pathway analyzer (v5.8)",
                               formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("compound", help="Compound/chemical name")
    ap.add_argument("--max-papers", type=int, default=1000)
    ap.add_argument("--include-metabolites", action="store_true")
    ap.add_argument("--min-evidence", type=int, default=1, help="Min papers per gene-disease to keep pathway")
    ap.add_argument("--since", type=str, default=None, help="YYYY/MM/DD")
    ap.add_argument("--until", type=str, default=None, help="YYYY/MM/DD")
    ap.add_argument("--sort", type=str, default="relevance", choices=["relevance","date"])
    ap.add_argument("--hgnc-cache", type=str, default="cache")
    ap.add_argument("--outdir", type=str, default="compound_analysis")
    ap.add_argument("--page-size", type=str, default="A4", choices=["A4","Letter"])
    ap.add_argument("--no-progress", action="store_true")
    ap.add_argument("--pubmed-cache-dir", type=str, default=DEFAULT_CACHE_DIR,
                    help="Directory to store PubMed cache files")
    ap.add_argument("--cache-ttl-days", type=int, default=7, help="Cache expiry in days (use 0 to always reuse)")
    ap.add_argument("--cache-mode", choices=["auto","refresh","off"], default="auto",
                    help="auto=use cache within TTL, refresh=ignore cache but overwrite, off=disable caching")
    ap.add_argument("--min-confidence", type=float, default=0.0,
                    help="Drop interactions below this confidence score (0.0~1.0)")
    ap.add_argument("--pdf-font", type=str, default=None,
                    help="Path to TrueType/OpenType font for PDF body text (set for full Unicode support).")
    ap.add_argument("--pdf-font-bold", type=str, default=None,
                    help="Optional font file for PDF headings (defaults to --pdf-font).")

    # scope/synonyms/window/human
    ap.add_argument("--scope", choices=["all","intervention"], default="intervention",
                    help="Require intervention/exposure context for compound mentions")
    ap.add_argument("--allow-ambiguous-syn", action="store_true",
                    help="Allow ambiguous synonyms (e.g., 'PC' for phosphatidylcholine)")
    ap.add_argument("--window", type=int, default=0, help="Sentence window to link compound and gene (0=same sentence)")
    ap.add_argument("--human-only", action="store_true", help="Keep only MeSH 'Humans' records")

    # coffee proxy
    ap.add_argument("--coffee-include-caffeine", action="store_true",
                    help="If compound is 'coffee', also treat 'caffeine' as proxy")

    args=ap.parse_args()

    if not (0.0 <= args.min_confidence <= 1.0):
        ap.error("--min-confidence must be between 0.0 and 1.0")

    if not getenv_str("NCBI_EMAIL"):
        tqdm.write("[WARN] NCBI_EMAIL env var not set. export NCBI_EMAIL=you@org")

    body_font_name, bold_font_name = prepare_pdf_fonts(args.pdf_font, args.pdf_font_bold)

    cache=None
    if args.cache_mode!="off":
        cache=PubmedCache(base_dir=args.pubmed_cache_dir,
                          ttl_days=args.cache_ttl_days,
                          enabled=True,
                          force_refresh=(args.cache_mode=="refresh"))

    analyzer=CompoundGeneAnalyzer(progress=not args.no_progress,
                                  require_intervention=(args.scope=="intervention"),
                                  allow_ambiguous_syn=args.allow_ambiguous_syn,
                                  window=args.window,
                                  human_only=args.human_only,
                                  coffee_include_caffeine=args.coffee_include_caffeine)
    res=analyzer.analyze(compound=args.compound, max_papers=args.max_papers,
                         include_metabolites=args.include_metabolites,
                         min_evidence_threshold=args.min_evidence,
                         since=args.since, until=args.until, sort=args.sort,
                         hgnc_cache=args.hgnc_cache,
                         pubmed_cache=cache,
                         min_confidence=args.min_confidence)

    print(f"\nTotal papers: {res.total_papers}")
    print(f"Direct interactions (strong): {res.direct_interactions}")
    print(f"Pathways: {len(res.pathways)}")

    # Top 10 preview - composition rule summary
    rows=[]
    for p in res.pathways[:10]:
        if not p.disease_associations: continue
        gd=p.disease_associations[0]
        conf=p.confidence_scores.get(gd.disease,0.0)
        gene_label=p.gene_display or p.gene
        species=gd.species or p.species
        rule=p.composition_rules.get(gd.disease, gd.rule)
        pmids=";".join(gd.pmids)
        snippet=""
        if p.gene_interactions:
            snippet=p.gene_interactions[0].context
            if len(snippet)>120: snippet=snippet[:120]+"..."
        rows.append([
            p.compound,
            gene_label + (f" [{species}]" if species and species!="unspecified" and species.lower() not in gene_label.lower() else ""),
            f"{p.main_interaction} ({p.main_strength})",
            f"{gd.relationship} ({gd.effect_direction})",
            p.overall_effects.get(gd.disease,"uncertain"),
            rule,
            f"{conf:.2f}",
            pmids,
            snippet
        ])
    if rows:
        print("\nTop 10 findings:")
        print("Compound\tGene(species)\tDrug->Gene\tGene->Disease\tDrug->Disease\tRule\tConf\tPMIDs\tSnippet")
        for r in rows: print("\t".join(r))

    if cache:
        cache_stats=cache.snapshot()
        if cache_stats:
            print("\nCache stats:")
            for k,v in sorted(cache_stats.items()):
                print(f"  {k}: {v}")

    rg=CompoundReportGenerator(res, pdf_font=body_font_name, pdf_font_bold=bold_font_name)
    paths=rg.generate(out_dir=args.outdir, page_size=args.page_size)
    print("\nGenerated files:")
    for k,v in paths.items(): print(f"  {k.upper()}: {v}")
    if "pdf" not in paths: print("  (Install 'reportlab' to enable PDF export)")

if __name__=="__main__":
    main()
