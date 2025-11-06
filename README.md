# NutriGenomics 데이터 파이프라인 설치 가이드

> “5단계만 따라 하면 원본 데이터 다운로드부터 분석 리포트까지 한 번에 끝납니다.”

## 📦 핵심 스크립트
- `setup_data_pipeline.py` – FDC/CTD/GWAS 다운로드 + 전처리를 한 번에 실행.
- `download_data.py` – 필요한 데이터셋만 골라서 다시 받을 때 사용.
- `preprocess.py` – `data/raw/*`를 읽어 `data/processed/*.parquet`를 생성.
- `data_structure_analyzer.py` – 전처리 결과를 시각화해 `food_health_report.html` 출력.
- `pubmed_compound_gene_disease_analyzer.py` – PubMed 기반 화합물→유전자→질환 체인을 추출/리포팅.

## ⚙️ 사전 준비
| 항목 | 요구 사항 |
| --- | --- |
| Python | 3.10 이상 (PEP 604 타입 표기 사용) |
| OS 도구 | `git`, `curl`(또는 `wget`), `unzip`, `gzip`, `tar`, `bash/zsh` |
| 디스크 | 원본 + 전처리 파일 포함 15GB 이상 여유 |
| 네트워크 | FDC/GWAS 전체 다운로드에 수 GB가 필요하므로 안정적인 고속 회선 권장 |

### OS 별 필수 패키지 설치 예시
```bash
# Ubuntu / Debian
sudo apt-get update
sudo apt-get install -y python3.11 python3.11-venv python3-pip git curl unzip

# macOS (Homebrew)
brew install python@3.11 git wget unzip
```

## 🚀 5분 설치 절차
1. **저장소로 이동**
   ```bash
   cd /path/to/NutriGenomics
   ```
2. **가상환경 만들기**
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate        # Windows: .venv\Scripts\activate
   python -m pip install --upgrade pip setuptools wheel
   ```
3. **필수 패키지 설치**
   ```bash
   pip install -r requirements.txt
   ```
4. **NCBI API 환경변수 설정 (PubMed 파이프라인용)**
   ```bash
   export NCBI_EMAIL="you@company.org"
   export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # 선택이지만 강력 권장
   ```
5. **전체 파이프라인 실행**
   ```bash
   python setup_data_pipeline.py --all
   ```
   실행이 끝나면 `data/processed/`에 Parquet, `food_health_report.html`에 시각화 결과가 생성됩니다.

### 다음 실행(재방문)
설치가 끝난 뒤에는 매번 아래 두 줄만 실행하면 됩니다.
```bash
source .venv/bin/activate           # Windows: .venv\Scripts\activate
python setup_data_pipeline.py --preprocess-only
```
원본 데이터를 다시 내려받아야 할 때만 `--download-only` 옵션을 추가하세요.

## 🔽 데이터 다운로드 & 전처리
| 명령 | 설명 |
| --- | --- |
| `python setup_data_pipeline.py --all` | 전체 다운로드 + 전처리(기본). |
| `python setup_data_pipeline.py --download-only` | 원본 데이터만 새로 받을 때. |
| `python setup_data_pipeline.py --preprocess-only --limit-rows 5000` | 로컬에 이미 받아둔 데이터를 빠르게 전처리하거나 스모크 테스트할 때. |
| `python download_data.py --datasets ctd gwas` | 특정 데이터셋만 갱신. (FDC 용량이 커서 건너뛰고 싶을 때 사용) |

> ⚠️ 네트워크/서버 이슈로 특정 파일을 받지 못하면 스크립트가 경고를 출력합니다. 이 경우 해당 링크에서 수동으로 파일을 내려받아 `data/raw/<dataset>/`에 넣어 두면 재실행 시 자동으로 감지합니다.

## 🔍 설치 검증 & 스모크 테스트
1. **전처리 로직 점검**
   ```bash
   python preprocess.py --diagnose all --limit-rows 1000
   ```
   → 파일 스키마/결측치 정보를 출력해 환경 문제를 빠르게 잡습니다.

2. **산출물 확인**
   ```bash
   ls data/processed
   python - <<'PY'
   import pandas as pd; print(pd.read_parquet("data/processed/fdc_foods.parquet").head())
   PY
   ```

3. **리포트 생성**
   ```bash
   python data_structure_analyzer.py
   open food_health_report.html      # Windows: start food_health_report.html
   ```

## 📊 PubMed 화합물-유전자-질환 파이프라인
- 기본 캐시 위치: `cache/pubmed` (`--pubmed-cache-dir`로 변경 가능).
- 최소 요구: `NCBI_EMAIL`. API key가 있으면 속도 제한이 완화됩니다 (`pubmed_compound_gene_disease_analyzer.py` 참조).
```bash
python pubmed_compound_gene_disease_analyzer.py \
  --compound "curcumin" \
  --max-papers 80 \
  --scope intervention \
  --human-only \
  --output-dir outputs/curcumin_run
```
- 캐시 정리: `--cache-mode refresh` 또는 `rm -rf cache/pubmed`.

## 📈 운영 팁 & 트러블슈팅
- **디스크 부족**: `data/raw/`(원본)과 `cache/`(PubMed 캐시)가 대부분 차지합니다. 오래된 릴리즈는 다른 위치로 옮겨 보관하세요.
- **중간 실패**: 동일 명령을 다시 실행하면 이미 받은 파일은 건너뜁니다. 필요 시 `data/raw/<dataset>` 폴더를 지우고 재실행하세요.
- **의존성 업데이트**: `pip install -r requirements.txt --upgrade` 후 `python -m pip check`로 호환성을 검증합니다.
- **Windows PowerShell**: 명령어 앞에 `python` 대신 `py`를 사용하고, 가상환경 활성화는 `.\.venv\Scripts\Activate.ps1`.

## 📄 Requirements
`requirements.txt`에는 실행에 필요한 모든 서드파티 패키지를 정리했습니다. 새 패키지가 필요하면 해당 파일을 수정하고 `pip install -r requirements.txt`를 다시 실행하세요.

---
질문이나 개선 아이디어가 있다면 README에 Issue/PR 가이드를 추가해 주세요. 즐거운 데이터 분석 되세요! 🎉
