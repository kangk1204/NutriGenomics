from pathlib import Path
import zipfile
import requests
from .io import sha256, utc_now, write_json

BASE = "https://www.bindingdb.org/rwd/bind/downloads/"
FILES = {
    "all": "BindingDB_All_202609_tsv.zip",
    "articles": "BindingDB_BindingDB_Articles_202609_tsv.zip",
    "assays": "BindingDB_Assays_202609_tsv.zip",
    "assay_map": "BindingDB_rsid_eaids_202609_tsv.zip",
}


def fetch(out, subset="all", include_assays=True):
    """Download dated native archives. CRC and SHA256 are checked before use."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    names = [FILES[subset]] + ([FILES["assays"], FILES["assay_map"]] if include_assays else [])
    records = []
    for name in names:
        destination = out / name
        url = BASE + name
        if not destination.exists():
            temporary = destination.with_suffix(".partial")
            with requests.get(url, stream=True, timeout=(30, 180)) as response:
                response.raise_for_status()
                with open(temporary, "wb") as stream:
                    for chunk in response.iter_content(1024 * 1024):
                        stream.write(chunk)
            temporary.replace(destination)
        with zipfile.ZipFile(destination) as archive:
            bad = archive.testzip()
            if bad:
                raise ValueError(f"ZIP CRC failed: {bad}")
            members = [{"name": info.filename, "bytes": info.file_size, "zip_date": list(info.date_time)} for info in archive.infolist()]
        records.append({"file": name, "url": url, "bytes": destination.stat().st_size,
                        "sha256": sha256(destination), "members": members})
        print(f"verified {name}: {destination.stat().st_size:,} bytes", flush=True)
    manifest = {"retrieved_utc": utc_now(), "provider": "BindingDB", "release_filename_version": "202609",
                "official_listing_file_updated_date": "2026-08-30", "files": records,
                "license": "BindingDB curated: CC BY 3.0; ChEMBL imported: CC BY-SA 3.0; preserve row origin",
                "official_listing": "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp"}
    write_json(out / "source_manifest.json", manifest)
    return manifest
