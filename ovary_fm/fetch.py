"""``ovary-fm fetch``: bring each enabled dataset's files to data/raw/<name>/.

Source types (config/corpus.yaml):

* ``local`` -- read in place; only checked.
* ``url`` -- plain downloads through the proxy (urllib honours https_proxy);
  a sha256 and the byte count go to download.json beside the files.
* ``cellxgene_collection`` -- the Census API (``pip install -e .[census]``):
  the datasets of one collection; with ``author_h5ad: true`` the curated
  h5ad that carries the authors' labels (Mouse_Compensatory's
  scripts/fetch_reference.py is the model). UNTESTED here: the Census
  hosts were unreachable from the machine this was written on.
* ``scbasecount`` -- per-sample h5ad files in the Arc bucket, resolved
  from the local copy of the sample metadata that
  scripts/scbasecount_metadata.py writes (``file_path`` column) and read
  over anonymous HTTPS.

On Biowulf run it as a batch job, not on the login node:
``bash biowulf/submit.sh fetch``.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import tarfile
import time
import urllib.request
from pathlib import Path

import pandas as pd

from .config import data_dir, write_provenance
from .io import dataset_path, datasets

GCS_HTTPS = "https://storage.googleapis.com/"


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: Path, overwrite: bool = False) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not overwrite:
        print(f"  = {dest.name} already on disk ({dest.stat().st_size:,} bytes)")
        return {"url": url, "path": str(dest), "bytes": dest.stat().st_size, "sha256": sha256_of(dest),
                "skipped": True}
    print(f"  > {url}")
    start = time.time()
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as handle:
        total = response.length or 0
        done = 0
        for chunk in iter(lambda: response.read(1 << 20), b""):
            handle.write(chunk)
            done += len(chunk)
            if total and done % (64 << 20) < (1 << 20):
                print(f"    {done / total:6.1%}  {done:,} / {total:,}", flush=True)
    tmp.rename(dest)
    record = {"url": url, "path": str(dest), "bytes": dest.stat().st_size, "sha256": sha256_of(dest),
              "seconds": round(time.time() - start, 1), "skipped": False}
    print(f"  < {dest.name}: {record['bytes']:,} bytes in {record['seconds']} s")
    return record


def _safe_extract(tar: tarfile.TarFile, dest: Path):
    try:
        tar.extractall(dest, filter="data")          # Python >= 3.12: no paths outside dest, no links
    except TypeError:
        members = [m for m in tar.getmembers() if not (m.name.startswith("/") or ".." in m.name.split("/"))]
        tar.extractall(dest, members=members)


def unpack(path: Path) -> Path | None:
    """A GEO RAW tar -> a folder beside it (per-sample tars inside it unpacked one level);
    an .h5ad.gz or .rds.gz -> the uncompressed file beside it. Returns what was produced."""
    name = path.name
    if name.endswith((".tar", ".tar.gz", ".tgz")):
        dest = path.parent / name.split(".tar")[0].removesuffix(".tgz")
        if dest.exists() and any(dest.iterdir()):
            print(f"  = {dest.name}/ already unpacked")
            return dest
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(path) as tar:
            _safe_extract(tar, dest)
        inner = [p for p in dest.rglob("*") if p.is_file() and p.name.endswith((".tar", ".tar.gz", ".tgz"))]
        for p in inner:
            sub = p.parent / p.name.split(".tar")[0].removesuffix(".tgz")
            sub.mkdir(exist_ok=True)
            with tarfile.open(p) as tar:
                _safe_extract(tar, sub)
        n = sum(1 for _ in dest.rglob("*") if _.is_file())
        print(f"  unpacked {name} -> {dest.name}/ ({n} files{', ' + str(len(inner)) + ' inner tars' if inner else ''})")
        return dest
    if name.endswith((".h5ad.gz", ".rds.gz", ".h5.gz")):
        out = path.with_suffix("")
        if out.exists():
            print(f"  = {out.name} already uncompressed")
            return out
        with gzip.open(path, "rb") as src, open(out, "wb") as dst:
            shutil.copyfileobj(src, dst, 1 << 24)
        print(f"  gunzipped {name} -> {out.name} ({out.stat().st_size:,} bytes)")
        return out
    return None


def _fetch_urls(entry, folder: Path) -> list[dict]:
    src = entry.get("source") or {}
    urls = src.get("urls") or []
    if not urls:
        print(f"  !! {entry['name']}: source.urls is empty (status {entry.get('status')}); "
              f"list the GEO supplementary folder and fill config/corpus.yaml")
        return []
    records = []
    for url in urls:
        rec = download(url, folder / url.rsplit("/", 1)[-1])
        if src.get("unpack", True):
            produced = unpack(Path(rec["path"]))
            if produced is not None:
                rec["unpacked_to"] = str(produced)
        records.append(rec)
    if any(str(r["path"]).endswith(".rds") for r in records):
        rds = next(r["path"] for r in records if str(r["path"]).endswith(".rds"))
        print(f"  a Seurat/SCE object: convert it before qc:\n"
              f"    module load R && Rscript R/rds_to_mtx.R {rds} {folder}/")
    return records


def _fetch_census(entry, dest: Path) -> list[dict]:
    try:
        import cellxgene_census
    except ImportError:
        raise SystemExit("cellxgene-census is not installed: pip install -e .[census]") from None
    src = entry["source"]
    collection = src.get("collection_id")
    dataset_id = src.get("dataset_id")
    version = src.get("census_version", "stable")
    if not collection and not dataset_id:
        raise SystemExit(f"{entry['name']}: a cellxgene_collection source needs collection_id or dataset_id")
    with cellxgene_census.open_soma(census_version=version) as census:
        table = census["census_info"]["datasets"].read().concat().to_pandas()
    if dataset_id:
        table = table[table["dataset_id"] == dataset_id]
    else:
        table = table[table["collection_id"] == collection]
    if table.empty:
        raise SystemExit(f"{entry['name']}: nothing in Census {version} for {collection or dataset_id}")
    # a collection may hold the same cells more than once (Isola 2024: the
    # atlas plus four lineage subsets); take the largest dataset unless told
    table = table.sort_values("dataset_total_cell_count", ascending=False)
    chosen = table.iloc[0]
    print(f"  census {version}: {len(table)} dataset(s) in the collection; taking "
          f"{chosen['dataset_id']} ({chosen['dataset_total_cell_count']:,} cells): {chosen.get('dataset_title')}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print(f"  = {dest} already on disk")
    elif src.get("author_h5ad", True):
        cellxgene_census.download_source_h5ad(chosen["dataset_id"], to_path=str(dest), census_version=version)
    else:
        adata = cellxgene_census.get_anndata(census=cellxgene_census.open_soma(census_version=version),
                                             organism="Mus musculus" if entry.get("species") == "mouse" else "Homo sapiens",
                                             obs_value_filter=f"dataset_id == '{chosen['dataset_id']}' and is_primary_data == True",
                                             X_name="raw")
        adata.write_h5ad(dest)
    table.to_csv(dest.parent / "census_datasets.csv", index=False)
    return [{"dataset_id": chosen["dataset_id"], "path": str(dest), "census_version": version,
             "sha256": sha256_of(dest) if dest.exists() else None}]


def _fetch_osf(entry, folder: Path) -> list[dict]:
    """An OSF project's files through the public API (https://api.osf.io/v2/nodes/<id>/files/osfstorage/).
    A project that asks for a sign-in answers 401/403: then download in a browser and drop the
    files under the dataset folder."""
    import urllib.error

    src = entry["source"]
    project = src.get("project")
    if not project:
        raise SystemExit(f"{entry['name']}: an osf source needs `project` (the 5-character id in the URL)")
    api = f"https://api.osf.io/v2/nodes/{project}/files/osfstorage/?page[size]=100"
    manual = (f"  !! OSF project {project} did not answer the public API; sign in at https://osf.io/{project}/files "
              f"(ORCID or institutional login), download the file(s) and put them under {folder}/, then run qc")
    try:
        with urllib.request.urlopen(api, timeout=120) as response:
            listing = json.loads(response.read().decode("utf-8"))
    except (urllib.error.HTTPError, urllib.error.URLError) as exc:
        print(f"  {api}: {exc}")
        print(manual)
        return []
    records = []
    for item in listing.get("data", []):
        attrs = item.get("attributes", {})
        if attrs.get("kind") != "file":
            continue
        url = (item.get("links") or {}).get("download")
        if not url:
            continue
        try:
            rec = download(url, folder / attrs["name"])
        except urllib.error.HTTPError as exc:
            print(f"  {attrs['name']}: {exc}")
            print(manual)
            continue
        if src.get("unpack", True):
            produced = unpack(Path(rec["path"]))
            if produced is not None:
                rec["unpacked_to"] = str(produced)
        records.append(rec)
    if not records:
        print(manual)
    return records


def _fetch_scbasecount(cfg, entry, folder: Path) -> list[dict]:
    src = entry["source"]
    species = entry.get("species", "mouse")
    meta_path = data_dir(cfg, "scbasecount") / f"{species}_sample_metadata.parquet"
    if not meta_path.exists():
        print(f"  !! {meta_path} is missing: run `python scripts/scbasecount_metadata.py --species {species}` first")
        return []
    meta = pd.read_parquet(meta_path)
    wanted = src.get("srx") or []
    if len(wanted) == 2 and all(str(w).startswith("SRX") for w in wanted) and src.get("range", True):
        lo, hi = (int(str(w)[3:]) for w in wanted)
        wanted = [f"SRX{i}" for i in range(lo, hi + 1)]
    rows = meta[meta["srx_accession"].isin(wanted)] if wanted else meta.iloc[0:0]
    if rows.empty:
        print(f"  !! none of {wanted[:3]}... in {meta_path}")
        return []
    records = []
    for _, row in rows.iterrows():
        gs = str(row["file_path"])
        url = GCS_HTTPS + gs[len("gs://"):] if gs.startswith("gs://") else gs
        records.append(download(url, folder / f"{row['srx_accession']}.h5ad"))
    rows.to_csv(folder / "scbasecount_samples.csv", index=False)
    return records


def run(cfg, dataset: str | None = None):
    write_provenance(cfg, "fetch", extra={"dataset": dataset})
    for entry in datasets(cfg, names=[dataset] if dataset else None):
        src = entry.get("source") or {}
        kind = src.get("type", "local")
        path = dataset_path(cfg, entry)
        folder = path if path.suffix == "" else path.parent
        print(f"== {entry['name']} ({kind}; {src.get('accession') or '-'})")
        if kind in ("local", "manual"):
            print(f"  {'found' if path.exists() else 'MISSING'}: {path}")
            if kind == "manual" and not path.exists():
                print(f"  download by hand from {src.get('url') or src.get('accession') or 'the source in the notes'} "
                      f"and put it under {folder}/")
            continue
        if kind == "url":
            records = _fetch_urls(entry, folder)
        elif kind == "osf":
            records = _fetch_osf(entry, folder)
        elif kind == "cellxgene_collection":
            records = _fetch_census(entry, path if path.suffix == ".h5ad" else folder / f"{entry['name']}.h5ad")
        elif kind == "scbasecount":
            records = _fetch_scbasecount(cfg, entry, folder)
        else:
            raise ValueError(f"{entry['name']}: unknown source type {kind!r}")
        if records:
            folder.mkdir(parents=True, exist_ok=True)
            (folder / "download.json").write_text(json.dumps(records, indent=2, default=str), encoding="utf-8")
