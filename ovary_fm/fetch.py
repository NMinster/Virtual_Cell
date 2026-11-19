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

import hashlib
import json
import os
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
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


ARCHIVE_SUFFIXES = (".zip", ".tar", ".tar.gz", ".tgz", ".tar.bz2")


def url_target_name(item) -> tuple[str, str]:
    """A `source.urls` item is a URL string or `{url: ..., name: ...}`.

    The file name is `name` when given, else the URL's last non-empty path
    segment. A URL whose last segment is an opaque id (OSF `/download/dsf4r/`)
    must carry `name`, or the download is refused rather than saved nameless.
    """
    if isinstance(item, dict):
        url, name = item["url"], item.get("name")
    else:
        url, name = str(item), None
    if not name:
        segment = urllib.parse.urlparse(url).path.rstrip("/").rsplit("/", 1)[-1]
        name = urllib.parse.unquote(segment)
    if not name or "." not in name:
        raise SystemExit(f"cannot name the download from {url!r}: give it as {{url: ..., name: <file.ext>}} in corpus.yaml")
    return url, name


def _safe_members(archive, folder: Path):
    """Reject archive members that would land outside `folder`."""
    root = folder.resolve()
    names = archive.namelist() if isinstance(archive, zipfile.ZipFile) else archive.getnames()
    for name in names:
        if not (root / name).resolve().is_relative_to(root):
            raise SystemExit(f"archive member escapes {folder}: {name!r}")


def extract_archive(path: Path, folder: Path) -> dict:
    """Unpack a .zip / .tar* next to itself, once (a marker file records it)."""
    marker = folder / f".{path.name}.extracted"
    if marker.exists():
        return {"archive": path.name, "extracted": True, "skipped": True}
    print(f"  ~ extracting {path.name}")
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as z:
            _safe_members(z, folder)
            z.extractall(folder)
            n = len(z.namelist())
    else:
        with tarfile.open(path) as t:
            _safe_members(t, folder)
            try:
                t.extractall(folder, filter="data")
            except TypeError:               # python < 3.12
                t.extractall(folder)
            n = len(t.getnames())
    marker.write_text(f"{n} members\n")
    return {"archive": path.name, "extracted": True, "members": n, "skipped": False}


def _fetch_urls(entry, folder: Path) -> list[dict]:
    items = (entry.get("source") or {}).get("urls") or []
    if not items:
        print(f"  !! {entry['name']}: source.urls is empty (status {entry.get('status')}); "
              f"list the repository's file folder and fill config/corpus.yaml")
        return []
    records = []
    for item in items:
        url, name = url_target_name(item)
        record = download(url, folder / name)
        if name.lower().endswith(ARCHIVE_SUFFIXES):
            record["extraction"] = extract_archive(folder / name, folder)
        records.append(record)
    return records


def _proxy_tiledb_config() -> dict:
    """TileDB-SOMA's S3 client ignores http(s)_proxy; hand it the proxy explicitly.

    On Biowulf compute nodes `biowulf/run.sh` exports https_proxy=http://dtn03-e0:3128.
    """
    proxy = os.environ.get("https_proxy") or os.environ.get("HTTPS_PROXY") or os.environ.get("http_proxy")
    if not proxy:
        return {}
    u = urllib.parse.urlparse(proxy if "://" in proxy else f"http://{proxy}")
    cfg = {"vfs.s3.proxy_host": u.hostname, "vfs.s3.proxy_port": str(u.port or 3128),
           "vfs.s3.proxy_scheme": u.scheme or "http"}
    if u.username:
        cfg["vfs.s3.proxy_username"] = u.username
    if u.password:
        cfg["vfs.s3.proxy_password"] = u.password
    return cfg


DISCOVER_API = "https://api.cellxgene.cziscience.com/curation/v1"


def _discover_collection(collection: str) -> pd.DataFrame:
    """The datasets of a CELLxGENE Discover collection, with their public h5ad URLs.

    A plain HTTPS GET, so it honours the proxy like every other download here;
    this is the route that works on a Biowulf compute node.
    """
    with urllib.request.urlopen(f"{DISCOVER_API}/collections/{collection}", timeout=120) as r:
        payload = json.load(r)
    rows = []
    for ds in payload.get("datasets", []):
        h5ad = [a for a in ds.get("assets", []) if str(a.get("filetype", "")).upper() == "H5AD"]
        if not h5ad:
            continue
        rows.append({"dataset_id": ds.get("dataset_id") or ds.get("id"),
                     "dataset_title": ds.get("title"),
                     "dataset_total_cell_count": int(ds.get("cell_count") or 0),
                     "h5ad_url": h5ad[0].get("url"), "h5ad_bytes": h5ad[0].get("filesize"),
                     "collection_id": collection, "collection_name": payload.get("name"),
                     "dataset_version_id": ds.get("dataset_version_id")})
    return pd.DataFrame(rows)


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
    # Route 1 (default when the authors' h5ad is wanted): the Discover API +
    # a plain HTTPS download of the public h5ad. Works through the proxy.
    if collection and src.get("author_h5ad", True) and src.get("route", "discover") == "discover":
        try:
            table = _discover_collection(collection)
        except Exception as exc:                       # noqa: BLE001 - fall through to SOMA
            print(f"  !! Discover API failed ({exc!r}); falling back to the Census SOMA route")
            table = pd.DataFrame()
        if dataset_id and not table.empty:
            table = table[table["dataset_id"] == dataset_id]
        if not table.empty:
            table = table.sort_values("dataset_total_cell_count", ascending=False)
            chosen = table.iloc[0]
            print(f"  discover: {len(table)} dataset(s) in the collection; taking "
                  f"{chosen['dataset_id']} ({chosen['dataset_total_cell_count']:,} cells): {chosen.get('dataset_title')}")
            dest.parent.mkdir(parents=True, exist_ok=True)
            record = download(chosen["h5ad_url"], dest)
            table.to_csv(dest.parent / "census_datasets.csv", index=False)
            return [{"dataset_id": chosen["dataset_id"], "path": str(dest), "route": "discover",
                     "dataset_version_id": chosen.get("dataset_version_id"), "sha256": record["sha256"]}]
    # Route 2: the Census SOMA object on S3, with the proxy handed to TileDB.
    from tiledbsoma import SOMATileDBContext
    context = SOMATileDBContext(tiledb_config=_proxy_tiledb_config())
    with cellxgene_census.open_soma(census_version=version, context=context) as census:
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
        adata = cellxgene_census.get_anndata(census=cellxgene_census.open_soma(census_version=version, context=context),
                                             organism="Mus musculus" if entry.get("species") == "mouse" else "Homo sapiens",
                                             obs_value_filter=f"dataset_id == '{chosen['dataset_id']}' and is_primary_data == True",
                                             X_name="raw")
        adata.write_h5ad(dest)
    table.to_csv(dest.parent / "census_datasets.csv", index=False)
    return [{"dataset_id": chosen["dataset_id"], "path": str(dest), "census_version": version,
             "sha256": sha256_of(dest) if dest.exists() else None}]


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
        if kind == "local":
            print(f"  {'found' if path.exists() else 'MISSING'}: {path}")
            continue
        if kind == "url":
            records = _fetch_urls(entry, folder)
        elif kind == "cellxgene_collection":
            records = _fetch_census(entry, path if path.suffix == ".h5ad" else folder / f"{entry['name']}.h5ad")
        elif kind == "scbasecount":
            records = _fetch_scbasecount(cfg, entry, folder)
        else:
            raise ValueError(f"{entry['name']}: unknown source type {kind!r}")
        if records:
            folder.mkdir(parents=True, exist_ok=True)
            (folder / "download.json").write_text(json.dumps(records, indent=2, default=str), encoding="utf-8")
