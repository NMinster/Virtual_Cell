#!/usr/bin/env python
"""The ovary samples of scBaseCount (Arc Institute, SRA-wide, uniformly processed).

    python scripts/scbasecount_metadata.py --species mouse
    python scripts/scbasecount_metadata.py --species human --release 2026-01-12

Downloads the species' sample-metadata parquet from the public bucket over
anonymous HTTPS (the bucket answers plain GET requests; gsutil is not
needed), keeps the ovary experiments (UBERON:0000992 by ontology id, or
"ovar" / "follic" / "granulosa" / "oocyte" / "cumulus" in the tissue
text), and writes data/scbasecount/<species>_ovary_candidates.csv with
the columns the manifest needs (srx_accession, obs_count, lib_prep,
tech_10x, cell_prep, tissue, disease, perturbation, cell_line,
file_path). It prints the HTTPS URL of one sample's h5ad so a single
download can be tested before `ovary-fm fetch` runs a whole series.

The bucket layout and the column names are those of the 2026-01-12
release as the research memo recorded them; if a column is missing the
script says which and writes what it has.
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import pandas as pd

BUCKET = "arc-institute-virtual-cell-atlas"
HTTPS = f"https://storage.googleapis.com/{BUCKET}/"          # the metadata parquet lives under the bucket
GCS_HTTPS = "https://storage.googleapis.com/"                  # a gs://bucket/key path maps to this + bucket/key
ORGANISM = {"mouse": "Mus_musculus", "human": "Homo_sapiens", "macaque": "Macaca_mulatta", "pig": "Sus_scrofa",
            "sheep": "Ovis_aries"}
OVARY_TERMS = ("ovar", "follic", "granulosa", "oocyte", "cumulus", "corpus luteum", "theca")


def metadata_url(release: str, feature: str, organism: str) -> str:
    return f"{HTTPS}scbasecount/{release}/metadata/{feature}/{organism}/sample_metadata.parquet"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--species", default="mouse", choices=sorted(ORGANISM))
    parser.add_argument("--release", default="2026-01-12")
    parser.add_argument("--feature", default="GeneFull_Ex50pAS", help="the STARsolo feature set (introns included)")
    parser.add_argument("--out", default="data/scbasecount")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    organism = ORGANISM[args.species]
    url = metadata_url(args.release, args.feature, organism)
    local = out / f"{args.species}_sample_metadata.parquet"
    if local.exists() and not args.overwrite:
        print(f"= {local} already on disk")
    else:
        print(f"> {url}")
        with urllib.request.urlopen(url, timeout=300) as response, open(local, "wb") as handle:
            for chunk in iter(lambda: response.read(1 << 20), b""):
                handle.write(chunk)
        print(f"< {local}: {local.stat().st_size:,} bytes")
    meta = pd.read_parquet(local)
    print(f"  {len(meta):,} samples, columns: {list(meta.columns)}")

    tissue = meta["tissue"].astype(str).str.lower() if "tissue" in meta else pd.Series("", index=meta.index)
    by_id = meta["tissue_ontology_term_id"].astype(str).eq("UBERON:0000992") if "tissue_ontology_term_id" in meta \
        else pd.Series(False, index=meta.index)
    by_text = tissue.apply(lambda t: any(term in t for term in OVARY_TERMS))
    cand = meta[by_id | by_text].copy()
    cand["match"] = ["ontology" if i else "text" for i in by_id[cand.index]]
    keep = [c for c in ("srx_accession", "entrez_id", "obs_count", "lib_prep", "tech_10x", "cell_prep", "tissue",
                        "tissue_ontology_term_id", "disease", "disease_ontology_term_id", "perturbation", "cell_line",
                        "czi_collection_id", "czi_collection_name", "file_path", "match") if c in cand]
    cand = cand[keep].sort_values("obs_count", ascending=False) if "obs_count" in cand else cand[keep]
    path = out / f"{args.species}_ovary_candidates.csv"
    cand.to_csv(path, index=False)
    n_cells = int(cand["obs_count"].sum()) if "obs_count" in cand else -1
    print(f"  ovary candidates: {len(cand)} experiments, {n_cells:,} cells -> {path}")
    if "cell_prep" in cand:
        print(cand.groupby("cell_prep")["obs_count"].agg(["size", "sum"]).to_string())
    if "perturbation" in cand:
        pert = cand[cand["perturbation"].astype(str).str.lower().ne("none") & cand["perturbation"].notna()]
        print(f"  with a perturbation annotation: {len(pert)} experiments")
        print(pert[["srx_accession", "obs_count", "perturbation"]].head(20).to_string(index=False))
    if len(cand) and "file_path" in cand:
        gs = str(cand["file_path"].iloc[0])
        print(f"\n  one h5ad to test a download:\n    {GCS_HTTPS + gs[5:] if gs.startswith('gs://') else gs}")


if __name__ == "__main__":
    main()
