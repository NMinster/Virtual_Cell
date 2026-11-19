"""The chunked csv reader, R-name unmangling, barcode regex, count floor, SCP metadata join."""
import numpy as np
import pandas as pd

from ovary_fm import io


def _write_csv(folder, n_genes=20, barcodes=None):
    barcodes = barcodes or [f"mary_1_1.1-1A_{i:06d}" for i in range(6)] + [f"mary_2_4.1-4D_{i:06d}" for i in range(6)]
    rng = np.random.default_rng(0)
    counts = rng.poisson(3.0, size=(n_genes, len(barcodes))).astype(float)
    counts[:, 0] = 0                                  # an empty droplet
    genes = ["X0610006L08Rik", "mt.Nd6"] + [f"Gene{i}" for i in range(n_genes - 2)]
    pd.DataFrame(counts, index=genes, columns=barcodes).to_csv(folder / "musMaryData.csv")
    return barcodes, counts


def test_chunked_csv_reader_matches_dense(tmp_path):
    barcodes, counts = _write_csv(tmp_path)
    adata = io.read_counts(tmp_path, {"csv_chunk_rows": 7})
    assert adata.shape == (len(barcodes), counts.shape[0])
    np.testing.assert_array_equal(adata.X.toarray(), counts.T)
    assert list(adata.var_names[:2]) == ["X0610006L08Rik", "mt.Nd6"]      # untouched without the option


def test_unmangle_barcode_regex_and_floor(tmp_path):
    barcodes, counts = _write_csv(tmp_path)
    adata = io.read_counts(tmp_path, {"unmangle_r_names": True, "csv_min_counts": 1,
                                      "barcode_regex": r"^mary_(?P<library>\d+)_(?P<sample>[^_]+)_\d+$"})
    assert list(adata.var_names[:2]) == ["0610006L08Rik", "mt-Nd6"]
    assert adata.n_obs == len(barcodes) - 1                               # the empty droplet is gone
    assert set(adata.obs["sample"]) == {"1.1-1A", "4.1-4D"}
    assert set(adata.obs["library"]) == {"1", "2"}


def test_scp_metadata_join_skips_type_row(tmp_path):
    barcodes, _ = _write_csv(tmp_path)
    meta = pd.DataFrame({"NAME": barcodes[:4], "stage": ["diestrus"] * 2 + ["estrus"] * 2,
                         "Level1": ["Granulosa", "Theca", "Stroma", "Stroma"]})
    lines = ["\t".join(meta.columns), "TYPE\tgroup\tgroup"] + ["\t".join(map(str, r)) for r in meta.to_numpy()]
    (tmp_path / "scp1914_metadata.txt").write_text("\n".join(lines) + "\n")
    adata = io.read_counts(tmp_path)
    assert adata.obs.loc[barcodes[0], "stage"] == "diestrus"
    assert adata.obs.loc[barcodes[3], "Level1"] == "Stroma"
    assert pd.isna(adata.obs.loc[barcodes[5], "stage"])
