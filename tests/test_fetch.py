"""Units for ovary_fm.fetch: URL naming, archive extraction, TileDB proxy config."""
import os
import zipfile

import pytest

from ovary_fm import fetch


def test_url_target_name_from_path_segment():
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE294nnn/GSE294534/suppl/GSE294534_10x_update.rds"
    assert fetch.url_target_name(url) == (url, "GSE294534_10x_update.rds")


def test_url_target_name_from_dict_name():
    item = {"url": "https://osf.io/download/dsf4r/", "name": "musCyclingOvaryData.zip"}
    assert fetch.url_target_name(item) == (item["url"], "musCyclingOvaryData.zip")


def test_url_target_name_refuses_nameless_opaque_url():
    with pytest.raises(SystemExit):
        fetch.url_target_name("https://osf.io/download/dsf4r/")


def test_extract_zip_once(tmp_path):
    archive = tmp_path / "data.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("sub/matrix.mtx", "%%MatrixMarket\n1 1 1\n1 1 3\n")
        z.writestr("sub/barcodes.tsv", "AAAC\n")
    first = fetch.extract_archive(archive, tmp_path)
    assert first["extracted"] and not first["skipped"] and first["members"] == 2
    assert (tmp_path / "sub" / "matrix.mtx").read_text().startswith("%%MatrixMarket")
    second = fetch.extract_archive(archive, tmp_path)
    assert second["skipped"]


def test_extract_zip_rejects_escape(tmp_path):
    archive = tmp_path / "evil.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("../outside.txt", "x")
    with pytest.raises(SystemExit):
        fetch.extract_archive(archive, tmp_path)


def test_proxy_tiledb_config(monkeypatch):
    monkeypatch.delenv("HTTPS_PROXY", raising=False)
    monkeypatch.delenv("http_proxy", raising=False)
    monkeypatch.setenv("https_proxy", "http://dtn03-e0:3128")
    cfg = fetch._proxy_tiledb_config()
    assert cfg == {"vfs.s3.proxy_host": "dtn03-e0", "vfs.s3.proxy_port": "3128", "vfs.s3.proxy_scheme": "http"}
    monkeypatch.delenv("https_proxy")
    assert fetch._proxy_tiledb_config() == {}
