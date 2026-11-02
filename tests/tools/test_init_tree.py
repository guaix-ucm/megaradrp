"""megaradrp-inittree installs the dataset of calibrations of drp.yaml as numina inittree does"""

import configparser
import logging
import zipfile

import pytest

import megaradrp.tools.init_tree as init_tree
from megaradrp.loader import load_drp

URL_2026_2 = "https://zenodo.org/records/18623771/files/guaix-ucm/megaradrp-calibrations-2026.2.zip"

FILES = {
    "guaix-ucm-megaradrp-calibrations-abc1234/README.md": b"readme",
    "guaix-ucm-megaradrp-calibrations-abc1234/.zenodo.json": b"{}",
    "guaix-ucm-megaradrp-calibrations-abc1234/LinesCatalog/LR-U_ThAr.ecsv": b"# %ECSV 1.0\n",
    "guaix-ucm-megaradrp-calibrations-abc1234/MasterBPM/master_bpm_9f3fecbf.fits.gz": b"bpm",
}


def test_declared_dataset():
    datasets = load_drp().calibrations
    assert init_tree.declared_dataset() == datasets[-1]
    assert init_tree.declared_dataset(URL_2026_2)["known_hash"] == "md5:3d98803089d4822ee00ed03a4dbba229"
    assert init_tree.declared_dataset("https://example.org/other.zip") is None


@pytest.fixture
def download(monkeypatch, tmp_path):
    """numina download_file, recording its arguments, returning a ZIP file like those of Zenodo"""
    zip_path = tmp_path / "dataset.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        for name, content in FILES.items():
            zf.writestr(name, content)
    calls = []

    def fake_download(url, known_hash=None):
        calls.append({"url": url, "known_hash": known_hash})
        return str(zip_path)

    monkeypatch.setattr("numina.util.download.download_file", fake_download)
    return calls


@pytest.mark.parametrize("dataset", [-1, 0])
def test_install_with_declared_hash(download, tmp_path, dataset):
    declared = load_drp().calibrations[dataset]
    init_tree.install_calibration_data(declared["url"], tmp_path / "calib", overwrite=False)
    assert download == [{"url": declared["url"], "known_hash": declared["known_hash"]}]


def test_install_other_url(download, tmp_path, caplog):
    """An URL that megaradrp does not declare is downloaded without checking its hash"""
    url = "https://example.org/other.zip"
    with caplog.at_level(logging.WARNING):
        init_tree.install_calibration_data(url, tmp_path / "calib", overwrite=False)
    assert download[0]["known_hash"] is None
    assert "its hash is not checked" in caplog.text


def test_install_given_hash(download, tmp_path):
    init_tree.install_calibration_data(URL_2026_2, tmp_path / "calib", overwrite=False, known_hash="md5:0")
    assert download[0]["known_hash"] == "md5:0"


def test_install_files(download, tmp_path):
    """The root directory of the ZIP file is removed, the hidden files are skipped"""
    base = tmp_path / "calib"
    plan = init_tree.install_calibration_data(URL_2026_2, base, overwrite=False)
    assert len(plan.new) == 3
    assert (base / "MEGARA/LinesCatalog/LR-U_ThAr.ecsv").read_bytes() == b"# %ECSV 1.0\n"
    assert (base / "MEGARA/MasterBPM/master_bpm_9f3fecbf.fits.gz").read_bytes() == b"bpm"
    assert not (base / "MEGARA/.zenodo.json").exists()
    assert sorted(path.name for path in base.iterdir()) == ["MEGARA"]


def test_install_again(download, tmp_path):
    """The files with the same contents are left as they are"""
    base = tmp_path / "calib"
    init_tree.install_calibration_data(URL_2026_2, base, overwrite=False)
    plan = init_tree.install_calibration_data(URL_2026_2, base, overwrite=False)
    assert (len(plan.new), len(plan.unchanged), len(plan.conflicts)) == (0, 3, 0)


def test_install_conflicts(download, tmp_path):
    """With a file with other contents, nothing is installed, unless overwrite"""
    base = tmp_path / "calib"
    readme = base / "MEGARA/README.md"
    readme.parent.mkdir(parents=True)
    readme.write_bytes(b"other")
    with pytest.raises(SystemExit):
        init_tree.install_calibration_data(URL_2026_2, base, overwrite=False)
    assert sorted(path.name for path in (base / "MEGARA").iterdir()) == ["README.md"]
    assert readme.read_bytes() == b"other"

    init_tree.install_calibration_data(URL_2026_2, base, overwrite=True)
    assert readme.read_bytes() == b"readme"
    assert (base / "MEGARA/LinesCatalog/LR-U_ThAr.ecsv").is_file()


def test_install_dry_run(download, tmp_path):
    base = tmp_path / "calib"
    plan = init_tree.install_calibration_data(URL_2026_2, base, overwrite=False, dry_run=True)
    assert len(plan.new) == 3
    assert not base.exists()


def test_numina_cfg(tmp_path, monkeypatch):
    """calibsdir is the base path, and the dataset is recorded as numina inittree does"""
    monkeypatch.chdir(tmp_path)
    dataset = init_tree.declared_dataset()
    init_tree.ensure_numina_cfg(base_path="calib-new", dataset=dataset)
    config = configparser.ConfigParser()
    config.read(tmp_path / ".numina.cfg")
    assert config["tool.run"]["calibsdir"] == "calib-new"
    assert config["tool.inittree.MEGARA"]["dataset"] == dataset["name"]
    assert config["tool.inittree.MEGARA"]["calibsdir"] == str(tmp_path / "calib-new")
    assert "megaradrp-version" in config["tool.megara.inittree"]


def test_numina_cfg_without_dataset(tmp_path, monkeypatch):
    """The call of previous versions, used in the notebooks of the cookbook"""
    monkeypatch.chdir(tmp_path)
    init_tree.ensure_numina_cfg()
    config = configparser.ConfigParser()
    config.read(tmp_path / ".numina.cfg")
    assert config["tool.run"]["calibsdir"] == "calibrations"
    assert not config.has_section("tool.inittree.MEGARA")


def test_main(download, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    init_tree.main(["--base-path", "calib-new"])
    assert download[0]["url"] == init_tree.declared_dataset()["url"]
    assert (tmp_path / "calib-new/MEGARA/LinesCatalog/LR-U_ThAr.ecsv").is_file()
    # the tree of directories for the products
    assert (tmp_path / "calib-new/MEGARA/9f3fecbf-376c-47ae-a188-313b7d829104/MasterBias").is_dir()
    config = configparser.ConfigParser()
    config.read(tmp_path / ".numina.cfg")
    assert config["tool.run"]["calibsdir"] == "calib-new"
    assert config["tool.inittree.MEGARA"]["dataset"] == init_tree.declared_dataset()["name"]
