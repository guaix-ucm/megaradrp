import os
import pathlib
import logging
import tarfile

from numina.user.cli import base_config
from numina.user.helpers import create_datamanager, load_observations
from numina.user.baserun import run_reduce
from numina.tests.testcache import download_cache


def main():

    logging.basicConfig(level=logging.DEBUG)

    basedir = pathlib.Path().resolve()

    tarball = "MEGARA-cookbook-M15_LCB_HR-R-v2.tar.gz"
    url = "https://guaix.fis.ucm.es/~spr/megara_test/{}".format(tarball)

    downloaded = download_cache(url)

    # Uncompress
    with tarfile.open(downloaded.name, mode="r:gz") as tar:
        tar.extractall()

    os.remove(downloaded.name)

    config = base_config()
    config["tool.run"]["basedir"] = str(basedir)
    # Registry of the reductions, the products are found there in later reductions
    config["tool.db"]["file"] = "numina-db.json"

    dm = create_datamanager(config, str(basedir / "control.yaml"))

    obsresults = [
        "0_bias.yaml",
        "2_M15_modelmap.yaml",
        "4_M15_fiberflat.yaml",
        "6_M15_Lcbadquisition.yaml",
        "8_M15_reduce_LCB.yaml",
        "1_M15_tracemap.yaml",
        "3_M15_wavecalib.yaml",
        "5_M15_twilight.yaml",
        "7_M15_Standardstar.yaml",
    ]

    sessions, loaded_obs = load_observations([str(basedir / name) for name in obsresults])
    dm.backend.add_obs(loaded_obs)

    for obsid in ["0_bias", "1_HR-R", "3_HR-R", "4_HR-R", "5_HR-R", "6_HR-R", "7_HR-R", "8_HR-R"]:
        run_reduce(dm, obsid)


if __name__ == "__main__":

    main()
