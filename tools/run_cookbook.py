import os
import pathlib
import logging
import tarfile

from numina.user.session import Session
from numina.testing.testcache import download_cache


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

    # Registry of the reductions, the products are found there in later reductions
    session = Session(basedir=basedir, control="control.yaml", db="numina-db.json")

    session.add_observations(
        "0_bias.yaml",
        "2_M15_modelmap.yaml",
        "4_M15_fiberflat.yaml",
        "6_M15_Lcbadquisition.yaml",
        "8_M15_reduce_LCB.yaml",
        "1_M15_tracemap.yaml",
        "3_M15_wavecalib.yaml",
        "5_M15_twilight.yaml",
        "7_M15_Standardstar.yaml",
    )

    for obsid in ["0_bias", "1_HR-R", "3_HR-R", "4_HR-R", "5_HR-R", "6_HR-R", "7_HR-R", "8_HR-R"]:
        session.run(obsid)


if __name__ == "__main__":

    main()
