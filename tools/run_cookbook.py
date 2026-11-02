"""Reduce the M15 LCB HR-R dataset of the MEGARA cookbook

The dataset is downloaded and extracted in the current directory, and
the calibration tree is created in the directory of M15, as
megaradrp-inittree does, with the dataset of calibrations declared by
megaradrp (the last one, or the one given with --dataset). The products
are found in the registry of the reductions, there is no need to copy
them to the calibration tree.

Usage: python run_cookbook.py [--dataset NAME] [tarball]

With a tarball, the local file is used instead of the download.
"""

import argparse
import logging
import os
import pathlib
import tarfile

from megaradrp.loader import load_drp
from megaradrp.tools.init_tree import ensure_numina_cfg, init_tree, install_calibration_data
from numina.testing.testcache import download_cache
from numina.user.inittree import find_dataset
from numina.user.session import Session

DATASET = "megara-cookbook-v4"
URL = f"https://guaix.fis.ucm.es/data/megaradrp/{DATASET}.tar.gz"

OBSERVATIONS = [
    ("0_Bias.yaml", "0_Bias"),
    ("1_TraceMap.yaml", "1_TraceMap_HR-R"),
    ("3_WaveCalib.yaml", "3_WaveCalib_HR-R"),
    ("4_FiberFlat.yaml", "4_FiberFlat_HR-R"),
    ("5_TwilightFlat.yaml", "5_TwilightFlat_HR-R"),
    ("6_LcbAdquisition.yaml", "6_LcbImage_HR-R"),
    ("7_StandardStar.yaml", "7_StandardStar_HR-R"),
    ("8_LcbImage.yaml", "8_LcbImage_HR-R_M15"),
]


def main(args=None):

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("tarball", nargs="?", help="local tarball of the cookbook, instead of the download")
    parser.add_argument("--dataset", help="dataset of calibrations of megaradrp, by default the last one")
    args = parser.parse_args(args)

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    if args.tarball:
        tarball = pathlib.Path(args.tarball)
    else:
        tarball = pathlib.Path(download_cache(URL).name)

    with tarfile.open(tarball, mode="r:gz") as tar:
        tar.extractall(filter="data")

    if not args.tarball:
        os.remove(tarball)

    basedir = pathlib.Path(DATASET, "M15_LCB_HR-R").resolve()
    # .numina.cfg and the calibration tree are in the directory of M15
    os.chdir(basedir)
    dataset = find_dataset(load_drp(), args.dataset)
    init_tree(base_path="calibrations", logger=logger)
    install_calibration_data(
        dataset["url"], "calibrations", overwrite=False, logger=logger, known_hash=dataset["known_hash"]
    )
    ensure_numina_cfg(logger=logger, base_path="calibrations", dataset=dataset)

    # Registry of the reductions, the products are found there in later reductions
    session = Session(basedir=basedir, db="numina-db.json")
    session.add_observations(*(obsfile for obsfile, _ in OBSERVATIONS))

    for _, obsid in OBSERVATIONS:
        session.run(obsid)


if __name__ == "__main__":

    main()
