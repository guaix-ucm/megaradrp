#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
"""Generate the calibration data tree structure.

The calibration data is a dataset declared by megaradrp in drp.yaml,
stored in a ZIP file available at Zenodo. The ZIP file is downloaded to
the cache of numina, its hash is checked, and its files are installed in
the directory MEGARA under the base path, as ``numina inittree MEGARA``
does (see :mod:`numina.user.inittree`): the directory at the top of the
ZIP file is removed from the paths, the hidden files are skipped, the
files that already exist with the same contents are left as they are,
and if some file exists with other contents, nothing is installed
unless --overwrite is given.

The tree of directories for the products of the reductions is created
too, and ``calibsdir`` is set to the base path in .numina.cfg.
"""

import argparse
import configparser
from datetime import datetime
import getpass
import logging
from pathlib import Path
import sys

from rich.logging import RichHandler
from rich_argparse import RichHelpFormatter

from megaradrp import __version__
from numina.user.console import NuminaConsole


def init_tree(base_path=None, dry_run=False, logger=None):
    """Create the calibration data tree structure.

    Parameters
    ----------
    base_path : str or Path, optional
        Base path where the calibration data tree will be created.
    dry_run : bool, optional
        If True, do not create directories, only print what would be done.
    logger : logging.Logger, optional
        Logger to use for logging messages. If None, a default logger will be used.
    """

    if base_path is None:
        raise ValueError("Base path must be provided.")

    if logger is None:
        logger = logging.getLogger(__name__)

    if not isinstance(base_path, (str, Path)):
        raise TypeError("Base path must be a string or a Path object.")
    if isinstance(base_path, str):
        base_path = Path(base_path)

    calibration_id_list = [
        "ca3558e3-e50d-4bbc-86bd-da50a0998a48",
        "9f3fecbf-376c-47ae-a188-313b7d829104",
    ]

    calibration_directory_list = [
        "LinesCatalog",
        "MasterBPM",
        "MasterBias",
        "MasterFiberFlat",
        "MasterSensitivity",
        "MasterSlitFlat",
        "MasterTwilightFlat",
        "ModelMap",
        "ReferenceExtinctionTable",
        "TraceMap",
        "WavelengthCalibration",
    ]

    ifu_mode_list = [
        "LCB",
        "MOS",
    ]

    arc_lamp_list = [
        "ThAr",
        "ThNe",
    ]

    # LineCatalogs are not available for all VPH gratings
    vph_grating_dict = {
        "HR-I": {"ThAr": True, "ThNe": True},
        "HR-R": {"ThAr": True, "ThNe": True},
        "LR-B": {"ThAr": True, "ThNe": False},
        "LR-I": {"ThAr": True, "ThNe": True},
        "LR-R": {"ThAr": True, "ThNe": True},
        "LR-U": {"ThAr": True, "ThNe": False},
        "LR-V": {"ThAr": True, "ThNe": False},
        "LR-Z": {"ThAr": True, "ThNe": True},
        "MR-B": {"ThAr": True, "ThNe": False},
        "MR-G": {"ThAr": True, "ThNe": False},
        "MR-I": {"ThAr": False, "ThNe": True},
        "MR-R": {"ThAr": True, "ThNe": True},
        "MR-RI": {"ThAr": True, "ThNe": True},
        "MR-U": {"ThAr": True, "ThNe": False},
        "MR-UB": {"ThAr": True, "ThNe": False},
        "MR-V": {"ThAr": True, "ThNe": True},
        "MR-VR": {"ThAr": False, "ThNe": True},
        "MR-Z": {"ThAr": False, "ThNe": True},
    }

    # for historical reasons, the directory tree includes "MEGARA"
    # as the top-level directory under the base path
    logger.info(f"Creating calibration data tree structure under base path: {base_path}")
    for calibration_id in calibration_id_list:
        for calibration_directory in calibration_directory_list:
            if calibration_directory in ["MasterBPM", "MasterBias", "MasterSlitFlat", "ReferenceExtinctionTable"]:
                path = base_path / "MEGARA" / calibration_id / calibration_directory
                logger.debug(f"Creating directory: {path}")
                if not dry_run:
                    path.mkdir(parents=True, exist_ok=True)
            elif calibration_directory == "LinesCatalog":
                for arc_lamp in arc_lamp_list:
                    for vph_grating in vph_grating_dict.keys():
                        if vph_grating_dict[vph_grating][arc_lamp]:
                            path = (
                                base_path / "MEGARA" / calibration_id / calibration_directory / arc_lamp / vph_grating
                            )
                            logger.debug(f"Creating directory: {path}")
                            if not dry_run:
                                path.mkdir(parents=True, exist_ok=True)
            else:
                for ifu_mode in ifu_mode_list:
                    for vph_grating in vph_grating_dict.keys():
                        path = base_path / "MEGARA" / calibration_id / calibration_directory / ifu_mode / vph_grating
                        logger.debug(f"Creating directory: {path}")
                        if not dry_run:
                            path.mkdir(parents=True, exist_ok=True)


def declared_dataset(url=None):
    """The dataset of calibrations declared by megaradrp in drp.yaml

    If `url` is None, the last dataset declared. Otherwise, the dataset
    with that URL, or None if megaradrp does not declare it.
    """
    from megaradrp.loader import load_drp

    datasets = load_drp().calibrations
    if url is None:
        return datasets[-1]
    for dataset in datasets:
        if dataset["url"] == url:
            return dataset
    return None


def install_calibration_data(url, base_path, overwrite, dry_run=False, logger=None, known_hash=None):
    """Install the calibration data from the given URL.

    The ZIP file is downloaded to the cache of numina, and it is not
    downloaded again. Its files are installed in the directory MEGARA
    under `base_path`, see :func:`numina.user.inittree.plan_install`.

    Parameters
    ----------
    url : str
        URL of the calibration data archive to download and extract.
    base_path : str or Path
        Base path where the calibration data will be installed.
    overwrite : bool
        If True, replace the existing files that have other contents.
        If False, nothing is installed if any such file is found.
    dry_run : bool, optional
        If True, download the ZIP file and print what would be done,
        without installing any file.
    logger : logging.Logger, optional
        Logger to use for logging messages. If None, a default logger will be used.
    known_hash : str, optional
        Hash of the ZIP file, as "md5:value". If None, the hash of the
        dataset declared by megaradrp with this URL is used, and if there
        is no such dataset, the hash is not checked.

    Returns
    -------
    numina.user.inittree.InstallPlan
        The files that are new, unchanged and with other contents.

    Raises
    ------
    SystemExit
        If some file exists with other contents and `overwrite` is False.
    """
    from numina.user.inittree import extract, plan_install
    from numina.util.download import download_file

    if logger is None:
        logger = logging.getLogger(__name__)

    if not isinstance(base_path, (str, Path)):
        raise TypeError("Base path must be a string or a Path object.")
    if isinstance(base_path, str):
        base_path = Path(base_path)

    if known_hash is None:
        dataset = declared_dataset(url)
        if dataset is None:
            logger.warning(f"The URL is not a dataset declared by megaradrp, its hash is not checked: {url}")
        else:
            known_hash = dataset["known_hash"]

    logger.info(f"Downloading calibration data from: {url}")
    zip_path = download_file(url, known_hash=known_hash)

    # for historical reasons, the files are installed under "MEGARA"
    destination = base_path / "MEGARA"
    plan = plan_install(zip_path, destination)
    logger.info(
        f"{len(plan.new)} new files, {len(plan.unchanged)} unchanged, "
        f"{len(plan.conflicts)} with other contents, in {destination}"
    )
    for _, target in plan.conflicts:
        logger.warning(f"File with other contents: {target}")
    if plan.conflicts and not overwrite:
        logger.warning("Nothing was installed. Use the --overwrite option to replace these files.")
        raise SystemExit(1)

    if dry_run:
        logger.info(f"Dry run: {len(plan.new) + len(plan.conflicts)} files would be installed in {destination}")
        return plan

    extract(zip_path, plan.new + plan.conflicts)
    logger.info(f"Installed {len(plan.new) + len(plan.conflicts)} files in {destination}")
    return plan


def ensure_numina_cfg(dry_run=False, logger=None, base_path="calibrations", dataset=None):
    """Ensure that the .numina.cfg file exists and contains the necessary configuration.

    The section ``[tool.megara.inittree]`` records the version of
    megaradrp, the user and the date, and ``calibsdir`` in ``[tool.run]``
    is set to `base_path`. If `dataset` is given, the dataset installed
    is recorded in ``[tool.inittree.MEGARA]``, as numina inittree does
    (see :func:`numina.user.inittree.update_config`), so that numina run
    warns if it is not the last dataset declared by megaradrp.

    Parameters
    ----------
    dry_run : bool, optional
        If True, do not actually write to the file. Only print what would be done.
    logger : logging.Logger, optional
        Logger to use for logging messages. If None, a default logger will be used.
    base_path : str or Path, optional
        Base path of the calibration data, the value of ``calibsdir``.
    dataset : dict, optional
        The dataset installed, with its name, url and known_hash.
    """
    from numina.user.inittree import update_config

    from megaradrp.loader import load_drp

    if logger is None:
        logger = logging.getLogger(__name__)

    # read the .numina.cfg file if it exists, or create a new one if it does not exist
    cfg_path = Path(".numina.cfg")
    config = configparser.ConfigParser()

    if dry_run:
        if not cfg_path.exists():
            logger.info(f"Creating new configuration file: {cfg_path.resolve()}")
        logger.info(f"calibsdir would be set to {base_path} in {cfg_path.resolve()}")
        return

    if cfg_path.exists():
        config.read(cfg_path)
    else:
        logger.info(f"Creating new configuration file: {cfg_path.resolve()}")

    # insert section with information
    if "tool.megara.inittree" not in config:
        config["tool.megara.inittree"] = {}
    config["tool.megara.inittree"]["description"] = "Auto-generated configuration file for calibration data."
    config["tool.megara.inittree"]["megaradrp-version"] = __version__
    config["tool.megara.inittree"]["username"] = getpass.getuser()
    config["tool.megara.inittree"]["creation-date"] = datetime.now().isoformat(timespec="seconds")

    # ensure the [tool.run] section exists
    if "tool.run" not in config:
        config["tool.run"] = {}

    if dataset is None:
        # set calibrations directory
        config["tool.run"]["calibsdir"] = str(base_path)

    # write the updated configuration back to the .numina.cfg file
    with cfg_path.open("w") as f:
        config.write(f)

    if dataset is not None:
        # set calibrations directory and record the dataset
        update_config(cfg_path, base_path, load_drp(), dataset)

    logger.info(f"Configuration file updated: {cfg_path.resolve()}")


def main(args=None):
    """Main function to create the calibration data tree structure."""
    # parse command-line options
    parser = argparse.ArgumentParser(
        description="Create the calibration data tree structure.",
        formatter_class=RichHelpFormatter,
    )
    # positional parameters
    parser.add_argument(
        "--base-path",
        help="Base path where the calibration data tree will be created.",
        type=str,
        default="./calibrations",
    )
    parser.add_argument(
        "--zip-url",
        help="URL of the calibration data ZIP file to download (default: the last dataset declared by megaradrp).",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="If set, overwrite existing files when moving extracted files."
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="If set, do not create directories, only print what would be done."
    )
    parser.add_argument(
        "--log-level", help="Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).", type=str, default="INFO"
    )
    parser.add_argument("--record", action="store_true", help="Record terminal output.")
    parser.add_argument("--echo", help="Display full command line", action="store_true")
    args = parser.parse_args(args)

    # initialize the console and logger
    console = NuminaConsole(record=args.record)
    logger = logging.getLogger(__name__)

    if args.echo:
        console.print(f"[bright red]Executing: {' '.join(sys.argv)}[/bright red]\n", end="")

    # configure logging
    if args.log_level in ["DEBUG", "WARNING", "ERROR", "CRITICAL"]:
        format_log = "%(name)s %(levelname)s: %(message)s"
        handlers = [RichHandler(console=console, show_time=False, markup=True)]
    else:
        format_log = "%(message)s"
        handlers = [RichHandler(console=console, show_time=False, markup=True, show_path=False, show_level=False)]
    logging.basicConfig(level=args.log_level, format=format_log, handlers=handlers)
    logging.getLogger("matplotlib").setLevel(logging.ERROR)  # supress matplotlib debug logs

    # welcome message
    console.rule("[bold magenta]Megara DRP - Calibration Data Tree Initialization[/bold magenta]")

    # create the directory tree structure
    init_tree(base_path=args.base_path, dry_run=args.dry_run, logger=logger)

    # download and install the calibration data, by default the last
    # dataset declared by megaradrp in drp.yaml
    if args.zip_url is None:
        dataset = declared_dataset()
    else:
        dataset = declared_dataset(args.zip_url)
        if dataset is None:
            dataset = {"name": Path(args.zip_url).stem, "url": args.zip_url, "known_hash": ""}
    install_calibration_data(
        url=dataset["url"],
        base_path=args.base_path,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        logger=logger,
        known_hash=dataset["known_hash"] or None,
    )

    # ensure the .numina.cfg file exists and contains the necessary configuration
    ensure_numina_cfg(dry_run=args.dry_run, logger=logger, base_path=args.base_path, dataset=dataset)

    # goodbye message
    logger.info(f"Calibration data tree structure created under base path: {args.base_path}")
    console.rule("[bold magenta]Calibration Data Tree Initialization Completed[/bold magenta]")
    if args.dry_run:
        logger.info("This was a dry run, no actual changes were made!")


if __name__ == "__main__":
    main()
