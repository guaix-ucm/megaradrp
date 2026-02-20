#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
"""Generate the calibration data tree structure.

The calibration data is stored in ZIP file available at Zenodo, which is
downloaded and extracted to the base path. The ZIP file contains a directory
named "guaix-ucm-megaradrp-calibrations-40fb7c (or something similar, depending
on the version), which is stripped from the paths of the extracted files when
moving them to the base calibration path. The ZIP file is then removed,
and the now empty extracted directory is also removed.

The code detects the name of the extracted directory by looking
for the first file in the extracted directory and extracting the
directory name from it.
"""

import argparse
import logging
from pathlib import Path
import pooch
import shutil
import sys
import tempfile
import zipfile

from rich.logging import RichHandler
from rich_argparse import RichHelpFormatter

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


def move_files_strip_fixed(files, base_path, to_strip=None, overwrite=False, dry_run=False, logger=None):
    """Move files to the base path, stripping a fixed part of the path.

    Parameters
    ----------
    files : list of str or Path
        List of file paths to move.
    base_path : str or Path
        Base path where the files will be moved to.
    to_strip : str, optional
        Fixed part of the path to strip from the file paths.
        If None, no stripping is done.
    overwrite : bool, optional
        If True, overwrite existing files when moving.
        If False, abort if any existing files are found.
    dry_run : bool, optional
        If True, do not actually move files, only print what would be done.
    logger : logging.Logger, optional
        Logger to use for logging messages. If None, a default logger will be used.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if to_strip is None:
        raise ValueError("to_strip must be provided.")

    logger.info(f"Moving files to base path: '{base_path}'")
    logger.info(f"Stripping fixed part of the path: '{to_strip}'")

    for f in files:
        if dry_run:
            if f[-1] == "/":
                # directory
                continue
            if "calibration_data.zip" in f:
                # ignore the ZIP file itself
                continue
            if ".gitattributes" in f:
                # ignore .gitattributes file, which is not part of the calibration data
                continue
            f = Path(f)
        else:
            f = Path(f)
            # ignore directories, we only want to move files
            if f.is_dir():
                continue
            # ignore the ZIP file itself
            if f.name == "calibration_data.zip":
                continue
            # ignore .gitattributes file, which is not part of the calibration data
            if f.name == ".gitattributes":
                continue
        # use as_posix to ensure we have a consistent string representation of the path
        rel_text = f.as_posix()
        marker = to_strip.lstrip("./")  # normaliza el comienzo
        # search for the marker in the path
        pos = rel_text.find(marker)
        if pos == -1:
            # the file is not under the uncompressed directory, we can ignore it
            continue
        # replace only the first occurrence of the marker in the path
        if dry_run:
            dest = base_path.as_posix() + rel_text.replace(marker, "", 1)
        else:
            dest = rel_text.replace(marker, "", 1)
        # insert "MEGARA" as the top-level directory under the base path
        dest = dest.replace(base_path.as_posix(), f"{base_path.as_posix()}/MEGARA", 1)
        if not dry_run:
            # ensure the destination directory exists
            dest = Path(dest)
            dest.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Moving: '{f}'  ->  '{dest}'")
        # if destination file already exists, do not overwrite it
        if not dry_run:
            if dest.exists():
                if not overwrite:
                    logger.warning(f"Destination file already exists, aborting: {dest}")
                    logger.warning("Use the --overwrite option to overwrite existing files.")
                    raise SystemExit(1)
                else:
                    logger.warning(f"Destination file already exists, overwriting: {dest}")
            try:
                f.replace(dest)
            except OSError:
                shutil.copy2(f, dest)
                f.unlink()


def install_calibration_data(url, base_path, overwrite, dry_run=False, logger=None):
    """Install the calibration data from the given URL.

    Parameters
    ----------
    url : str
        URL of the calibration data archive to download and extract.
    base_path : str or Path
        Base path where the calibration data will be installed.
    overwrite : bool
        If True, overwrite existing files when moving extracted files.
        If False, abort if any existing files are found.
    dry_run : bool, optional
        If True, donwload the ZIP file in a temporary location to
        extract the list of files that would be extracted, but do not actually
        create directories or move files. This is useful to check what would be done
        without making any changes to the file system.
    logger : logging.Logger, optional
        Logger to use for logging messages. If None, a default logger will be used.
    """
    if logger is None:
        logger = logging.getLogger(__name__)

    if not isinstance(base_path, (str, Path)):
        raise TypeError("Base path must be a string or a Path object.")
    if isinstance(base_path, str):
        base_path = Path(base_path)

    local_zip_filename = "calibration_data.zip"
    logger.info(f"Downloading and installing calibration data from: {url}")
    if dry_run:
        logger.debug("Dry run mode: downloading ZIP file to temporary location.")
        # download the ZIP file but do not extract it, since we only want to get 
        # the list of files that would be extracted, without actually creating 
        # directories or moving files. Note that the ZIP file is downloaded to
        # a temporary location, and it is removed after leaving the with block, 
        # so it does not affect the file system.
        with tempfile.TemporaryDirectory() as tmpdir:
            zip_path = pooch.retrieve(
                url=url,
                known_hash="md5:3d98803089d4822ee00ed03a4dbba229",
                fname=local_zip_filename,
                path=tmpdir,
                progressbar=True,
            )
            # get the list of files that would be extracted from the ZIP file
            files = zipfile.ZipFile(zip_path, "r").namelist()
    else:
        # download and extract the ZIP file to the base path
        files = pooch.retrieve(
            url=url,
            known_hash="md5:3d98803089d4822ee00ed03a4dbba229",
            fname=local_zip_filename,
            path=base_path,
            progressbar=True,
            processor=pooch.Unzip(extract_dir="."),
        )

    # the extracted files appear under a directory named
    # "guaix-ucm-megaradrp-calibrations-40fb7c2"
    # (or something similar, depending on the version),
    # so we need to strip that part of the path when moving the files
    # to the base path. First, we need to find the actual name of the
    # extracted directory, which we can do by looking for the first file in the list
    # that contains the expected marker in its path.
    dir_to_strip = None
    for f in files:
        p = Path(f).resolve()
        for part in p.parts:
            if part.startswith("guaix-ucm-megaradrp-calibrations-"):
                dir_to_strip = part + '/'  # add the trailing slash to ensure we match only the directory name
                break
    if dir_to_strip is None:
        raise ValueError("Could not find the expected marker in the extracted files.")
    else:
        logger.info(f"Found directory to strip: '{dir_to_strip}'")

    # now we can move the files to the base path,
    # stripping the fixed part of the path
    move_files_strip_fixed(files, base_path=base_path, to_strip=dir_to_strip, overwrite=overwrite, dry_run=dry_run, logger=logger)

    # remove the now empty extracted directory
    extracted_dir = base_path / dir_to_strip
    if extracted_dir.is_dir():
        logger.info(f"Removing extracted directory: {extracted_dir}")
        shutil.rmtree(extracted_dir)
    else:
        if dry_run:
            logger.info(f"Removing extracted directory: {extracted_dir} (dry run, not actually removing)")
        else:
            logger.warning(f"Expected extracted directory not found: {extracted_dir}")

    # remove the downloaded zip file 
    # (this is also valid for the dry run, since in that case
    #  the ZIP file is downloaded to a temporary location)
    zip_path = base_path / local_zip_filename
    if zip_path.is_file():
        logger.info(f"Removing downloaded zip file: {zip_path}")
        zip_path.unlink()
    else:
        if dry_run:
            logger.info(f"Removing downloaded zip file: {zip_path} (dry run, not actually removing)")
        else:
            logger.warning(f"Expected zip file not found: {zip_path}")


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
        help="URL of the calibration data ZIP file to download.",
        type=str,
        default="https://zenodo.org/records/18623771/files/guaix-ucm/megaradrp-calibrations-2026.2.zip",
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

    # download and install the calibration data
    install_calibration_data(
        url=args.zip_url,
        base_path=args.base_path,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        logger=logger,
    )

    # goodbye message
    logger.info(f"Calibration data tree structure created under base path: {args.base_path}")
    console.rule("[bold magenta]Calibration Data Tree Initialization Completed[/bold magenta]")
    if args.dry_run:
        logger.info("This was a dry run, no actual changes were made!")


if __name__ == "__main__":
    main()
