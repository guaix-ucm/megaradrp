#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
"""Generate the calibration data tree structure."""

import argparse
from pathlib import Path
import pooch
import shutil


def init_tree(base_path=None, dry_run=False):
    """Create the calibration data tree structure.

    Parameters
    ----------
    base_path : str or Path, optional
        Base path where the calibration data tree will be created.
    dry_run : bool, optional
        If True, do not create directories, only print what would be done.
    """

    if base_path is None:
        raise ValueError("Base path must be provided.")

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

    for calibration_id in calibration_id_list:
        for calibration_directory in calibration_directory_list:
            if calibration_directory in ["MasterBPM", "MasterBias", "MasterSlitFlat", "ReferenceExtinctionTable"]:
                path = base_path / calibration_id / calibration_directory
                print(f"Creating directory: {path}")
                if not dry_run:
                    path.mkdir(parents=True, exist_ok=True)
            elif calibration_directory == "LinesCatalog":
                for arc_lamp in arc_lamp_list:
                    for vph_grating in vph_grating_dict.keys():
                        if vph_grating_dict[vph_grating][arc_lamp]:
                            path = base_path / calibration_id / calibration_directory / arc_lamp / vph_grating
                            print(f"Creating directory: {path}")
                            if not dry_run:
                                path.mkdir(parents=True, exist_ok=True)
            else:
                for ifu_mode in ifu_mode_list:
                    for vph_grating in vph_grating_dict.keys():
                        path = base_path / calibration_id / calibration_directory / ifu_mode / vph_grating
                        print(f"Creating directory: {path}")
                        if not dry_run:
                            path.mkdir(parents=True, exist_ok=True)



def move_files_strip_fixed(files, to_strip=None):
    """Move files to the base path, stripping a fixed part of the path.
    
    Parameters
    ----------
    files : list of str or Path
        List of file paths to move.
    to_strip : str, optional
        Fixed part of the path to strip from the file paths. 
        If None, no stripping is done.
    """
    if to_strip is None:
        raise ValueError("to_strip must be provided.")
    
    for f in files:
        f = Path(f)
        if f.is_dir():
            continue
        # use as_posix to ensure we have a consistent string representation of the path
        rel_text = f.as_posix()
        marker = to_strip.lstrip("./")  # normaliza el comienzo
        # search for the marker in the path
        pos = rel_text.find(marker)
        if pos == -1:
            print(f"Marker '{marker}' not found in path: {rel_text}")
            continue
        # replace only the first occurrence of the marker in the path
        dest = rel_text.replace(marker, "", 1)
        dest = Path(dest)
        # ensure the destination directory exists
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"Moving: {f}  ->  {dest}")
        try:
            f.replace(dest)
        except OSError:
            shutil.copy2(f, dest)
            f.unlink()
        print(f"Moved: {f}  ->  {dest}")


def install_calibration_data(url, base_path):
    """Install the calibration data from the given URL.

    Parameters
    ----------
    url : str
        URL of the calibration data archive to download and extract.
    base_path : str or Path
        Base path where the calibration data will be installed.
    """
    if not isinstance(base_path, (str, Path)):
        raise TypeError("Base path must be a string or a Path object.")
    if isinstance(base_path, str):
        base_path = Path(base_path)

    local_zip_filename = "calibration_data.zip"
    print(f"Downloading and installing calibration data from: {url}")
    files = pooch.retrieve(
        url=url,
        known_hash="md5:3d98803089d4822ee00ed03a4dbba229",
        fname=local_zip_filename,
        path=base_path,
        progressbar=True,
        processor=pooch.Unzip(extract_dir='.'),
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
                dir_to_strip = part
                break
    if dir_to_strip is None:
        raise ValueError("Could not find the expected marker in the extracted files.")
    else:
        print(f"Found directory to strip: {dir_to_strip}")
    
    # now we can move the files to the base path, 
    # stripping the fixed part of the path
    move_files_strip_fixed(files, to_strip=dir_to_strip)

    # remove the now empty extracted directory
    extracted_dir = base_path / dir_to_strip
    if extracted_dir.is_dir():
        print(f"Removing extracted directory: {extracted_dir}")
        shutil.rmtree(extracted_dir)
    else:
        print(f"Expected extracted directory not found: {extracted_dir}")

    # remove the downloaded zip file
    zip_path = base_path / local_zip_filename
    if zip_path.is_file():
        print(f"Removing downloaded zip file: {zip_path}")
        zip_path.unlink()
    else:
        print(f"Expected zip file not found: {zip_path}")


def main(args=None):
    """Main function to create the calibration data tree structure."""
    # parse command-line options
    parser = argparse.ArgumentParser(description="Create the calibration data tree structure.")
    # positional parameters
    parser.add_argument(
        "--base-path",
        help="Base path where the calibration data tree will be created.",
        type=str,
        default="./calibrations",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="If set, do not create directories, only print what would be done."
    )
    args = parser.parse_args(args)

    # create the directory tree structure
    init_tree(base_path=args.base_path, dry_run=args.dry_run)

    # download and install the calibration data
    if not args.dry_run:
        install_calibration_data(
            url="https://zenodo.org/records/18623771/files/guaix-ucm/megaradrp-calibrations-2026.2.zip",
            base_path=args.base_path,
        )


if __name__ == "__main__":
    main()
