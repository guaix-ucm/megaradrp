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


def main(args=None):
    """Main function to create the calibration data tree structure."""
    # parse command-line options
    parser = argparse.ArgumentParser(description="Create the calibration data tree structure.")
    # positional parameters
    parser.add_argument(
        "--base-path",
        help="Base path where the calibration data tree will be created.",
        type=str,
        default="./calibration_data",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="If set, do not create directories, only print what would be done."
    )
    args = parser.parse_args(args)
    init_tree(base_path=args.base_path, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
