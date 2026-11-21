#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Measure the offset of the traces of a TraceMap or ModelMap in an image"""

import argparse
import json
import logging
import sys

from astropy.io import fits
import numina.types.structured as structured

from megaradrp.processing.traceoffset import measure_extraction_offset

_logger = logging.getLogger(__name__)


def plot_offset(result, path):
    """chi-square of the global search and offset of each band"""
    from matplotlib.figure import Figure

    fig = Figure(figsize=(10, 4), layout="constrained")
    ax1, ax2 = fig.subplots(1, 2)
    ax1.plot(result.grid, (result.curve - result.curve.min()) / result.reduced_chi2, ".-")
    ax1.axvline(result.grid[result.curve.argmin()], color="gray", ls=":")
    ax1.set_yscale("symlog")
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel("c0 with c1 = 0 (pixels)")
    ax1.set_ylabel("chi-square - minimum (reduced units)")
    centres = [0.5 * (band["x0"] + band["x1"]) for band in result.bands]
    ax2.errorbar(centres, [band["c0"] for band in result.bands], [band["error"] for band in result.bands], fmt="o")
    ax2.axhline(result.coef[0], color="gray", ls=":")
    ax2.set_xlabel("column")
    ax2.set_ylabel("c0 of the band, with the fitted c1 (pixels)")
    fig.suptitle(
        f"offset [{result.coef[0]:.3f}, {result.coef[1]:.3g}], {'reliable' if result.reliable else 'NOT reliable'}"
    )
    fig.savefig(path)


def main(args=None):
    parser = argparse.ArgumentParser(
        description="Measure the offset of the traces of a TraceMap or ModelMap in an image (extraction_offset)"
    )
    parser.add_argument("image", help="FITS image after the basic processing (e.g. reduced_image.fits), in electrons")
    parser.add_argument("traces", help="JSON file with a TraceMap or a ModelMap")
    parser.add_argument("--nbands", type=int, default=16, help="number of bands of columns (default: 16)")
    parser.add_argument("--band-width", type=int, default=64, help="columns of each band (default: 64)")
    parser.add_argument(
        "--search", type=float, nargs=2, default=(-8.0, 8.0), metavar=("MIN", "MAX"), help="range of the search of c0"
    )
    parser.add_argument("--step", type=float, default=0.25, help="step of the search of c0 (default: 0.25)")
    parser.add_argument("--ron", type=float, help="readout noise in electrons (default: RDNOISE1 of the image, or 3.4)")
    parser.add_argument("--no-slope", action="store_true", help="fit only c0")
    parser.add_argument("--plot", help="PDF or PNG file with the chi-square of the search and the offset of each band")
    parser.add_argument("--json", help="write the result to this JSON file")
    args = parser.parse_args(args)

    logging.basicConfig(level=logging.WARNING, format="%(levelname)s: %(message)s")
    with fits.open(args.image) as hdul:
        data = hdul[0].data.astype(float)
        header = hdul[0].header
    if header.get("BUNIT", "ELECTRON").upper() != "ELECTRON":
        _logger.warning("the image is in %s, not in electrons: the weights of the fit are approximate", header["BUNIT"])
    ron = args.ron if args.ron is not None else header.get("RDNOISE1", 3.4)
    traces = structured.open(args.traces)

    result = measure_extraction_offset(
        data,
        traces,
        nbands=args.nbands,
        band_width=args.band_width,
        search=tuple(args.search),
        step=args.step,
        ron=ron,
        fit_slope=not args.no_slope,
    )

    print(f"extraction_offset: [{result.coef[0]:.4f}, {result.coef[1]:.6g}]")
    print(f"uncertainty: [{result.error[0]:.4f}, {result.error[1]:.3g}]")
    print(f"reliable: {result.reliable}" + ("" if result.reliable else f" ({result.message})"))
    print(f"alias contrast: {result.alias_contrast:.1f}, bootstrap: {result.alias_bootstrap:.1%}")
    print(f"fibers with signal: {result.nsignal}, width scale: {result.sigma_scale:.2f}")
    print("offset of each band (columns: c0 +- error):")
    for band in result.bands:
        print(f"  {band['x0']:4d}-{band['x1']:4d}: {band['c0']:+.3f} +- {band['error']:.3f}")

    if args.plot:
        plot_offset(result, args.plot)
    if args.json:
        state = {
            "coef": result.coef,
            "error": result.error,
            "reliable": result.reliable,
            "message": result.message,
            "alias_contrast": result.alias_contrast,
            "alias_bootstrap": result.alias_bootstrap,
            "nsignal": result.nsignal,
            "sigma_scale": result.sigma_scale,
            "reduced_chi2": result.reduced_chi2,
            "bands": result.bands,
        }
        with open(args.json, "w") as fd:
            json.dump(state, fd, indent=1)
    return 0 if result.reliable else 1


if __name__ == "__main__":
    sys.exit(main())
