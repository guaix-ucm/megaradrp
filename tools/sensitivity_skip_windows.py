"""Compute the windows skipped by the smoothing of the sensitivity of the standard star recipes

The smoothing of the sensitivity (megaradrp.processing.fluxcalib.smooth_sensitivity)
skips the telluric bands of the atmosphere and the strong lines of the standard
stars. The windows are written to
src/megaradrp/instrument/configs/sensitivity_skip_windows.ecsv

Telluric bands: from a model transmission of the atmosphere (for example, from
the ESO SkyCalc sky model, in nm, vacuum), degraded to the resolution of the LR
gratings (R = 6000) and to the broadening of the spectrum of the star by the
recipe (sigma_resolution, about 7 Angstrom in LR-Z), after dividing it by its
continuum (an upper envelope over 500 Angstrom): the continuum extinction is
corrected with the extinction curve. A window covers the wavelengths where the
transmission is below 0.97; windows closer than 15
Angstrom are merged, and the wavelengths are converted to air.

Bluer than 6000 Angstrom there are no telluric bands deeper than the threshold:
the envelope of the steep continuum of the blue (Rayleigh) gives false windows
there, so the search starts at 6000 Angstrom.

The broad H2O band at the red end (beyond 9250 Angstrom, the end of LR-Z) is
not skipped: masking it leaves the spline without data at the end of the
grating, and the fall of the response of the instrument there cannot be
separated from the absorption with a single standard. The sensitivity keeps
the telluric absorption of the standard in that band.

Stellar lines: Balmer, Ca II H and K, the Ca II triplet and the Paschen lines
out of the H2O band, wide enough for the broadened lines.

usage: python tools/sensitivity_skip_windows.py TRANSMISSION
    TRANSMISSION: two columns, wavelength (nm, vacuum) and transmission
"""

import sys

import numpy
from astropy.table import Table
from scipy.ndimage import gaussian_filter1d, maximum_filter1d

R_LR = 6000.0
BROADENING_AA = 7.0
THRESHOLD = 0.97
MERGE_AA = 15.0
RED_LIMIT_AA = 9250.0
CONTINUUM_AA = 500.0
WL_RANGE = (6000.0, 10000.0)
STELLAR = [
    (3920.0, 3945.0, "Ca II K"),
    (3955.0, 3985.0, "Ca II H, H epsilon"),
    (4070.0, 4135.0, "H delta"),
    (4310.0, 4370.0, "H gamma"),
    (4830.0, 4895.0, "H beta"),
    (6533.0, 6593.0, "H alpha"),
    (8490.0, 8510.0, "Ca II 8498"),
    (8530.0, 8555.0, "Ca II 8542"),
    (8590.0, 8606.0, "Pa 14"),
    (8650.0, 8675.0, "Ca II 8662"),
    (8740.0, 8760.0, "Pa 12"),
    (8855.0, 8872.0, "Pa 11"),
]


def vac_to_air(w):
    """Morton (2000), wavelengths in Angstrom"""
    s2 = (1e4 / w) ** 2
    return w / (1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2))


def telluric_windows(wl_nm, trans, input_resolution):
    wl = wl_nm * 10.0
    dlog = numpy.median(numpy.diff(numpy.log(wl)))
    # only the bands: the continuum extinction (Rayleigh, aerosols, ozone) is corrected with the extinction curve;
    # computed with all the data, and the range is selected after, to avoid the effects of the borders
    npix = int(CONTINUUM_AA / numpy.median(wl) / dlog)
    continuum = gaussian_filter1d(maximum_filter1d(trans, npix), npix / 4)
    trans = trans / numpy.maximum(continuum, 1e-6)
    sigma_res = numpy.sqrt(1 / R_LR**2 - 1 / input_resolution**2) / 2.3548
    smooth = gaussian_filter1d(trans, sigma_res / dlog)
    smooth = gaussian_filter1d(smooth, BROADENING_AA / numpy.median(wl) / dlog)
    sel = (wl > WL_RANGE[0]) & (wl < WL_RANGE[1])
    wl, smooth = wl[sel], smooth[sel]
    low = smooth < THRESHOLD
    windows = []
    i = 0
    while i < low.size:
        if low[i]:
            j = i
            while j < low.size and low[j]:
                j += 1
            windows.append([wl[i], wl[j - 1]])
            i = j
        else:
            i += 1
    merged = []
    for w1, w2 in windows:
        if merged and w1 - merged[-1][1] < MERGE_AA:
            merged[-1][1] = w2
        else:
            merged.append([w1, w2])
    result = []
    for w1, w2 in merged:
        a1, a2 = vac_to_air(w1), vac_to_air(w2)
        if a1 < RED_LIMIT_AA:
            result.append((round(a1, 1), round(min(a2, RED_LIMIT_AA), 1)))
    return result


def main(argv):
    path = argv[1]
    data = numpy.loadtxt(path)
    wl_nm, trans = data[:, 0], data[:, 1]
    input_resolution = 20000.0
    rows = [(w1, w2, "telluric", "O2/H2O") for w1, w2 in telluric_windows(wl_nm, trans, input_resolution)]
    rows += [(w1, w2, "stellar", name) for w1, w2, name in STELLAR]
    rows.sort()
    table = Table(rows=rows, names=["wl_min", "wl_max", "origin", "feature"], dtype=[float, float, str, str])
    table["wl_min"].unit = table["wl_max"].unit = "Angstrom"
    table["wl_min"].description = "First wavelength of the window, in air"
    table["wl_max"].description = "Last wavelength of the window, in air"
    table["origin"].description = (
        "telluric: absorption bands of the atmosphere; stellar: strong lines of the standard stars"
    )
    table["feature"].description = "Band or line"
    table.meta["description"] = "Windows skipped by the smoothing of the sensitivity of the standard star recipes"
    table.meta["method"] = (
        f"telluric: transmission of {path.split('/')[-1]} divided by its continuum ({CONTINUUM_AA} A), "
        f"degraded to R = {R_LR:.0f} and a broadening of "
        f"{BROADENING_AA} A, below {THRESHOLD}, merged closer than {MERGE_AA} A, in air, not beyond {RED_LIMIT_AA} A "
        "(the broad H2O band at the red end stays in the sensitivity); stellar: Balmer, Ca II, Paschen. "
        "Generated by tools/sensitivity_skip_windows.py"
    )
    out = "src/megaradrp/instrument/configs/sensitivity_skip_windows.ecsv"
    table.write(out, format="ascii.ecsv", overwrite=True)
    for row in table:
        print(f"{row['wl_min']:8.1f} {row['wl_max']:8.1f} {row['origin']:9s} {row['feature']}")


if __name__ == "__main__":
    main(sys.argv)
