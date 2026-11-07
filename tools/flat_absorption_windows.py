"""Compute the windows of the absorptions in the light of the calibration lamps of MEGARA

The light of the halogen lamps crosses about 4 m of air before the fibers
(measured in the commissioning flats of 2017 with the O2 A and B bands), so
the fiber flats show the O2 and H2O lines, and the lamp has the Na D lines.
The fiber flat recipe skips these windows. They are written to
src/megaradrp/instrument/configs/flat_absorption_windows.ecsv

The transmission of the air path is computed line by line with HITRAN (Voigt
profiles, homogeneous air at 775 hPa and 10 C, 4 m of air with 0.03 mm of
precipitable water), convolved with the resolution of each grating class
(FWHM = wavelength / R) and converted to air wavelengths. A window covers the
wavelengths where the absorption is deeper than 0.1 %, grown by one FWHM at
each side; windows closer than one FWHM are merged.

usage: python tools/flat_absorption_windows.py [O2.par H2O.par]

Without arguments, the lines are downloaded from https://hitran.org/lbl/api
"""

import sys
import urllib.request

import numpy
from astropy.table import Table
from scipy.ndimage import binary_dilation, gaussian_filter1d
from scipy.special import voigt_profile

P_ATM = 0.765
T_K = 283.0
PATH_M = 4.0
WATER_MM = 0.03
THRESHOLD = 0.001
WL_RANGE = (6150.0, 9900.0)
RESOLUTION = {"LR": 6000.0, "MR": 12000.0, "HR": 20000.0}
# The Na D lines of the lamp: +-2 A around each line
NAD = [(5887.95, 5891.95), (5893.92, 5897.92)]

KB = 1.380649e-23
C2 = 1.4387769
AMU = 1.66053907e-27
DNU = 0.002
STEP = 0.02
URL = "https://hitran.org/lbl/api?iso_ids_list={}&numin={:.0f}&numax={:.0f}"


def vac_to_air(w):
    """Morton (2000), wavelengths in Angstrom"""
    s2 = (1e4 / w) ** 2
    return w / (1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2))


def read_par(lines, smin):
    rows = []
    for line in lines:
        if len(line) < 100:
            continue
        s = float(line[15:25])
        if s >= smin:
            rows.append(
                (
                    float(line[3:15]),
                    s,
                    float(line[35:40]),
                    float(line[40:45]),
                    float(line[45:55]),
                    float(line[55:59]),
                    float(line[59:67]),
                )
            )
    return numpy.array(rows)


def tau_of(lines, mass, qexp, column, p_self, nu):
    """Optical depth of a column (molecules / cm2) on the wavenumber grid nu"""
    sigma = numpy.zeros_like(nu)
    for nu0, s296, gair, gself, elow, nair, dair in lines:
        s = (
            s296
            * (296.0 / T_K) ** qexp
            * numpy.exp(-C2 * elow * (1 / T_K - 1 / 296.0))
            * (1 - numpy.exp(-C2 * nu0 / T_K))
            / (1 - numpy.exp(-C2 * nu0 / 296.0))
        )
        nuc = nu0 + dair * P_ATM
        gam = (296.0 / T_K) ** nair * (gair * (P_ATM - p_self) + gself * p_self)
        sig_d = nuc / 2.99792458e8 * numpy.sqrt(KB * T_K / (mass * AMU))
        i0, i1 = numpy.searchsorted(nu, (nuc - 1.0, nuc + 1.0))
        sigma[i0:i1] += s * voigt_profile(nu[i0:i1] - nuc, sig_d, gam)
    return sigma * column


def windows(mask, wl, min_gap):
    edges = numpy.flatnonzero(numpy.diff(numpy.r_[0, mask.astype(int), 0]))
    out = []
    for a, b in zip(edges[::2], edges[1::2] - 1):
        if out and wl[a] - out[-1][1] < min_gap:
            out[-1][1] = wl[b]
        else:
            out.append([wl[a], wl[b]])
    return out


def main():
    numin, numax = 1e8 / (WL_RANGE[1] + 50), 1e8 / (WL_RANGE[0] - 50)
    if len(sys.argv) == 3:
        o2 = open(sys.argv[1]).readlines()
        h2o = open(sys.argv[2]).readlines()
    else:
        o2 = urllib.request.urlopen(URL.format("36,37,38", numin, numax)).read().decode().splitlines()
        h2o = urllib.request.urlopen(URL.format("1,2,3,4", numin, numax)).read().decode().splitlines()
    n_air = P_ATM * 101325 / (KB * T_K) * 1e-6  # cm-3
    o2_column = 0.2095 * n_air * PATH_M * 100.0
    water_column = WATER_MM / (18 * AMU) * 1e-4  # kg/m2 -> molecules/cm2
    es = 611.2 * numpy.exp(17.62 * (T_K - 273.15) / (243.12 + T_K - 273.15))
    p_self = water_column / (PATH_M * 100.0) * KB * T_K * 1e6 / 101325
    p_self = min(p_self, es / 101325)

    nu = numpy.arange(numin, numax, DNU)
    tau = tau_of(read_par(o2, 1e-29), 32.0, 1.0, o2_column, 0.0, nu)
    tau += tau_of(read_par(h2o, 1e-27), 18.0, 1.5, water_column, p_self, nu)
    wl_air = vac_to_air(1e8 / nu[::-1])
    fine = numpy.arange(WL_RANGE[0], WL_RANGE[1], STEP)
    trans = numpy.exp(-numpy.interp(fine, wl_air, tau[::-1]))

    rows = []
    for name, res in RESOLUTION.items():
        absorbed = numpy.zeros_like(fine)
        # the FWHM changes with wavelength: convolve in segments of 250 A
        for lo in numpy.arange(WL_RANGE[0], WL_RANGE[1], 250.0):
            sel = (fine >= lo) & (fine < lo + 250.0)
            fwhm = (lo + 125.0) / res
            conv = gaussian_filter1d(trans, fwhm / 2.3548 / STEP)
            absorbed[sel] = 1 - conv[sel]
        mask = absorbed > THRESHOLD
        fwhm_mean = numpy.mean(fine) / res
        mask = binary_dilation(mask, iterations=int(round(fwhm_mean / STEP)))
        for a, b in windows(mask, fine, fwhm_mean):
            rows.append((name, round(float(a), 2), round(float(b), 2), "O2/H2O"))
        for a, b in NAD:
            rows.append((name, a, b, "NaD"))

    table = Table(rows=rows, names=("resolution", "wl_min", "wl_max", "origin"))
    table.sort(["resolution", "wl_min"])
    table["wl_min"].unit = table["wl_max"].unit = "Angstrom"
    table["resolution"].description = "Resolution class of the gratings: LR, MR or HR"
    table["wl_min"].description = "First wavelength of the window, in air"
    table["wl_max"].description = "Last wavelength of the window, in air"
    table["origin"].description = "O2/H2O: air along the path of the calibration light; NaD: the lamp"
    table.meta["description"] = (
        "Windows of the absorptions in the light of the calibration lamps of MEGARA, skipped by the fiber flat"
    )
    table.meta["method"] = (
        f"HITRAN line by line, {PATH_M} m of air at {P_ATM * 1013.25:.0f} hPa and {T_K - 273.15:.0f} C with "
        f"{WATER_MM} mm of precipitable water, absorption > {100 * THRESHOLD} % at FWHM = wavelength / R "
        f"(R = {RESOLUTION}), grown by one FWHM; Na D +-2 A. Generated by tools/flat_absorption_windows.py"
    )
    out = "src/megaradrp/instrument/configs/flat_absorption_windows.ecsv"
    table.write(out, format="ascii.ecsv", overwrite=True)
    for name in RESOLUTION:
        sub = table[table["resolution"] == name]
        print(f"{name}: {len(sub)} windows, {numpy.sum(sub['wl_max'] - sub['wl_min']):.0f} A in total")


if __name__ == "__main__":
    main()
