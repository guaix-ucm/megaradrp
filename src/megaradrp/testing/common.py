#
# Copyright 2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Simulated detector, headers and observation results for the tests"""

import astropy.io.fits as fits
from numina.core import DataFrame, ObservationResult
from numina.instrument.assembly import configure_obresult

from megaradrp.instrument.components.detector import MegaraDetector, MegaraDetectorSat, ReadParams
from megaradrp.loader import load_drp

# The instrument configuration of the test component, and a date where it is valid
TEST_INSCONF = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
TEST_DATE_OBS = "2017-11-09T11:00:00.0"

# Shape of the data of the detector, and size of the overscan and prescan
DSHAPE = (2056 * 2, 2048 * 2)
OSCAN = 50
PSCAN = 50


def make_detector(qe=1.0, dark=3.0 / 3600.0, ron=2.0, gain=1.0, bias=(1000.0, 1000.0), saturation=True):
    """A simulated detector of MEGARA

    Parameters
    ----------
    qe : float or array
        Quantum efficiency, an array of shape DSHAPE or a number.
    dark : float
        Dark current, in counts per second.
    ron, gain : float
        Readout noise and gain of both amplifiers.
    bias : tuple of float
        Bias level of the two amplifiers.
    saturation : bool
        A detector that saturates (MegaraDetectorSat) or not (MegaraDetector).
    """
    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias[0])
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias[1])
    detclass = MegaraDetectorSat if saturation else MegaraDetector
    return detclass(
        "megara_test_detector",
        DSHAPE,
        OSCAN,
        PSCAN,
        qe=qe,
        dark=dark,
        readpars1=readpars1,
        readpars2=readpars2,
        bins="11",
    )


def raw_header(vph="LR-U", insmode="MOS", insconf=TEST_INSCONF, date_obs=TEST_DATE_OBS):
    """The keywords of a raw image used to select the instrument configuration"""
    header = fits.Header()
    header["DATE-OBS"] = date_obs
    header["INSCONF"] = insconf
    header["INSTRUME"] = "MEGARA"
    header["VPH"] = vph
    header["INSMODE"] = insmode
    return header


def make_obresult(mode, filenames=(), header=None):
    """An observation result of MEGARA with the images in filenames

    The instrument configuration is selected and configured as in numina
    run, with the header (raw_header() by default) as the first image.
    """
    if header is None:
        header = raw_header()
    # The images may not exist, the configuration is taken from the header
    sample = ObservationResult(instrument="MEGARA")
    sample.frames = [DataFrame(frame=fits.HDUList([fits.PrimaryHDU(header=header)]))]
    configure_obresult(sample, load_drp())

    ob = ObservationResult()
    ob.instrument = "MEGARA"
    ob.mode = mode
    ob.configuration = sample.configuration
    ob.profile = sample.profile
    ob.frames = [DataFrame(filename=str(name)) for name in filenames]
    return ob
