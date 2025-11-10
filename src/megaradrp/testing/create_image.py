#
# Copyright 2015-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
import numina.core
import numpy
from astropy.io import fits as fits
import numpy as np
from numina.core import ObservationResult
from numina.instrument import assembly as asb
from numina.types.dataframe import DataFrame

from megaradrp.datamodel import create_default_fiber_header
from megaradrp.instrument.components.detector import (
    ReadParams,
    MegaraDetectorSat,
    MegaraDetector,
)
from megaradrp.simulation.actions import simulate_flat
from megaradrp.testing.create_ob import generate_bias


def create_empty_img(insmode):
    return create_simple_img(insmode)


def create_simple_img(insmode="LCB", insconf="ca3558e3-e50d-4bbc-86bd-da50a0998a48"):
    hdr = {
        "INSTRUME": "MEGARA",
        "INSCONF": insconf,
        "INSMODE": insmode,
        # "DATE-OBS": "2019-02-21T01:02:02.2",
        "UUID": "410f6ea0-c3df-481c-9820-5cf9e0ed3d91",
        "VPH": "LR-B",
    }

    hdu = fits.PrimaryHDU(data=[[1]])
    for k, v in hdr.items():
        hdu.header[k] = v

    hdrf = create_default_fiber_header(insmode)
    fibers = fits.ImageHDU(header=hdrf, name="FIBERS")

    hdulist = fits.HDUList([hdu, fibers])
    return hdulist


def create_simple_hdul():
    """Create a simple image for testing"""
    prim = fits.PrimaryHDU()
    prim.header["instrume"] = "MEGARA"
    prim.header["VPH"] = "LR-B"
    prim.header["DATE-OBS"] = "2019-02-21T01:02:02.2"
    prim.header["insmode"] = "LCB"
    prim.header["UUID"] = "410f6ea0-c3df-481c-9820-5cf9e0ed3d91"
    fibers = fits.ImageHDU(name="FIBERS")
    fibers.header["CONFID"] = "a908e9d2-7599-41f3-9e9d-91042df01da5"
    simple_img = fits.HDUList([prim, fibers])
    return simple_img


def create_simple_frame(insmode="LCB", insconf="ca3558e3-e50d-4bbc-86bd-da50a0998a48"):
    img = create_simple_img(insmode, insconf)
    frame = numina.core.DataFrame(frame=img)
    return frame


def create_scene_1212(value):
    shape = (623, 4300)
    data = value + numpy.zeros(shape, dtype="float32")
    return data


def create_rss(scene, wlmap):
    hdu = fits.PrimaryHDU(scene)
    hdrf = create_default_fiber_header("LCB")
    fibers = fits.ImageHDU(header=hdrf, name="FIBERS")
    rss_map = fits.ImageHDU(wlmap, name="WLMAP")
    return fits.HDUList([hdu, fibers, rss_map])


def generate_bias_file():
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    ron = 2.0
    gain = 1.0
    bias = 1000.0

    qe = 0.8 * np.ones(DSHAPE)
    qe[0:15, 0:170] = 0.0

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(
        "megara_test_detector",
        DSHAPE,
        OSCAN,
        PSCAN,
        qe=qe,
        dark=(3.0 / 3600.0),
        readpars1=readpars1,
        readpars2=readpars2,
        bins="11",
    )

    return simulate_flat(detector, exposure=1.0, source=5000.0)


def generate_imgs(nimages, nfibers, nsamples):
    imgs = []
    for _ in range(nimages):
        data = numpy.empty((nfibers, nsamples))
        hdu1 = fits.PrimaryHDU(data)
        hdu2 = fits.ImageHDU(name="FIBERS")
        img = fits.HDUList([hdu1, hdu2])
        imgs.append(img)
    return imgs


##############################################################################################


def create_detector():
    DSHAPE = (2056 * 2, 2048 * 2)
    PSCAN = 50
    OSCAN = 50
    qe = 1.0
    dcurrent = 3.0 / 3600

    readpars1 = ReadParams(gain=1.0, ron=2.0, bias=1000.0)
    readpars2 = ReadParams(gain=1.0, ron=2.0, bias=1005.0)

    detector = MegaraDetector(
        "Detector",
        DSHAPE,
        OSCAN,
        PSCAN,
        qe=qe,
        dark=dcurrent,
        readpars1=readpars1,
        readpars2=readpars2,
        bins="11",
    )
    return detector


def create_sample_files(temporary_path, number=5):
    from megaradrp.simulation.actions import simulate_flat
    from megaradrp.instrument.components.detector import ReadParams, MegaraDetectorSat
    from megaradrp.recipes.calibration.bpm import BadPixelsMaskRecipe

    config_uuid = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
    date_obs = "2017-11-09T11:00:00.0"
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50
    ron = 2.0
    gain = 1.0
    bias = 1000.0

    qe = 0.8 * np.ones(DSHAPE)
    qe[0:15, 0:170] = 0.0

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(
        "megara_test_detector",
        DSHAPE,
        OSCAN,
        PSCAN,
        qe=qe,
        dark=(3.0 / 3600.0),
        readpars1=readpars1,
        readpars2=readpars2,
        bins="11",
    )

    source2 = 1.0

    fs = [
        simulate_flat(detector, exposure=1.0, source=5000 * source2)
        for i in range(number)
    ]

    header = fits.Header()
    header["DATE-OBS"] = date_obs
    header["INSCONF"] = config_uuid
    header["INSTRUME"] = "MEGARA"
    header["VPH"] = "LR-U"
    header["INSMODE"] = "MOS"
    for aux in range(len(fs)):
        fits.writeto(
            f"{temporary_path}/flat_{aux}.fits", fs[aux], header=header, overwrite=True
        )

    result = generate_bias(detector, number, temporary_path)
    result.master_bias.frame.writeto(
        f"{temporary_path}/master_bias_data0.fits", overwrite=True
    )  # Master Bias

    ob = ObservationResult()
    ob.instrument = "MEGARA"
    ob.mode = "bias_image"
    pkg_paths = ["megaradrp.instrument.configs"]
    store = asb.load_paths_store(pkg_paths)
    insmodel = asb.assembly_instrument(store, config_uuid, date_obs, by_key="uuid")
    insmodel.configure_with_header(header)
    ob.configuration = insmodel

    names = []
    for aux in range(number):
        names.append(f"{temporary_path}/flat_{aux}.fits")
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = BadPixelsMaskRecipe()
    ri = recipe.create_input(
        obresult=ob,
        master_bias=DataFrame(
            filename=open(temporary_path + "/master_bias_data0.fits").name
        ),
    )
    aux = recipe.run(ri)
    aux.master_bpm.frame.writeto(f"{temporary_path}/master_bpm.fits", overwrite=True)

    return names
