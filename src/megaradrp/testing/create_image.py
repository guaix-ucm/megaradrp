#
# Copyright 2015-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
import warnings

import numina.core
import numpy
from astropy.io import fits as fits
import numpy as np
from numina.types.dataframe import DataFrame

from megaradrp.datamodel import create_default_fiber_header
from megaradrp.simulation.actions import simulate_flat
from megaradrp.testing.common import DSHAPE, make_detector, make_obresult, raw_header
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
    """A simulated flat, of a detector with dead pixels

    .. deprecated:: 0.23
            The tests use the fixture `simulated_flats`. It will be
            removed in 0.24

    """
    warnings.warn("'generate_bias_file' is deprecated, it will be removed in 0.24", DeprecationWarning, stacklevel=2)
    qe = 0.8 * np.ones(DSHAPE)
    qe[0:15, 0:170] = 0.0
    detector = make_detector(qe=qe)
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
    return make_detector(qe=1.0, bias=(1000.0, 1005.0), saturation=False)


def create_sample_files(temporary_path, number=5):
    """Simulated flats, master bias and mask of bad pixels, in temporary_path

    .. deprecated:: 0.23
            The tests use the fixtures `simulated_flats` and `master_bpm`.
            It will be removed in 0.24

    """
    warnings.warn("'create_sample_files' is deprecated, it will be removed in 0.24", DeprecationWarning, stacklevel=2)
    from megaradrp.recipes.calibration.bpm import BadPixelsMaskRecipe

    qe = 0.8 * np.ones(DSHAPE)
    qe[0:15, 0:170] = 0.0
    detector = make_detector(qe=qe)

    header = raw_header()
    names = [f"{temporary_path}/flat_{idx}.fits" for idx in range(number)]
    for name in names:
        fits.writeto(name, simulate_flat(detector, exposure=1.0, source=5000.0), header=header, overwrite=True)

    result = generate_bias(detector, number, temporary_path)
    result.master_bias.frame.writeto(f"{temporary_path}/master_bias_data0.fits", overwrite=True)  # Master Bias

    ob = make_obresult("bias_image", names, header)
    recipe = BadPixelsMaskRecipe()
    ri = recipe.create_input(
        obresult=ob,
        master_bias=DataFrame(filename=f"{temporary_path}/master_bias_data0.fits"),
    )
    aux = recipe.run(ri)
    aux.master_bpm.frame.writeto(f"{temporary_path}/master_bpm.fits", overwrite=True)

    return names
