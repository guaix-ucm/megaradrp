#
# Copyright 2015-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Tests for the bpm mode recipe module."""

import astropy.io.fits as fits
from numina.array import combine
from numina.core import DataFrame, ObservationResult
from numina.core.requirements import ObservationResultRequirement
import numina.instrument.assembly as asb

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.testing.create_image import create_sample_files


class DerivedRecipe(MegaraBaseRecipe):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_bpm = MasterBPMRequirement()

    def __init__(self, directorio):
        self.directorio = directorio
        super(DerivedRecipe, self).__init__()

    def run(self, rinput):
        import copy

        N = len(rinput.obresult.frames)
        obresult1 = copy.copy(rinput.obresult)
        obresult1.frames = rinput.obresult.frames[:N]

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        reduced1 = img

        fits.writeto(
            self.directorio + "/reduced_flat.fits", reduced1[0].data, overwrite=True
        )

        fits.writeto(
            self.directorio + "/reduced_flat_bpm.fits", reduced1[0].data, overwrite=True
        )

        return self.create_result()


def test_bpm_corrector():
    import shutil
    from tempfile import mkdtemp

    config_uuid = "4fd05b24-2ed9-457b-b563-a3c618bb1d4c"
    date_obs = "2017-11-09T11:00:00.0"
    directorio = mkdtemp()
    names = create_sample_files(directorio, number=4)

    ob = ObservationResult()
    ob.instrument = "MEGARA"
    ob.mode = "MegaraBiasImage"
    pkg_paths = ["megaradrp.instrument.configs"]
    store = asb.load_paths_store(pkg_paths)
    insmodel = asb.assembly_instrument(store, config_uuid, date_obs, by_key="uuid")
    header = fits.Header()
    header["DATE-OBS"] = date_obs
    header["INSCONF"] = config_uuid
    header["INSTRUME"] = "MEGARA"
    header["VPH"] = "LR-U"
    header["INSMODE"] = "MOS"
    insmodel.configure_with_header(header)
    ob.configuration = insmodel
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = DerivedRecipe(directorio)
    ri = recipe.create_input(
        obresult=ob,
        master_bias=DataFrame(
            filename=open(directorio + "/master_bias_data0.fits").name
        ),
        master_bpm=DataFrame(filename=open(directorio + "/master_bpm.fits").name),
    )

    recipe.run(ri)
    shutil.rmtree(directorio)
