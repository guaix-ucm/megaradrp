#
# Copyright 2015 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#

"""Tests for the bpm mode recipe module."""

from numina.core import DataFrame, ObservationResult
import astropy.io.fits as fits
import numpy as np
import os
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import BadPixelCorrector

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.processing.trimover import OverscanCorrector, TrimImage


class TestRecipe(MegaraBaseRecipe):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_bpm = MasterBPMRequirement()

    def __init__(self):
        super(TestRecipe, self).__init__(version="0.1.0")
        self._MegaraBaseRecipe__flow['TestRecipe'] = [OverscanCorrector,
                                                      TrimImage,
                                                      BadPixelCorrector]

    def run(self, rinput):
        import copy

        directorio = os.path.dirname(os.path.abspath(__file__))

        N = len(rinput.obresult.frames)
        obresult1 = copy.copy(rinput.obresult)
        obresult1.frames = rinput.obresult.frames[:N]

        params = {}

        with rinput.master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        params['biasmap'] = mbias

        reduced1 = self.bias_process_common(obresult1, params)
        fits.writeto(directorio + '/tmp2/reduced_flat.fits', reduced1[0].data, clobber=True)

        try:
            with rinput.master_bpm.open() as hdul:
                bpm = hdul[0].data.copy()
            params['bpm'] = bpm
        except:
            pass

        reduced1 = self.bias_process_common(obresult1, params)
        fits.writeto(directorio + '/tmp2/reduced_flat_bpm.fits', reduced1[0].data, clobber=True)

        return True


def test_bpm_corrector():
    number = 5
    DSHAPE = (2056 * 2, 2048 * 2)

    eq = 0.8 * np.ones(DSHAPE)
    eq[0:15, 0:170] = 0.0

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    names = []

    directorio = os.path.dirname(os.path.abspath(__file__))

    for aux in range(number):
        names.append('%s/tmp2/flat_%s.fits' % (directorio,aux))
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = TestRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(
        filename=open(directorio + '/tmp2/master_bias_data0.fits').name),
                             master_bpm=DataFrame(filename=open(
                                 directorio+ '/tmp2/master_bpm.fits').name))

    recipe.run(ri)


if __name__ == "__main__":
    test_bpm_corrector()
