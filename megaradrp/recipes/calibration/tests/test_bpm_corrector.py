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
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import BadPixelCorrector

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterBPMRequirement
from megaradrp.processing.trimover import OverscanCorrector, TrimImage

from megaradrp.recipes.calibration.tests.test_bpm_common import crear_archivos

class TestRecipe(MegaraBaseRecipe):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_bpm = MasterBPMRequirement()

    def __init__(self, directorio):
        self.directorio = directorio
        super(TestRecipe, self).__init__(version="0.1.0")
        self._MegaraBaseRecipe__flow['TestRecipe'] = [OverscanCorrector,
                                                      TrimImage,
                                                      BadPixelCorrector]

    def run(self, rinput):
        import copy

        N = len(rinput.obresult.frames)
        obresult1 = copy.copy(rinput.obresult)
        obresult1.frames = rinput.obresult.frames[:N]

        params = {}

        with rinput.master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        params['biasmap'] = mbias

        reduced1 = self.bias_process_common(obresult1, params)
        fits.writeto(self.directorio + '/reduced_flat.fits', reduced1[0].data,
                     clobber=True)

        try:
            with rinput.master_bpm.open() as hdul:
                bpm = hdul[0].data.copy()
            params['bpm'] = bpm
        except:
            pass

        reduced1 = self.bias_process_common(obresult1, params)
        fits.writeto(self.directorio + '/reduced_flat_bpm.fits',
                     reduced1[0].data, clobber=True)

        return True


def test_bpm_corrector():
    import shutil
    from tempfile import mkdtemp

    directorio = mkdtemp()
    names = crear_archivos(directorio)

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = TestRecipe(directorio)
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(
        filename=open(directorio + '/master_bias_data0.fits').name),
                             master_bpm=DataFrame(filename=open(
                                 directorio + '/master_bpm.fits').name))

    recipe.run(ri)
    shutil.rmtree(directorio)

if __name__ == "__main__":
    test_bpm_corrector()
