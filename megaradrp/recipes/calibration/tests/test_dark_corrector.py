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
import numpy as np
from numina.core import DataFrame, ObservationResult
import astropy.io.fits as fits
from numina.core.requirements import ObservationResultRequirement
from numina.flow.processing import DarkCorrector, BiasCorrector

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.requirements import MasterBiasRequirement, MasterDarkRequirement
from megaradrp.processing.trimover import OverscanCorrector, TrimImage


class TestRecipe(MegaraBaseRecipe):
    obresult = ObservationResultRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()

    def __init__(self, directorio):
        self.directorio = directorio
        super(TestRecipe, self).__init__(version="0.1.0")
        self._MegaraBaseRecipe__flow['TestRecipe'] = [OverscanCorrector,
                                                      TrimImage, BiasCorrector,
                                                      DarkCorrector]

    def run(self, rinput):
        params = {}

        with rinput.master_bias.open() as hdul:
            mbias = hdul[0].data.copy()

        params['biasmap'] = mbias

        reduced1 = self.bias_process_common(rinput.obresult, params)
        fits.writeto(self.directorio + '/reduced_flat.fits', reduced1[0].data,
                     clobber=True)

        try:
            with rinput.master_dark.open() as hdul:
                dark = hdul[0].data.copy()
            params['dark'] = dark
        except:
            pass

        reduced2 = self.bias_process_common(rinput.obresult, params)

        return reduced1, reduced2


def test_dark_corrector():
    import shutil
    from tempfile import mkdtemp
    from megaradrp.tests.simulation import simulate_flat_fits
    from megaradrp.tests.simulation import ReadParams, MegaraDetectorSat
    from megaradrp.tests.simulation import MegaraImageFactory
    from megaradrp.recipes.calibration.tests.test_bpm_common import \
        generate_bias

    directorio = mkdtemp()
    print ('Path: %s' % directorio)

    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50

    # ron = 2.0
    ron = 2.0
    gain = 1.0
    bias = 1000.0
    eq = 0.8 * np.ones(DSHAPE)

    exptime = (3.0 / 3600.0)
    exposure = 3600
    dark_image = exptime * np.ones(DSHAPE)
    pheader = fits.Header()
    pheader['EXPTIME'] = exptime
    pheader['EXPOSED'] = exptime
    dark_image = fits.PrimaryHDU(dark_image, header=pheader)
    dark_image.writeto('%s/master_dark.fits' % (directorio), clobber=True)

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, eq=eq,
                                 dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    factory = MegaraImageFactory()
    flat_image = simulate_flat_fits(factory, detector, exposure=exposure,
                                    source=5000)
    flat_image.writeto('%s/flat.fits' % (directorio), clobber=True)

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.frames = [DataFrame(filename=open(directorio + '/flat.fits').name)]

    master_bias = generate_bias(detector, 1, directorio)
    master_bias_data = master_bias.biasframe.frame[0].data

    fits.writeto('%s/master_bias.fits' % directorio, master_bias_data,
                 clobber=True)  # Master Bias

    recipe = TestRecipe(directorio)
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(
        filename=open(directorio + '/master_bias.fits').name),
                             master_dark=DataFrame(filename=open(
                                 directorio + '/master_dark.fits').name))

    reduced1, reduced2 = recipe.run(ri)
    shutil.rmtree(directorio)
    assert np.all(
        (reduced1[0].data - (exptime * exposure)) == reduced2[0].data)


if __name__ == "__main__":
    test_dark_corrector()
