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


def generate_bias(detector, number, temporary_path):
    from megaradrp.tests.simulation import simulate_bias
    from megaradrp.recipes.calibration.bias import BiasRecipe

    fs = [simulate_bias(detector) for i in range(number)]
    for aux in range(len(fs)):
        fits.writeto('%s/bias_%s.fits' % (temporary_path, aux), fs[aux],
                     clobber=True)

    fs = ["%s/bias_%s.fits" % (temporary_path, i) for i in range(number)]

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    ob.frames = [DataFrame(filename=f) for f in fs]

    recipe = BiasRecipe()
    ri = recipe.create_input(obresult=ob)
    return recipe.run(ri)

def crear_archivos(temporary_path):
    from megaradrp.tests.simulation import simulate_flat
    from megaradrp.tests.simulation import ReadParams, MegaraDetectorSat
    from megaradrp.recipes.calibration.bpm import BadPixelsMaskRecipe

    number = 5
    PSCAN = 50
    DSHAPE = (2056 * 2, 2048 * 2)
    OSCAN = 50
    ron = 2.0
    gain = 1.0
    bias = 1000.0

    eq = 0.8 * np.ones(DSHAPE)
    eq[0:15, 0:170] = 0.0

    readpars1 = ReadParams(gain=gain, ron=ron, bias=bias)
    readpars2 = ReadParams(gain=gain, ron=ron, bias=bias)

    detector = MegaraDetectorSat(DSHAPE, OSCAN, PSCAN, eq=eq,
                                 dark=(3.0 / 3600.0),
                                 readpars1=readpars1, readpars2=readpars2,
                                 bins='11')

    source2 = 1.0

    fs = [simulate_flat(detector, exposure=1.0, source=5000 * source2) for i in
          range(number)]

    for aux in range(len(fs)):
        fits.writeto('%s/flat_%s.fits' % (temporary_path, aux), fs[aux],
                     clobber=True)

    master_bias = generate_bias(detector, number, temporary_path)
    master_bias_data = master_bias.master_bias.frame[0].data

    fits.writeto('%s/master_bias_data0.fits' % temporary_path,
                 master_bias_data, clobber=True)  # Master Bias

    ob = ObservationResult()
    ob.instrument = 'MEGARA'
    ob.mode = 'bias_image'
    names = []

    for aux in range(number):
        names.append('%s/flat_%s.fits' % (temporary_path, aux))
    ob.frames = [DataFrame(filename=open(nombre).name) for nombre in names]

    recipe = BadPixelsMaskRecipe()
    ri = recipe.create_input(obresult=ob, master_bias=DataFrame(
        filename=open(temporary_path + '/master_bias_data0.fits').name))
    aux = recipe.run(ri)
    fits.writeto('%s/master_bpm.fits' % temporary_path, aux.master_bpm.frame[0].data[1], clobber=True)

    return names

