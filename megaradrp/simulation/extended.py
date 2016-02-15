#
# Copyright 2016 Universidad Complutense de Madrid
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

"""Extended multiwavelength simulation"""


from __future__ import print_function


import numpy as np


def create_dummy_spectrum(wl_in):
    wlc = [0.37, 0.38, 0.39, 0.40, 0.41, 0.42]
    s = 8e-6

    photons_in = np.zeros_like(wl_in)

    for c in wlc:
        photons_in +=  4e4 * np.exp(-0.5*((wl_in-c) / s)**2)

    photons_in += 4.0e4
    return photons_in


# A ThAR arc...
def create_th_ar_arc_spectrum(wl_in):

    lines = [(3719.41400, 58.52137),
             (3803.10800, 54.91138),
             (3828.37100, 56.63661),
             (3839.72400, 46.62416),
             (4019.13100, 40.06343),
             (4071.99600, 42.12776),
             (4131.76200, 50.57991),
             (4158.58400, 68.03227),
             (4200.65300, 48.00637),
             (4277.55800, 78.81951),
             (4348.11900, 80.24873)
             ]

    s = 8e-6
    #wl_in = vph.wltable_interp()
    photons_in = 0.0 * wl_in

    for cwl, flux in lines:
        c = cwl / 1e4
        photons_in +=  flux * 1.5e4 * np.exp(-0.5*((wl_in-c) / s)**2)

    #photons_in += 1.0e2

    return photons_in
