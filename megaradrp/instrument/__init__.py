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


# Values for recipe Trace
# Relative threshold for each VPH in LCB
vph_thr = {
    'LCB': {
        'LR-U': 0.02,
        'LR-B': 0.02,
        'LR-V': 0.27,
        'LR-R': 0.10,
        'LR-I': 0.27,
        'LR-Z': 0.27,
        'HR-I': 0.20,
    },
    'MOS': {
        'LR-U': 0.02,
        'LR-B': 0.02,
        'LR-V': 0.27,
        'LR-R': 0.10,
        'LR-I': 0.27,
        'LR-Z': 0.27,
        'HR-I': 0.20,
    },
}


# Values for recipe Arc
vph_thr_arc = {
    'LCB':
        {'LR-U': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-B': {'min_distance': 10,
                  'threshold': 0.05},
         'LR-V': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-I': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-Z': {'min_distance': 10,
                  'threshold': 0.00},
         'HR-R': {'min_distance': 10,
                  'threshold': 0.00},
         'HR-I': {'min_distance': 10,
                  'threshold': 0.00}
         },
    'MOS':
        {'LR-U': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-B': {'min_distance': 10,
                  'threshold': 0.05},
         'LR-V': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-I': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-Z': {'min_distance': 10,
                  'threshold': 0.00},
         'HR-R': {'min_distance': 10,
                  'threshold': 0.00},
         'HR-I': {'min_distance': 10,
                  'threshold': 0.00}
         },
}


# FIXED values for arc calibration
WLCALIB_PARAMS = {
    'LCB': {
        'LR-U': {
            'crval': 3640,
            'cdelt': 0.186,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-B': {
            'crval': 4280.0,
            'cdelt': 0.23,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-V': {
            'crval': 5060.0,
            'cdelt': 0.27,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-R': {
            'crval': 6030.0,
            'cdelt': 0.31,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-I': {
            'crval': 7140.0,
            'cdelt': 0.37,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-Z': {
            'crval': 7960.0,
            'cdelt': 0.41,
            'crpix': 1.0,
            'npix': 4300
        },
        'HR-R': {
            'crval': 6390.0,
            'cdelt': 0.0974,
            'crpix': 1.0,
            'npix': 4300
        },
        'HR-I': {
            'crval': 8350.0,
            'cdelt': 0.13,
            'crpix': 1.0,
            'npix': 4300
        },
    },
    'MOS': {
        'LR-U': {
            'crval': 3640,
            'cdelt': 0.186,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-B': {
            'crval': 4280.0,
            'cdelt': 0.23,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-V': {
            'crval': 5060.0,
            'cdelt': 0.27,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-R': {
            'crval': 6030.0,
            'cdelt': 0.31,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-I': {
            'crval': 7140.0,
            'cdelt': 0.37,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-Z': {
            'crval': 7960.0,
            'cdelt': 0.41,
            'crpix': 1.0,
            'npix': 4300
        },
        'HR-R': {
            'crval': 6390.0,
            'cdelt': 0.0974,
            'crpix': 1.0,
            'npix': 4300
        },
        'HR-I': {
            'crval': 8350.0,
            'cdelt': 0.13,
            'crpix': 1.0,
            'npix': 4300
        },
    },
}
