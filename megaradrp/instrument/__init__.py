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


# Number of fibers in the boxes of the pseudo-slit of the LCB
MEGARA_PSEUDOSLIT_BOXES_LCB = [
    {'nfibers': 21,
     'name': '8b'},
    {'nfibers': 21,
     'name': '7b'},
    {'nfibers': 21,
     'name': '6b'},
    {'nfibers': 28,
     'name': '5b'},
    {'nfibers': 28,
     'name': '4b'},
    {'nfibers': 35,
     'name': '3b'},
    {'nfibers': 42,
     'name': '2b'},
    {'nfibers': 77,
     'name': '1b'},
    {'nfibers': 77,
     'name': '0'},
    {'nfibers': 77,
     'name': '1a'},
    {'nfibers': 42,
     'name': '2a'},
    {'nfibers': 35,
     'name': '3a'},
    {'nfibers': 28,
     'name': '4a'},
    {'nfibers': 28,
     'name': '5a'},
    {'nfibers': 21,
     'name': '6a'},
    {'nfibers': 21,
     'name': '7a'},
    {'nfibers': 21,
     'name': '8a'}
]

# Number of fibers in the boxes of the pseudo-slit of the MOS
MEGARA_PSEUDOSLIT_BOXES_MOS = [
    {'nfibers': 7,
     'name': '9b'},
    {'nfibers': 21,
     'name': '8b'},
    {'nfibers': 21,
     'name': '7b'},
    {'nfibers': 21,
     'name': '6b'},
    {'nfibers': 28,
     'name': '5b'},
    {'nfibers': 28,
     'name': '4b'},
    {'nfibers': 35,
     'name': '3b'},
    {'nfibers': 42,
     'name': '2b'},
    {'nfibers': 77,
     'name': '1b'},
    {'nfibers': 77,
     'name': '0'},
    {'nfibers': 77,
     'name': '1a'},
    {'nfibers': 42,
     'name': '2a'},
    {'nfibers': 35,
     'name': '3a'},
    {'nfibers': 28,
     'name': '4a'},
    {'nfibers': 28,
     'name': '5a'},
    {'nfibers': 21,
     'name': '6a'},
    {'nfibers': 21,
     'name': '7a'},
    {'nfibers': 21,
     'name': '8a'},
    {'nfibers': 14,
     'name': '9a'}
]

MEGARA_PSEUDOSLIT_BOXES = {
    'LCB': MEGARA_PSEUDOSLIT_BOXES_LCB,
    'MOS': MEGARA_PSEUDOSLIT_BOXES_MOS
}


# Values for recipe Trace
# Relative threshold for each VPH in LCB
vph_thr = {
    'LCB': {
        'LR-I': 0.27,
        'LR-R': 0.10,
        'LR-V': 0.27,
        'LR-Z': 0.27,
        'LR-U': 0.02,
        'HR-I': 0.20,
    },
    'MOS': {
        'LR-I': 0.27,
        'LR-R': 0.10,
        'LR-V': 0.27,
        'LR-Z': 0.27,
        'LR-U': 0.02,
        'HR-I': 0.20,
    },
}


# Values for recipe Arc
vph_thr_arc = {
    'LCB':
        {'LR-I': {'min_distance': 10,
                  'threshold': 0.06},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-V': {'min_distance': 30,
                  'threshold': 0.19},
         'LR-Z': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-U': {'min_distance': 10,
                  'threshold': 0.02, }
         },
    'MOS':
        {'LR-I': {'min_distance': 10,
                  'threshold': 0.06},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-V': {'min_distance': 30,
                  'threshold': 0.19},
         'LR-Z': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-U': {'min_distance': 10,
                  'threshold': 0.02, }
         },
}


# FIXED values for arc calibration
WLCALIB_PARAMS = {
    'LCB': {
        'LR-I': {
            'crval': 7140.0,
            'cdelt': 0.37,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-Z': {
            'crval': 7985.0,
            'cdelt': 0.41,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-R': {
            'crval': 6030.0,
            'cdelt': 0.31,
            'crpix': 1.0,
            'npix': 4300
        },
        'LR-U': {
            'crval': 3620,
            'cdelt': 0.186,
            'crpix': 1.0,
            'npix': 4300
        }
    },
    'MOS': {
        'LR-U': {
            'crval': 3610,
            'cdelt': 0.186,
            'crpix': 1.0,
            'npix': 4400
        }
    },
}
