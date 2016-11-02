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
     'id': '8b'},
    {'nfibers': 21,
     'id': '7b'},
    {'nfibers': 21,
     'id': '6b'},
    {'nfibers': 28,
     'id': '5b'},
    {'nfibers': 28,
     'id': '4b'},
    {'nfibers': 35,
     'id': '3b'},
    {'nfibers': 42,
     'id': '2b'},
    {'nfibers': 77,
     'id': '1b'},
    {'nfibers': 77,
     'id': '0'},
    {'nfibers': 77,
     'id': '1a'},
    {'nfibers': 42,
     'id': '2a'},
    {'nfibers': 35,
     'id': '3a'},
    {'nfibers': 28,
     'id': '4a'},
    {'nfibers': 28,
     'id': '5a'},
    {'nfibers': 21,
     'id': '6a'},
    {'nfibers': 21,
     'id': '7a'},
    {'nfibers': 21,
     'id': '8a'}
]

# Number of fibers in the boxes of the pseudo-slit of the MOS
MEGARA_PSEUDOSLIT_BOXES_MOS = [
    {'nfibers': 7,
     'id': '9b'},
    {'nfibers': 21,
     'id': '8b'},
    {'nfibers': 21,
     'id': '7b'},
    {'nfibers': 21,
     'id': '6b'},
    {'nfibers': 28,
     'id': '5b'},
    {'nfibers': 28,
     'id': '4b'},
    {'nfibers': 35,
     'id': '3b'},
    {'nfibers': 42,
     'id': '2b'},
    {'nfibers': 77,
     'id': '1b'},
    {'nfibers': 77,
     'id': '0'},
    {'nfibers': 77,
     'id': '1a'},
    {'nfibers': 42,
     'id': '2a'},
    {'nfibers': 35,
     'id': '3a'},
    {'nfibers': 28,
     'id': '4a'},
    {'nfibers': 28,
     'id': '5a'},
    {'nfibers': 21,
     'id': '6a'},
    {'nfibers': 21,
     'id': '7a'},
    {'nfibers': 21,
     'id': '8a'},
    {'nfibers': 14,
     'id': '9a'}
]

MEGARA_PSEUDOSLIT_BOXES = {
    'LCB': MEGARA_PSEUDOSLIT_BOXES_LCB,
    'MOS': MEGARA_PSEUDOSLIT_BOXES_MOS
}

# Values for recipe Trace
# Relative threshold for each VPH in LCB
vph_thr = {'LR-I': 0.27,
           'LR-R': 0.37,
           'LR-V': 0.27,
           'LR-Z': 0.27,
           'LR-U': 0.02,
           'HR-I': 0.20,
           }

# Values for recipe Arc
vph_thr_arc = {
    'default':
        {'LR-I': {'min_distance': 10,
                  'threshold': 0.06},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.20},
         'LR-V': {'min_distance': 30,
                  'threshold': 0.19},
         'LR-Z': {'min_distance': 60,
                  'threshold': 0.02},
         'LR-U': {'min_distance': 10,
                  'threshold': 0.02, }
         },
}


# FIXED values for arc calibration
vph_thr_wl_calib = {
    'default': {
        'LR-I': {
            'crval': 7140.0,
            'cdelt': 0.37,
            'crpix': 1.0,
            'npix': 4300
        },
    },
}
