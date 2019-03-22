#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#


MEGARA_PLATESCALE = 1.2120 # arcsec / mm

MEGARA_IAA = -163.854 # deg


# Values for recipe Trace
# Relative threshold for each VPH in LCB
vph_thr = {
    'LCB': {
        'LR-U': 0.02,
        'LR-B': 0.02,
        'LR-V': 0.02,
        'LR-R': 0.10,
        'LR-I': 0.27,
        'LR-Z': 0.27,
        'MR-U': 0.05,
        'MR-UB': 0.05,
        'MR-B': 0.05,
        'MR-VR': 0.05,
        'MR-G': 0.05,
        'MR-V': 0.02,
        'MR-R': 0.05,
        'MR-RI': 0.05,
        'MR-I': 0.05,
        'MR-Z': 0.05,
        'HR-I': 0.20,
    },
    'MOS': {
        'LR-U': 0.02,
        'LR-B': 0.02,
        'LR-V': 0.02,
        'LR-R': 0.10,
        'LR-I': 0.27,
        'LR-Z': 0.27,
        'MR-U': 0.05,
        'MR-UB': 0.05,
        'MR-B': 0.05,
        'MR-G': 0.05,
        'MR-V': 0.02,
        'MR-VR': 0.05,
        'MR-R': 0.10,
        'MR-RI': 0.10,
        'MR-I': 0.05,
        'MR-Z': 0.05,
        'HR-I': 0.20,
    },
}


# Values for recipe Arc
vph_thr_arc = {
    'LCB':
        {'LR-U': {'min_distance': 10,
                  'threshold': 0.02},
         'LR-B': {'min_distance': 10,
                  'threshold': 0.00},
         'LR-V': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.00},
         'LR-I': {'min_distance': 10,
                  'threshold': 0.00},
         'LR-Z': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-U': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-UB': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-B': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-G': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-V': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-VR': {'min_distance': 10,
                   'threshold': 0.00},
         'MR-R': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-RI': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-I': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-Z': {'min_distance': 10,
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
                  'threshold': 0.00},
         'LR-V': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-R': {'min_distance': 10,
                  'threshold': 0.00},
         'LR-I': {'min_distance': 10,
                  'threshold': 0.01},
         'LR-Z': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-U': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-UB': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-B': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-G': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-V': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-VR': {'min_distance': 10,
                   'threshold': 0.00},
         'MR-R': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-RI': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-I': {'min_distance': 10,
                  'threshold': 0.00},
         'MR-Z': {'min_distance': 10,
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
            'crval': 3620,
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
        'MR-U': {
            'crval': 3905.0,
            'cdelt': 0.089,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-UB': {
            'crval': 4210.0,
            'cdelt': 0.10,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-B': {
            'crval': 4568.0,
            'cdelt': 0.11,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-G': {
            'crval': 4944.0,
            'cdelt': 0.122,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-V': {
            'crval': 5375.0,
            'cdelt': 0.132,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-VR': {
            'crval': 5850.0,
            'cdelt': 0.145,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-R': {
            'crval': 6210.0,
            'cdelt': 0.16,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-RI': {
            'crval': 6735.0,
            'cdelt': 0.17,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-I': {
            'crval': 7360.0,
            'cdelt': 0.1845,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-Z': {
            'crval': 8770.0,
            'cdelt': 0.225,
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
            'crval': 3620,
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
        'MR-U': {
            'crval': 3905.0,
            'cdelt': 0.089,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-UB': {
            'crval': 4210.0,
            'cdelt': 0.10,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-B': {
            'crval': 4568.0,
            'cdelt': 0.11,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-G': {
            'crval': 4944.0,
            'cdelt': 0.122,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-V': {
            'crval': 5375.0,
            'cdelt': 0.132,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-VR': {
            'crval': 5850.0,
            'cdelt': 0.145,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-R': {
            'crval': 6210.0,
            'cdelt': 0.16,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-RI': {
            'crval': 6735.0,
            'cdelt': 0.17,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-I': {
            'crval': 7360.0,
            'cdelt': 0.1845,
            'crpix': 1.0,
            'npix': 4300
        },
        'MR-Z': {
            'crval': 8770.0,
            'cdelt': 0.225,
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
