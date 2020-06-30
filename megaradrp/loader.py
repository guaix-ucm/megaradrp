#
# Copyright 2014-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Load MEGARA DRP"""

from numina.core import drp_load
import numina.core.config as cfg


class MegaraDrpLoader(object):
    """Custom loader class

    This class modifies the rawimage field of the observing modes
    of MEGARA
    """
    @staticmethod
    def mode_loader(mode_node):
        import megaradrp.datamodel as DM
        if 'rawimage' in mode_node:
            rname = mode_node['rawimage']
            mode_node['rawimage'] = DM.MegaraDataType[rname]
        return mode_node


def load_drp():
    """Entry point to load MEGARA DRP."""
    return drp_load('megaradrp', 'drp.yaml', confclass=MegaraDrpLoader)


def is_fits_megara(pathname):
    "Check is any FITS"
    import astropy.io.fits as fits
    # FIXME: incomplete
    if pathname.endswith('.fits') or pathname.endswith('.fits.gz'):
        with fits.open(pathname) as hdulist:
            prim = hdulist[0].header
            instrument = prim.get("INSTRUME", "unknown")
            if instrument == "MEGARA":
                return True
    else:
        return False


@cfg.describe.register('image/fits', is_fits_megara, priority=15)
def describe_fits_megara(pathname):
    import megaradrp.datamodel as DM
    import astropy.io.fits as fits
    with fits.open(pathname) as hdulist:
        return DM.describe_hdulist_megara(hdulist)


@cfg.check.register('MEGARA')
def check_obj_megara(obj, astype=None, level=None):
    import megaradrp.datamodel as DM
    return DM.check_obj_megara(obj, astype=astype, level=level)
