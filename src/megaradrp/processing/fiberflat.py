#
# Copyright 2011-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import datetime
import logging
import sys

if sys.version_info[:2] <= (3, 10):
    datetime.UTC = datetime.timezone.utc

from astropy.io import fits
import numpy
from numina.processing import Corrector
from numina.frame.utils import copy_img

_logger = logging.getLogger(__name__)

# name of the extension with the mask of the fiber flat, written by the fiber flat recipe: the
# pixels of each fiber where the fiber flat is interpolated, or where the light is not corrected
FLAT_MASK_EXTNAME = "FLATMASK"
# name of the mask in the fiber flats made before, with only the bits of the zero order
ZERO_ORDER_MASK_EXTNAME = "ZOMASK"

# bits of the mask of the fiber flat
# zero order of the VPH, interpolated in the fiber flat
FLAT_MASK_ZERO_ORDER = 1
# zero order of the VPH, weak, not interpolated
FLAT_MASK_ZERO_ORDER_WEAK = 2
# absorption windows of the calibration light, interpolated in the fiber flat
FLAT_MASK_ABSORPTION = 4


def flat_mask_hdu(mask):
    """Extension with the mask of the fiber flat, with the meaning of its bits in the header"""
    hdu = fits.ImageHDU(mask, name=FLAT_MASK_EXTNAME)
    hdu.header["BIT1"] = ("ZEROINT", "Zero order of the VPH, interpolated in the fiber flat")
    hdu.header["BIT2"] = ("ZEROWEAK", "Zero order of the VPH, weak, not interpolated")
    hdu.header["BIT3"] = ("ABSINT", "Absorption of the calibration light, interpolated in the fiber flat")
    hdu.header["COMMENT"] = "Pixels of each fiber where the fiber flat is interpolated or does not correct the light"
    return hdu


def find_flat_mask(hdul):
    """The extension with the mask of the fiber flat, also with the name of the older fiber flats, or None"""
    for name in (FLAT_MASK_EXTNAME, ZERO_ORDER_MASK_EXTNAME):
        if name in hdul:
            return hdul[name]
    return None


class CommonFlatCorrector(Corrector):
    """A Node that corrects from fiber flat."""

    def __init__(self, flatdata, datamodel=None, calibid="calibid-unknown", dtype="float32"):

        super(CommonFlatCorrector, self).__init__(datamodel=datamodel, calibid=calibid, dtype=dtype)
        if isinstance(flatdata, fits.HDUList):
            self.corr = flatdata[0].data
        elif isinstance(flatdata, fits.ImageHDU):
            self.corr = flatdata.data
        else:
            self.corr = numpy.asarray(flatdata)

        self.corrmean = self.corr.mean()
        self.flattag = "flat"

    def run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug("correct %s in image %s", self.flattag, imgid)

        # Avoid nan values when divide
        my_mask = self.corr == 0.0
        self.corr[my_mask] = 1.0

        img["primary"].data /= self.corr
        hdr = img["primary"].header

        self.header_update(hdr, imgid)

        return img

    def header_update(self, hdr, imgid):
        hdr["NUM-FLT"] = self.calibid
        hdr["history"] = f"Flat correction {imgid}"
        tnow = datetime.datetime.now(datetime.UTC)
        hdr["history"] = f"Flat correction time {tnow.isoformat()}"


class FiberFlatCorrector(CommonFlatCorrector):
    """A Node that corrects from fiber flat.

    The fibers invalid in the fiber flat are marked invalid in the corrected image.
    The mask of the fiber flat (extension FLATMASK, or ZOMASK in the fiber flats
    made before, with only the zero order) is added to the corrected image as
    FLATMASK, or combined with OR with the one it already has.
    """

    def __init__(self, fiberflat, datamodel=None, dtype="float32"):
        calibid = datamodel.get_imgid(fiberflat)
        super(FiberFlatCorrector, self).__init__(flatdata=fiberflat, datamodel=datamodel, calibid=calibid, dtype=dtype)
        self.flattag = "fiberflat"
        self.invalid_fibers = []
        if isinstance(fiberflat, fits.HDUList) and "FIBERS" in fiberflat:
            hdr = fiberflat["FIBERS"].header
            self.invalid_fibers = sorted(
                int(key[3:6])
                for key, value in hdr.items()
                if key.startswith("FIB") and key.endswith("_V") and value is False
            )
        self.flat_mask = None
        if isinstance(fiberflat, fits.HDUList):
            mask_hdu = find_flat_mask(fiberflat)
            if mask_hdu is not None:
                self.flat_mask = mask_hdu.copy()
                self.flat_mask.name = FLAT_MASK_EXTNAME

    def run(self, img):
        img = super(FiberFlatCorrector, self).run(img)
        if self.invalid_fibers and "FIBERS" in img:
            hdr = img["FIBERS"].header
            for fibid in self.invalid_fibers:
                hdr[f"FIB{fibid:03d}_V"] = (False, "Fiber is invalid")
        if self.flat_mask is not None:
            mask = self.flat_mask.data
            current = find_flat_mask(img)
            if mask.shape != img["primary"].data.shape:
                _logger.warning(
                    "mask of the fiber flat %s does not match the image %s",
                    mask.shape,
                    img["primary"].data.shape,
                )
            elif current is not None:
                current.name = FLAT_MASK_EXTNAME
                current.data = current.data | mask
            else:
                img.append(self.flat_mask.copy())
        return img

    def header_update(self, hdr, imgid):
        hdr["NUM-FIBF"] = self.calibid
        hdr["history"] = f"Fiber flat correction {imgid}"
        tnow = datetime.datetime.now(datetime.UTC)
        hdr["history"] = f"Fiber flat correction time {tnow.isoformat()}"


class Splitter(Corrector):
    """Make a copy of its input"""

    def __init__(self):
        super(Splitter, self).__init__()
        self.out = None

    def run(self, img):
        self.out = self.copy_img(img)
        return img

    def copy_img(self, img):
        return copy_img(img)


class FlipLR(Corrector):
    """Make a copy of its input"""

    def __init__(self):
        super(FlipLR, self).__init__()

    def run(self, img):
        img[0].data = numpy.fliplr(img[0].data)
        return img
