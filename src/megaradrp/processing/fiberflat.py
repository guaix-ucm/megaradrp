#
# Copyright 2011-2023 Universidad Complutense de Madrid
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

# name of the extension with the mask of the zero order of the VPH, written by the fiber flat recipe
ZERO_ORDER_MASK_EXTNAME = "ZOMASK"


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
    The mask of the zero order of the VPH of the fiber flat (extension ZOMASK) is
    added to the corrected image, or combined with OR with the one it already has.
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
        self.zero_order_mask = None
        if isinstance(fiberflat, fits.HDUList) and ZERO_ORDER_MASK_EXTNAME in fiberflat:
            self.zero_order_mask = fiberflat[ZERO_ORDER_MASK_EXTNAME].copy()

    def run(self, img):
        img = super(FiberFlatCorrector, self).run(img)
        if self.invalid_fibers and "FIBERS" in img:
            hdr = img["FIBERS"].header
            for fibid in self.invalid_fibers:
                hdr[f"FIB{fibid:03d}_V"] = (False, "Fiber is invalid")
        if self.zero_order_mask is not None:
            mask = self.zero_order_mask.data
            if mask.shape != img["primary"].data.shape:
                _logger.warning(
                    "zero order mask of the fiber flat %s does not match the image %s",
                    mask.shape,
                    img["primary"].data.shape,
                )
            elif ZERO_ORDER_MASK_EXTNAME in img:
                img[ZERO_ORDER_MASK_EXTNAME].data = img[ZERO_ORDER_MASK_EXTNAME].data | mask
            else:
                img.append(self.zero_order_mask.copy())
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
