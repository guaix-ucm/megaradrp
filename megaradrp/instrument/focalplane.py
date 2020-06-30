#
# Copyright 2016-2020 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Focal plane description for MEGARA"""

from __future__ import division

import re
import enum
import warnings


class FocalPlaneConf(object):
    """Configuration of focal plane"""
    def __init__(self, name=""):
        self.name = name
        self.conf_id = 1
        self.nbundles = 0
        self.nfibers = 0
        self.bundles = {}
        self.fibers = {}
        self.funit = "mm"

    @classmethod
    def from_header(cls, hdr):
        """Create a FocalPlaneConf object from map-like header"""
        conf = FocalPlaneConf()
        defaults = {}
        defaults['LCB'] = (9, 623)
        defaults['MOS'] = (92, 644)

        conf.name = hdr.get('INSMODE')
        insmode = conf.name
        conf.conf_id = hdr.get('CONFID', 1)
        conf.nbundles = hdr.get('NBUNDLES', defaults[insmode][0])
        conf.nfibers = hdr.get('NFIBERS', defaults[insmode][1])
        conf.funit = funit = hdr.get("FUNIT", "arcsec")
        # Read bundles

        bun_ids = []
        fib_ids = []
        bundles = conf.bundles
        fibers = conf.fibers

        # loop over everything, count BUN%03d_P and FIB%03d_B
        pattern1 = re.compile(r"BUN(\d+)_P")
        pattern2 = re.compile(r"FIB(\d+)_B")
        for key in hdr:
            bun_re = pattern1.match(key)
            fib_re = pattern2.match(key)
            if bun_re:
                bun_idx = int(bun_re.group(1))
                bun_ids.append(bun_idx)
            elif fib_re:
                fib_idx = int(fib_re.group(1))
                fib_ids.append(fib_idx)

        for i in bun_ids:
            bb = BundleConf()
            bb.id = i
            bb.target_priority = hdr["BUN%03d_P" % i]
            bb.target_name = hdr["BUN%03d_I" % i]
            bb.target_type = TargetType[hdr["BUN%03d_T" % i]]
            bb.enabled = hdr.get("BUN%03d_E" % i, True)
            bb.x = hdr.get("BUN%03d_X" % i, 0.0)
            bb.y = hdr.get("BUN%03d_Y" % i, 0.0)
            bb.pa = hdr.get("BUN%03d_O" % i, 0.0)
            bb.fibers = {}
            bundles[i] = bb

        for fibid in fib_ids:
            ff = FiberConf()
            ff.fibid = fibid

            # Coordinates
            ff.d = hdr["FIB%03d_D" % fibid]
            ff.r = hdr["FIB%03d_R" % fibid]
            ff.o = 0  # hdr["FIB%03d_O" % fibid]
            # Active
            ff.inactive = not hdr["FIB%03d_A" % fibid]

            # Coordinates XY
            ff.x = hdr["FIB%03d_X" % fibid]
            ff.y = hdr["FIB%03d_Y" % fibid]

            ff.bundle_id = hdr["FIB%03d_B" % fibid]
            ff.name = hdr.get("FIB%03d_N" % fibid, 'unknown')

            ff.w1 = hdr.get("FIB%03dW1" % fibid, None)
            ff.w2 = hdr.get("FIB%03dW2" % fibid, None)

            # Validity
            if ff.inactive:
                ff.valid = False
            else:
                ff.valid = hdr.get("FIB%03d_V" % fibid, True)

            bundles[ff.bundle_id].fibers[ff.fibid] = ff
            fibers[ff.fibid] = ff

        return conf

    @classmethod
    def from_img(cls, img):
        """Create a FocalPlaneConf object from a FITS image"""
        return cls.from_header(img['FIBERS'].header)

    def sky_fibers(self, valid_only=False, ignored_bundles=None):
        result = []
        if ignored_bundles is None:
            ignored_bundles = []

        for bundle in self.bundles.values():
            if bundle.id in ignored_bundles:
                continue
            if bundle.target_type is TargetType.SKY:
                if valid_only:
                    for fib in bundle.fibers.values():
                        if fib.valid:
                            result.append(fib.fibid)
                else:
                    result.extend(bundle.fibers.keys())
        return result

    def connected_fibers(self, valid_only=False):

        if self.name == 'MOS':
            raise ValueError('not working for MOS')

        result = []
        for bundle in self.bundles.values():
            if bundle.target_type is not TargetType.SKY:
                if valid_only:
                    for fib in bundle.fibers.values():
                        if fib.valid:
                            result.append(fib)
                else:
                    result.extend(bundle.fibers.values())
        return result

    def inactive_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if fiber.inactive:
                result.append(fiber.fibid)
        return result

    def active_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if not fiber.inactive:
                result.append(fiber.fibid)
        return result

    def valid_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if fiber.valid:
                result.append(fiber.fibid)
        return result

    def invalid_fibers(self):
        result = []
        for fiber in self.fibers.values():
            if not fiber.valid:
                result.append(fiber.fibid)
        return result

    def spectral_coverage(self):
        lowc = []
        upperc = []
        for fibid, r in self.fibers.items():
            if r.w1:
                lowc.append(r.w1)
            if r.w2:
                upperc.append(r.w2)

        mn = max(lowc)
        nn = min(lowc)

        mx = min(upperc)
        nx = max(upperc)
        return (mn, mx), (nn, nx)

    def bundles_to_table(self):
        """Convert bundles to a Table"""

        import astropy.table

        attrnames = ['id', 'x', 'y', 'pa', 'enabled',
                  'target_type', 'target_priority', 'target_name']
        cnames = ['bundle_id', 'x', 'y', 'pa', 'enabled',
                  'target_type', 'target_priority', 'target_name']
        obj_data = {}
        for a, c in zip(attrnames, cnames):
            obj_data[c] = [getattr(ob, a) for ob in self.bundles.values()]
        result = astropy.table.Table(obj_data, names=cnames)
        result['x'].unit = self.funit
        result['y'].unit = self.funit
        result['pa'].unit = 'deg'
        return result

    def fibers_to_table(self):
        """Convert fibers to a Table"""
        import astropy.table
        attrnames = ['fibid', 'name', 'x', 'y', 'inactive', 'valid',
                     'bundle_id']
        cnames = ['fibid', 'name', 'x', 'y', 'inactive', 'valid',
                  'bundle_id']
        obj_data = {}

        for a, c in zip(attrnames, cnames):
            obj_data[c] = [getattr(ob, a) for ob in self.fibers.values()]
        result = astropy.table.Table(obj_data, names=cnames)
        result['x'].unit = self.funit
        result['y'].unit = self.funit
        return result


class FiberConfs(FocalPlaneConf):
    """Configuration of focal plane

    .. deprecated:: 0.10
            `FiberConfs` is replaced by `FocalPlaneConf`. It will
            be removed in 1.0

    """
    def __init__(self):
        super(FiberConfs, self).__init__()
        warnings.warn("The 'FiberConfs' class was renamed to 'FocalPlaneConf'", DeprecationWarning, stacklevel=2)


class TargetType(enum.Enum):
    """Possible targest in a fiber bundle"""
    SOURCE = 1
    UNKNOWN = 2
    UNASSIGNED = 3
    SKY = 4
    REFERENCE = 5
    # aliases for the other fields
    STAR = 5
    BLANK = 4


class BundleConf(object):
    """Description of a bundle"""
    def __init__(self):
        self.id = 0
        self.target_type = TargetType.UNASSIGNED
        self.target_priority = 0
        self.target_name = 'unknown'
        self.x_fix = 0
        self.y_fix = 0
        self.pa_fix = 0
        self.x = 0
        self.y = 0
        self.pa = 0
        self.enabled = True


class FiberConf(object):
    """Description of the fiber"""
    def __init__(self):
        self.fibid = 0
        self.name = 'unknown'
        self.bundle_id = None
        self.inactive = False
        self.valid = True
        self.x = 0.0
        self.y = 0.0
