#
# Copyright 2016-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Focal plane description for MEGARA"""

from __future__ import division


import math
import re
import warnings

import astropy.coordinates

from .ienums import TargetType, BundleType
from megaradrp.processing.hexgrid import connected6


class FocalPlaneConf(object):
    """Configuration of focal plane"""
    def __init__(self, name='LCB'):
        self.name = name
        self.conf_id = "00000000-0000-0000-0000-000000000000"
        bundles = dict()
        if name == 'LCB':
            # 1 LCB + 8 SKY
            bundles[0] = LcbBundleConf()
            for i in range(93, 100 + 1):
                bundles[i] = SkyBundleConf(i)

        elif name == 'MOS':
            # 92 RPs
            for i in range(1, 92 + 1):
                bundles[i] = BundleConf(i, BundleType.RP)
        else:
            raise ValueError(f"name {name} is invalid")
        self.nbundles = len(bundles)
        self.nfibers = sum(bundle.nfibers for bundle in bundles.values())
        self.bundles = bundles
        self.fibers = {}
        self.funit = "mm"

    def attach_fibers(self, fibers):
        bun_fib = {}
        # Fibers in each bundle
        for fibid, ff in fibers.items():
            if ff.bundle_id not in bun_fib:
                bun_fib[ff.bundle_id] = []
            bun_fib[ff.bundle_id].append(ff.fibid)

        for bid, bundle in self.bundles.items():
            # Attach the fibers
            bun_fibers = {fibid: fibers[fibid] for fibid in bun_fib[bid]}
            bundle.attach_fibers(bun_fibers)
        self.fibers = fibers

    @classmethod
    def from_header(cls, hdr):
        """Create a FocalPlaneConf object from map-like header"""

        # defaults = dict()
        # defaults['LCB'] = (9, 623)
        # defaults['MOS'] = (92, 644)

        insmode = hdr.get('INSMODE')
        confid = hdr.get('CONFID', "00000000-0000-0000-0000-000000000000")

        conf = FocalPlaneConf(insmode)
        conf.conf_id = confid

        if conf.nbundles != hdr.get('NBUNDLES'):
            raise ValueError(f'checking NBUNDLES != {conf.nbundles}')
        if conf.nfibers != hdr.get('NFIBERS'):
            raise ValueError(f'checking NFIBERS != {conf.nfibers}')

        conf.funit = hdr.get("FUNIT", "arcsec")
        # Read bundles

        bun_ids = []
        fib_ids = []
        fibers = {}

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

        for fibid in fib_ids:
            ff = FiberConf.from_header(hdr, fibid)
            # if ff.bundle_id not in bun_fib:
            #     bun_fib[ff.bundle_id] = []
            # bun_fib[ff.bundle_id].append(ff.fibid)
            fibers[ff.fibid] = ff

        for bid in bun_ids:
            # Get bundle
            bb = conf.bundles[bid]
            bb.target_priority = hdr["BUN%03d_P" % bid]
            bb.target_name = hdr["BUN%03d_I" % bid]
            bb.target_type = TargetType[hdr["BUN%03d_T" % bid]]
            bb.enabled = hdr.get("BUN%03d_E" % bid, True)
            bb.x = hdr.get("BUN%03d_X" % bid, 0.0)
            bb.y = hdr.get("BUN%03d_Y" % bid, 0.0)
            bb.pa = hdr.get("BUN%03d_O" % bid, 0.0)

        conf.attach_fibers(fibers)

        # Double check
        if conf.name == 'LCB':
            refid = 614
            ref_fiber = conf.fibers[refid]
            if ref_fiber.x < -6:
                # arcsec
                if conf.funit != 'arcsec':
                    print('warning')
            else:
                # mm
                if conf.funit != 'mm':
                    print('warning')
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
        """Return the fibers connected in the IFU"""
        if self.name == 'MOS':
            return []

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

    def nearby_fibers(self, fibid):
        """Obtain the fiber IDs of fibers around fibid"""

        # Find the bundle the fiber belongs to
        try:
            fiber_conf = self.fibers[fibid]
            bundle = self.bundles[fiber_conf.bundle_id]
            result = bundle.nearby(fibid)
            return result
        except KeyError as err:
            # equivalent to
            # raise ValueError from None
            ex = ValueError(f'fibid {fibid} does not exist')
            raise ex from err

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


class BundleConf(object):
    """Description of a bundle"""
    def __init__(self, bundle_id, bundle_type, target_type=TargetType.UNASSIGNED,
                 target_priority=0, target_name='unknown', enabled=True, movable=True):
        self.id = bundle_id
        self.name = "unknown"
        self.bundle_type = bundle_type
        self.target_type = target_type
        self.target_priority = target_priority
        self.target_name = target_name
        self.enabled = enabled
        self.movable = movable
        self.nfibers = 7
        if self.bundle_type == BundleType.LCB:
            self.nfibers = 567
        self.x_fix = 0
        self.y_fix = 0
        self.pa_fix = 0
        self.x = 0
        self.y = 0
        self.pa = 0
        # Experimental
        self.fibers = {}
        self._map1 = {}
        self._map2 = {}

    def attach_fibers(self, fibers):
        self.fibers = fibers
        # Experimental
        cos_30 = math.sqrt(3) / 2
        # This is ~= SPAXEL_SCALE
        dd = 0.44225029050703585
        self._map1 = {}
        self._map2 = {}
        for r in self.fibers.values():
            ii = r.x / (dd * cos_30)
            jj = (r.y / dd - 0.5) - ii * 0.5
            uv = (round(ii), round(jj))
            self._map1[r.fibid] = uv
            self._map2[uv] = r.fibid

    def nearby(self, fibid):
        """Find fiber ids around a central fibid"""
        uv1 = self._map1[fibid]
        result = []
        for con in connected6(*uv1):
            try:
                result.append(self._map2[con])
            except KeyError:
                pass
        return result


class LcbBundleConf(BundleConf):
    """Description of the LCB bundle"""
    def __init__(self, bundle_id=0):
        super(LcbBundleConf, self).__init__(
            bundle_id=bundle_id, bundle_type=BundleType.LCB,
            movable=False
        )
        self.name = "LCB"
        self.nrows = 21
        self.ncols = 27


class SkyBundleConf(BundleConf):
    """Description of a sky bundle"""
    def __init__(self, bundle_id):
        super(SkyBundleConf, self).__init__(
            bundle_id=bundle_id, bundle_type=BundleType.SKY,
            target_type=TargetType.SKY, movable=False
        )


class FiberConf(object):
    """Description of the fiber"""
    def __init__(self, fibid=0, bundle_id=None, inactive=False):
        self.fibid = fibid
        self.name = 'unknown'
        self.bundle_id = bundle_id
        self.inactive = inactive
        self.valid = True
        self.x = 0.0
        self.y = 0.0
        self.coord = None

    @classmethod
    def from_header(cls, hdr, fibid):
        ff = FiberConf(fibid=fibid)
        # Coordinates
        dec = hdr["FIB%03d_D" % fibid]
        ra = hdr["FIB%03d_R" % fibid]
        frame = 'icrs'
        ff.d = dec
        ff.r = ra
        ff.coord = astropy.coordinates.SkyCoord(
            ra, dec, frame=frame, unit='deg'
        )
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
        return ff
