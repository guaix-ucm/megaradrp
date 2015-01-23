#
# Copyright 2011-2014 Universidad Complutense de Madrid
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


from astropy.io import fits
import numpy as np

from numina.core import BaseRecipeAutoQC as MegaraBaseRecipe  # @UnusedImport


# row / column
_binning = {'11': [1, 1], '21': [1, 2], '12': [2, 1], '22': [2, 2]}
_direc = ['normal', 'mirror']


def create(image, direction='normal', bins='11'):
    '''Create a image with overscan for testing.'''

    if direction not in _direc:
        raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

    if direction == 'normal':
        direcfun = lambda x: x
    else:
        direcfun = np.fliplr

    if bins not in _binning:
        raise ValueError("%s must be one if '11', '12', '21, '22'" % bins)

    bng = _binning[bins]

    nr = 2056 / bng[0]
    nc = 2048 / bng[1]

    nr2 = 2 * nr
    nc2 = 2 * nc

    oscan1 = 50 / bng[0]
    oscan2 = oscan1 * 2

    psc1 = 50 / bng[0]
    psc2 = 2 * psc1

    fshape = (nr2 + oscan2, nc2 + psc2)

    # Row block 1
    rb1 = slice(0, nr)
    rb1m = slice(nr, nr + oscan1)
    # Row block 2
    rb2 = slice(nr + oscan2, nr2 + oscan2)
    rb2m = slice(nr + oscan1, nr + oscan2)
    # Col block
    cb = slice(psc1, nc2 + psc1)
    # Col block left
    cbl = slice(0, psc1)
    # Col block right
    cbr = slice(nc2 + psc1, nc2 + psc2)

    # Mode normal
    trim1 = (rb1, cb)
    pcol1 = (rb1, cbl)
    ocol1 = (rb1, cbr)
    orow1 = (rb1m, cb)
    print trim1, ocol1, orow1, pcol1

    trim2 = (rb2, cb)
    pcol2 = (rb2, cbr)
    ocol2 = (rb2, cbl)
    orow2 = (rb2m, cb)
    print trim2, ocol2, orow2, pcol2

    finaldata = np.zeros(fshape, dtype='float32')

    finaldata[trim1] = direcfun(np.atleast_2d(np.arange(0, nc2)))
    finaldata[trim2] = direcfun(np.atleast_2d(np.arange(0, nc2)))

    finaldata[orow1] = 3
    finaldata[orow2] = 4

    finaldata[pcol1] = 5
    finaldata[pcol2] = 6

    finaldata[ocol1] = 7
    finaldata[ocol2] = 8

    hdu = fits.PrimaryHDU(data=finaldata)
    hdu.writeto(image, clobber=True)


def trim_and_o(image, out='trimmed.fits', direction='normal', bins='11'):
    '''Trim a MEGARA image with overscan.'''

    with fits.open(image) as hdul:
        hdu = trim_and_o_hdu(hdul[0])
        hdu.writeto(out, clobber=True)


def trim_and_o_hdu(hdu):
    '''Trim a MEGARA HDU with overscan.'''

    # FIXME: this should come from the header
    direction = 'normal'
    bins = '11'

    finaldata = trim_and_o_array(hdu.data, direction=direction, bins=bins)

    hdu.data = finaldata
    return hdu


def trim_and_o_array(array, direction='normal', bins='11'):
    '''Trim a MEGARA array with overscan.'''

    if direction not in _direc:
        raise ValueError("%s must be either 'normal' or 'mirror'" % direction)

    if direction == 'normal':
        direcfun = lambda x: x
    else:
        direcfun = np.fliplr

    if bins not in _binning:
        raise ValueError("%s must be one if '11', '12', '21, '22'" % bins)

    OSCANW = 100
    PSCANW = 50
    H_X_DIM = 2048
    H_Y_DIM = 2056

    bng = _binning[bins]

    nr2 = H_Y_DIM * 2 / bng[0]
    nc2 = H_X_DIM * 2 / bng[1]

    nr = H_Y_DIM / bng[0]
    nc = H_X_DIM / bng[1]

    oscan2 = OSCANW / bng[0]
    psc1 = PSCANW / bng[0]

    finaldata = np.empty((nr2, nc2), dtype='float32')
    finaldata[:nr, :] = direcfun(array[:nr, psc1:nc2 + psc1])
    finaldata[nr:, :] = direcfun(array[nr + oscan2:, psc1:nc2 + psc1])
    return finaldata

from numina.flow.processing import TagOptionalCorrector, TagFits
import logging

_logger = logging.getLogger('numina.recipes.megara')


class OverscanCorrector(TagOptionalCorrector):

    '''A Node that corrects a frame from overscan.'''

    def __init__(self, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        # FIXME: these should come from the header
        bng = [1, 1]
        nr = 2056 / bng[0]
        nc = 2048 / bng[1]
        nr2 = 2 * nr
        nc2 = 2 * nc
        oscan1 = 50 / bng[0]
        oscan2 = oscan1 * 2
        psc1 = 50 / bng[0]
        psc2 = 2 * psc1
        fshape = (nr2 + oscan2, nc2 + psc2)
        # Row block 1
        rb1 = slice(0, nr)
        rb1m = slice(nr, nr + oscan1)
        # Row block 2
        rb2 = slice(nr + oscan2, nr2 + oscan2)
        rb2m = slice(nr + oscan1, nr + oscan2)
        # Col block
        cb = slice(psc1, nc2 + psc1)

        # Col block left
        cbl = slice(0, psc1)
        # Col block right
        cbr = slice(nc2 + psc1, nc2 + psc2)

        # Mode normal
        self.trim1 = (rb1, cb)
        self.pcol1 = (rb1, cbl)
        self.ocol1 = (rb1, cbr)
        self.orow1 = (rb1m, cb)

        self.trim2 = (rb2, cb)
        self.pcol2 = (rb2, cbr)
        self.ocol2 = (rb2, cbl)
        self.orow2 = (rb2m, cb)

        if tagger is None:
            tagger = TagFits('NUM-OVPE', 'Over scan/prescan')

        super(OverscanCorrector, self).__init__(datamodel=datamodel,
                                                tagger=tagger,
                                                mark=mark,
                                                dtype=dtype)

    def _run(self, img):
        data = img[0].data

        p1 = data[self.pcol1].mean()
        _logger.debug('prescan1 is %f', p1)
        or1 = data[self.orow1].mean()
        _logger.debug('row overscan1 is %f', or1)
        oc1 = data[self.ocol1].mean()
        _logger.debug('col overscan1 is %f', oc1)
        avg = (p1 + or1 + oc1) / 3.0
        _logger.debug('average scan1 is %f', avg)
        data[self.trim1] -= avg

        p2 = data[self.pcol2].mean()
        _logger.debug('prescan2 is %f', p2)
        or2 = data[self.orow2].mean()
        _logger.debug('row overscan2 is %f', or2)
        oc2 = data[self.ocol2].mean()
        _logger.debug('col overscan2 is %f', oc2)
        avg = (p2 + or2 + oc2) / 3.0
        _logger.debug('average scan2 is %f', avg)
        data[self.trim2] -= avg
        return img


class TrimImage(TagOptionalCorrector):

    '''A Node that trims images.'''

    def __init__(self, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-TRIM', 'Trimming')

        super(TrimImage, self).__init__(datamodel=datamodel,
                                        tagger=tagger,
                                        mark=mark,
                                        dtype=dtype)

    def _run(self, img):
        _logger.debug('trimming image %s', img)

        img[0] = trim_and_o_hdu(img[0])

        return img


class ApertureExtractor(TagOptionalCorrector):

    '''A Node that extracts apertures.'''

    def __init__(self, trace, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-MAE', 'MEGARA Aperture extractor')

        super(ApertureExtractor, self).__init__(datamodel=datamodel,
                                                tagger=tagger,
                                                mark=mark,
                                                dtype=dtype)
        self.trace = trace

    def _run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('extracting apertures in image %s', imgid)
        rss = apextract(img[0].data, self.trace)
        img[0].data = rss

        return img


class FiberFlatCorrector(TagOptionalCorrector):

    '''A Node that corrects from fiber flat.'''

    def __init__(self, fiberflat, datamodel=None, mark=True,
                 tagger=None, dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-MFF', 'MEGARA Fiber flat correction')

        super(FiberFlatCorrector, self).__init__(datamodel=datamodel,
                                                 tagger=tagger,
                                                 mark=mark,
                                                 dtype=dtype)

        if isinstance(fiberflat, fits.HDUList):
            self.corr = fiberflat[0].data
        elif isinstance(fiberflat, np.ndarray):
            self.corr = fiberflat

        self.corrid = self.get_imgid(fiberflat)

    def _run(self, img):
        imgid = self.get_imgid(img)
        _logger.debug('correct from fiber flat in image %s', imgid)

        return img


def apextract(data, trace):
    '''Extract apertures.'''
    rss = np.empty((trace.shape[0], data.shape[1]), dtype='float32')
    for idx, r in enumerate(trace):
        l = r[0]
        r = r[2] + 1
        sl = (slice(l, r), )
        m = data[sl].sum(axis=0)
        rss[idx] = m
    return rss


def peakdet(v, delta, x=None, back=0.0):
    '''Basic peak detection.'''

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)

    if len(v) != len(x):
        raise ValueError('Input vectors v and x must have same length')

    if not np.isscalar(delta):
        raise TypeError('Input argument delta must be a scalar')

    if delta <= 0:
        raise ValueError('Input argument delta must be positive')

    maxtab = []
    mintab = []

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i, this in enumerate(v):
        if this > mx and this > back:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx - delta and this > back:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn + delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab), np.array(mintab)
