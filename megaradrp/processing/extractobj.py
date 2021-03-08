
# Copyright 2011-2021 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Extract objects from RSS images"""

import math
import uuid

import numpy
import astropy.wcs
import astropy.io.fits as fits
import astropy.units as u
from scipy.spatial import KDTree
from scipy.ndimage.filters import gaussian_filter
from numina.array.wavecalib.crosscorrelation import periodic_corr1d

#import megaradrp.datamodel as dm
import megaradrp.instrument.focalplane as fp
from megaradrp.processing.fluxcalib import update_flux_limits


def coverage_det(arr):
    """Compute coverage"""

    ma = numpy.max(arr)
    mn = numpy.min(arr)

    maxcover = len(arr) + 1

    if ma == 0:
        return slice(0, 0)

    if mn == 1:
        return slice(0, maxcover)

    # We have 0s and 1s
    diffa = numpy.diff(arr)

    # For the left, find first +1
    val_f, = numpy.where(diffa==1)
    pix_f = 0
    if len(val_f) > 0:
        pix_f = val_f[0] + 1

    # For the rigth, find first -1
    val_l, = numpy.where(diffa==-1)
    pix_l = maxcover
    if len(val_l) > 0:
        pix_l = val_l[0] + 1

    return slice(pix_f, pix_l)


def extract_star(rssimage, position, npoints, fiberconf, logger=None):
    """Extract a star given its center and the number of fibers to extract"""

    # FIXME: handle several positions

    logger.info('extracting star')

    # fiberconf = dm.get_fiberconf(rssimage)
    logger.debug("Configuration UUID is %s", fiberconf.conf_id)

    rssdata = rssimage[0].data
    pdata = rssimage['wlmap'].data

    points = [position]
    fibers = fiberconf.connected_fibers(valid_only=True)
    grid_coords = []
    for fiber in fibers:
        grid_coords.append((fiber.x, fiber.y))
    # setup kdtree for searching
    kdtree = KDTree(grid_coords)

    # Other posibility is
    # query using radius instead
    # radius = 1.2
    # kdtree.query_ball_point(points, k=7, r=radius)

    dis_p, idx_p = kdtree.query(points, k=npoints)

    logger.info('Using %d nearest fibers', npoints)
    totals = []
    for diss, idxs, point in zip(dis_p, idx_p, points):
        # For each point
        logger.info('For point %s', point)
        colids = []
        coords = []
        for dis, idx in zip(diss, idxs):
            fiber = fibers[idx]
            colids.append(fiber.fibid - 1)
            logger.debug('adding fibid %d', fiber.fibid)

            coords.append((fiber.x, fiber.y))

        colids.sort()
        flux_fiber = rssdata[colids]
        flux_total = rssdata[colids].sum(axis=0)
        coverage_total = pdata[colids].sum(axis=0)

        max_cover = coverage_total.max()
        some_value_region = coverage_total > 0
        max_value_region = coverage_total == max_cover
        valid_region = max_value_region

        # Interval with maximum coverage
        #nz_max, = numpy.nonzero(numpy.diff(max_value_region))
        # Interval with at least 1 fiber
        #nz_some, = numpy.nonzero(numpy.diff(some_value_region))
        nz_max_slice = coverage_det(max_value_region.astype('int'))
        nz_some_slice = coverage_det(some_value_region.astype('int'))

        # Collapse the flux in the optimal region
        perf = flux_fiber[:, nz_max_slice].sum(axis=1)
        # Contribution of each fiber to the total flux, 1D
        perf_norm = perf / perf.sum()
        contributions = numpy.zeros(shape=(rssdata.shape[0],))
        contributions[colids] = perf_norm
        # Contribution of each fiber to the total flux, 2D
        flux_per_fiber = pdata * contributions[:, numpy.newaxis]
        flux_sum = flux_per_fiber.sum(axis=0)
        # In the region max_value_region, flux_sum == 1
        # In some_value_region 0 < flux_sum < 1
        # Outside is flux_sum == 0
        flux_correction = numpy.zeros_like(flux_sum)
        flux_correction[some_value_region] = 1.0 / flux_sum[some_value_region]

        # Limit to 10
        flux_correction = numpy.clip(flux_correction, 0, 10)

        flux_total_c = flux_total * flux_correction
        # plt.axvspan(nz_some[0], nz_some[1], alpha=0.2, facecolor='r')
        # plt.axvspan(nz_max[0], nz_max[1], alpha=0.2, facecolor='r')
        # plt.plot(flux_correction)
        # plt.show()
        #
        # plt.axvspan(nz_some[0], nz_some[1], alpha=0.2)
        # plt.axvspan(nz_max[0], nz_max[1], alpha=0.2)
        # plt.plot(flux_total)
        # plt.plot(flux_total * flux_correction, 'g')
        # ax2 = plt.gca().twinx()
        # ax2.plot(coverage_total)
        #
        # plt.show()
        pack = flux_total_c, colids, nz_max_slice, nz_some_slice
        # pack = flux_total_c
        totals.append(pack)

    return totals[0]


def compute_centroid(rssdata, fiberconf, c1, c2, point, logger=None):
    """Compute centroid near coordinates given by 'point'"""

    logger.debug("LCB configuration is %s", fiberconf.conf_id)

    fibers = fiberconf.connected_fibers(valid_only=True)
    grid_coords = []
    for fiber in fibers:
        grid_coords.append((fiber.x, fiber.y))
    # setup kdtree for searching
    kdtree = KDTree(grid_coords)

    # Other posibility is
    # query using radius instead
    # radius = 1.2
    # kdtree.query_ball_point(points, k=7, r=radius)

    npoints = 19
    # 1 + 6  for first ring
    # 1 + 6  + 12  for second ring
    # 1 + 6  + 12  + 18 for third ring
    points = [point]
    dis_p, idx_p = kdtree.query(points, k=npoints)

    logger.info('Using %d nearest fibers', npoints)
    for diss, idxs, point in zip(dis_p, idx_p, points):
        # For each point
        logger.info('For point %s', point)
        colids = []
        coords = []
        for dis, idx in zip(diss, idxs):
            fiber = fibers[idx]
            colids.append(fiber.fibid - 1)
            coords.append((fiber.x, fiber.y))

        coords = numpy.asarray(coords)
        flux_per_cell = rssdata[colids, c1:c2].mean(axis=1)
        flux_per_cell_total = flux_per_cell.sum()
        flux_per_cell_norm = flux_per_cell / flux_per_cell_total
        # centroid
        scf = coords.T * flux_per_cell_norm
        centroid = scf.sum(axis=1)
        logger.info('centroid: %s', centroid)
        # central coords
        c_coords = coords - centroid
        scf0 = scf - centroid[:, numpy.newaxis] * flux_per_cell_norm
        mc2 = numpy.dot(scf0, c_coords)
        logger.info('2nd order moments, x2=%f, y2=%f, xy=%f', mc2[0,0], mc2[1,1], mc2[0,1])
        return centroid


def compute_dar(img, logger=None, debug_plot=False):
    """Compute Diferencial Atmospheric Refraction"""

    fp_conf = fp.FocalPlaneConf.from_img(img)
    wlcalib = astropy.wcs.WCS(img[0].header)

    rssdata = img[0].data
    cut1 = 500
    cut2 = 3500
    colids = []
    x = []
    y = []
    for fiber in fp_conf.fibers.values():
        colids.append(fiber.fibid - 1)
        x.append(fiber.x)
        y.append(fiber.y)

    cols = []
    xdar = []
    ydar = []
    delt = 50

    point = [2.0, 2.0]
    # Start in center of range
    ccenter = (cut2 + cut1) // 2
    # Left
    for c in range(ccenter, cut1, -delt):
        c1 = c - delt // 2
        c2 = c + delt // 2

        z = rssdata[colids, c1:c2].mean(axis=1)
        centroid = compute_centroid(rssdata, fp_conf, c1, c2, point, logger=logger)
        cols.append(c)
        xdar.append(centroid[0])
        ydar.append(centroid[1])
        point = centroid

    cols.reverse()
    xdar.reverse()
    ydar.reverse()

    point = [2.0, 2.0]
    # Star over
    # Right
    for c in range(ccenter, cut2, delt):
        c1 = c - delt // 2
        c2 = c + delt // 2
        z = rssdata[colids, c1:c2].mean(axis=1)
        centroid = compute_centroid(rssdata, fp_conf, c1, c2, point)
        cols.append(c)
        xdar.append(centroid[0])
        ydar.append(centroid[1])
        point = centroid

    rr = [[col, 0] for col in cols]
    world = wlcalib.wcs_pix2world(rr, 0)

    if debug_plot:
        import matplotlib.pyplot as plt
        import megaradrp.visualization as vis
        import megaradrp.simulation.refraction as r
        from astropy import units as u

        plt.subplots_adjust(hspace=0.5)
        plt.subplot(111)
        ax = plt.gca()
        plt.xlim([-8, 8])
        plt.ylim([-8, 8])
        col = vis.hexplot(ax, x, y, z, cmap=plt.cm.YlOrRd_r)
        plt.title(f"Fiber map, {c1} {c2}")
        cb = plt.colorbar(col)
        cb.set_label('counts')
        plt.show()

        zenith_distance = 60 * u.deg
        press = 79993.2 * u.Pa
        rel = 0.013333333
        temp = 11.5 * u.deg_C

        ll = r.differential_p(
                zenith_distance,
                wl=world[:,0] * u.AA,
                wl_reference=4025 * u.AA,
                temperature=temp,
                pressure=press,
                relative_humidity=rel,
        )
        plt.plot(world[:,0], xdar, '*-')
        plt.plot(world[:,0], ydar, '*-')
        plt.plot(world[:, 0], 2.0 * u.arcsec + ll.to(u.arcsec), '-')
        plt.show()

    return world[:, 0], xdar, ydar


def mix_values(wcsl, spectrum, star_interp):

    r1 = numpy.arange(spectrum.shape[0])
    r2 = r1 * 0.0
    lm = numpy.array([r1, r2])
    # Values are 0-based
    wavelen_ = wcsl.all_pix2world(lm.T, 0)
    if wcsl.wcs.cunit[0] == u.dimensionless_unscaled:
        # CUNIT is empty, assume Angstroms
        wavelen = wavelen_[:, 0] * u.AA
    else:
        wavelen = wavelen_[:, 0] * wcsl.wcs.cunit[0]

    wavelen_aa = wavelen.to(u.AA)

    response_0 = spectrum
    mag_ref = star_interp(wavelen_aa) * u.ABmag
    response_1 = mag_ref.to(u.Jy).value

    return wavelen_aa, response_0, response_1


def compute_broadening(flux_low, flux_high, sigmalist,
                       remove_mean=False, frac_cosbell=None, zero_padding=None,
                       fminmax=None, naround_zero=None, nfit_peak=None):

    # normalize each spectrum dividing by its median
    flux_low /= numpy.median(flux_low)
    flux_high /= numpy.median(flux_high)

    offsets = []
    fpeaks = []
    sigmalist = numpy.asarray(sigmalist)
    for sigma in sigmalist:
        # broaden reference spectrum
        flux_ref_broad = gaussian_filter(flux_high, sigma)
        # plot the two spectra

        # periodic correlation between the two spectra
        offset, fpeak = periodic_corr1d(
            flux_ref_broad, flux_low,
            remove_mean=remove_mean,
            frac_cosbell=frac_cosbell,
            zero_padding=zero_padding,
            fminmax=fminmax,
            naround_zero=naround_zero,
            nfit_peak=nfit_peak,
            norm_spectra=True,
        )
        offsets.append(offset)
        fpeaks.append(fpeak)

    fpeaks = numpy.asarray(fpeaks)
    offsets = numpy.asarray(offsets)

    # import matplotlib.pyplot as plt
    # #
    # plt.plot(sigmalist, offsets, color='r')
    # ax2 = plt.gca().twinx()
    # ax2.plot(sigmalist, fpeaks, color='b')
    # plt.show()
    #
    offset_broad = offsets[numpy.argmax(fpeaks)]
    sigma_broad = sigmalist[numpy.argmax(fpeaks)]

    return offset_broad, sigma_broad


def generate_sensitivity(final, spectrum, star_interp, extinc_interp,
                         wl_coverage1, wl_coverage2, sigma=20.0):
    """Generate sensitivity response"""

    wcsl = astropy.wcs.WCS(final[0].header)

    r1 = numpy.arange(final[0].shape[1])
    r2 = r1 * 0.0
    lm = numpy.array([r1, r2])
    # Values are 0-based
    wavelen_ = wcsl.all_pix2world(lm.T, 0)
    if wcsl.wcs.cunit[0] == u.dimensionless_unscaled:
        # CUNIT is empty, assume Angstroms
        wavelen = wavelen_[:, 0] * u.AA
    else:
        wavelen = wavelen_[:, 0] * wcsl.wcs.cunit[0]

    wavelen_aa = wavelen.to(u.AA)

    airmass = final[0].header['AIRMASS']
    exptime = final[0].header['EXPTIME']

    response_0 = spectrum / exptime
    r0max = response_0.max()
    valid = response_0 > 0
    if r0max <= 0:
        raise ValueError("maximum of 'spectrum' is <= 0")
    if sum(valid) <= 0:
        raise ValueError("'spectrum' is <= 0")

    # In magAB
    # f(Jy) = 3631 * 10^-0.4 mAB
    mag_ref = (star_interp(wavelen_aa) + extinc_interp(wavelen_aa) * airmass) * u.ABmag
    response_1 = mag_ref.to(u.Jy).value

    r1max = response_1.max()
    if r1max <= 0:
        raise ValueError("maximum of 'star_interp' is <= 0")

    r0 = response_0 / r0max
    r1 = response_1 / r1max

    pixm1 = wl_coverage1.start
    pixm2 = wl_coverage1.stop
    pixr1 = wl_coverage2.start
    pixr2 = wl_coverage2.stop

    pixlims = {}
    pixlims['PIXLIMR1'] = pixr1 + 1  # Convert to 1-ref
    pixlims['PIXLIMR2'] = pixr2
    pixlims['PIXLIMM1'] = pixm1 + 1  # Convert to 1-ref
    pixlims['PIXLIMM2'] = pixm2

    max_valid = numpy.zeros_like(valid)
    max_valid[wl_coverage1] = True

    partial_valid = numpy.zeros_like(valid)
    partial_valid[wl_coverage2] = True

    valid = numpy.zeros_like(response_0)
    valid[wl_coverage1] = 1

    pixf1, pixf2 = int(math.floor(pixm1 +  2* sigma)), int(math.ceil(pixm2 - 2 * sigma))

    pixlims['PIXLIMF1'] = pixf1 + 1
    pixlims['PIXLIMF2'] = pixf2

    flux_valid = numpy.zeros_like(valid, dtype='bool')
    flux_valid[pixf1:pixf2] = True

    if sigma > 0:
        r0_ens = gaussian_filter(r0, sigma=sigma)
    else:
        r0_ens = r0

    ratio2 = r0_ens / r1
    s_response = ratio2 * (r0max / r1max)

    # FIXME: add history
    sens = fits.PrimaryHDU(s_response, header=final[0].header)
    # delete second axis keywords
    # FIXME: delete axis with wcslib
    for key in ['CRPIX2', 'CRVAL2', 'CDELT2', 'CTYPE2']:
        if key in sens.header:
            del sens.header[key]

    sens.header['uuid'] = str(uuid.uuid1())
    sens.header['tunit'] = ('Jy', "Final units")

    update_flux_limits(sens.header, pixlims, wcs=wcsl, ref=0)

    return sens


