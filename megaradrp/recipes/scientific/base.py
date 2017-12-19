#
# Copyright 2011-2017 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Base scientific recipe for MEGARA"""


import uuid
import math

from astropy.io import fits
import numpy as np
import numpy
import matplotlib.pyplot as plt

from numina.core import Parameter
from numina.core.requirements import ObservationResultRequirement
from numina.flow import SerialFlow
from numina.array import combine

from megaradrp.core.recipe import MegaraBaseRecipe
import megaradrp.requirements as reqs
from megaradrp.processing.combine import basic_processing_with_combination
from megaradrp.utils import copy_img
from megaradrp.processing.aperture import ApertureExtractor
from megaradrp.processing.wavecalibration import WavelengthCalibrator
from megaradrp.processing.fiberflat import Splitter, FlipLR, FiberFlatCorrector
from megaradrp.processing.twilight import TwilightCorrector


class ImageRecipe(MegaraBaseRecipe):
    """Base Image."""

    # Requirements  
    obresult = ObservationResultRequirement()
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()
    master_slitflat = reqs.MasterSlitFlatRequirement()
    master_wlcalib = reqs.WavelengthCalibrationRequirement()
    master_fiberflat = reqs.MasterFiberFlatRequirement()
    master_twilight = reqs.MasterTwilightRequirement()
    master_traces = reqs.MasterAperturesRequirement()
    extraction_offset = Parameter([0.0], 'Offset traces for extraction')
    ignored_sky_bundles = Parameter([], 'Ignore these sky bundles')
    master_sensitivity = reqs.SensitivityRequirement()
    reference_extinction = reqs.ReferenceExtinction()
    relative_threshold = Parameter(0.3, 'Threshold for peak detection')

    def base_run(self, rinput):

        # 2D reduction
        flow1 = self.init_filters(rinput, rinput.obresult.configuration)
        img = basic_processing_with_combination(rinput, flow1, method=combine.median)
        hdr = img[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(img, 'reduced_image.fits')

        reduced2d = copy_img(img)

        # 1D, extraction, Wl calibration, Flat fielding
        reduced_rss = self.run_reduction_1d(img,
            rinput.master_traces, rinput.master_wlcalib,
            rinput.master_fiberflat, rinput.master_twilight,
            offset=rinput.extraction_offset
        )
        self.save_intermediate_img(reduced_rss, 'reduced_rss.fits')

        return reduced2d, reduced_rss

    def run_reduction_1d(self, img, tracemap, wlcalib, fiberflat, twflat=None, offset=None):
        # 1D, extraction, Wl calibration, Flat fielding
        correctors = []
        correctors.append(ApertureExtractor(tracemap, self.datamodel, offset=offset))
        correctors.append(FlipLR())
        correctors.append(WavelengthCalibrator(wlcalib, self.datamodel))
        correctors.append(FiberFlatCorrector(fiberflat.open(), self.datamodel))

        if twflat:
            correctors.append(TwilightCorrector(twflat.open(), self.datamodel))

        flow2 = SerialFlow(correctors)

        reduced_rss =  flow2(img)
        return reduced_rss

    def run_sky_subtraction(self, img, ignored_sky_bundles):
        # Sky subtraction
        self.logger.info('obtain fiber information')
        sky_img = copy_img(img)
        final_img = copy_img(img)
        fiberconf = self.datamodel.get_fiberconf(sky_img)
        # Sky fibers
        skyfibs = fiberconf.sky_fibers(valid_only=True,
                                       ignored_bundles=ignored_sky_bundles)
        self.logger.debug('sky fibers are: %s', )
        # Create empty sky_data
        target_data = img[0].data

        target_map = img['WLMAP'].data
        sky_data = numpy.zeros_like(img[0].data)
        sky_map = numpy.zeros_like(img['WLMAP'].data)
        sky_img[0].data = sky_data

        for fibid in skyfibs:
            rowid = fibid - 1
            sky_data[rowid] = target_data[rowid]
            sky_map[rowid] = target_map[rowid]
            if False:

                plt.plot(sky_data[rowid])
                plt.title("%d" % fibid)
                plt.show()
        # Sum
        coldata = sky_data.sum(axis=0)
        colsum = sky_map.sum(axis=0)

        # Divide only where map is > 0
        mask = colsum > 0
        avg_sky = numpy.zeros_like(coldata)
        avg_sky[mask] = coldata[mask] / colsum[mask]

        # This should be done only on valid fibers
        # The information of which fiber is valid
        # is in the tracemap, not in the header
        for fibid in fiberconf.valid_fibers():
            rowid = fibid - 1
            final_img[0].data[rowid, mask] = img[0].data[rowid, mask] - avg_sky[mask]

        return final_img, img, sky_img

    def read_wcs(self, hdr):
        crpix = hdr['CRPIX1']
        wlr0 = hdr['CRVAL1']
        delt = hdr['CDELT1']
        return crpix, wlr0, delt

    def get_wcallib(self, lambda1, lambda2, fibras, traces, rss, neigh_info, grid):

        # Take a look at == []
        indices = []
        wlcalib = []
        for elem in traces:
            if elem['aperture']['function']['coefficients']:
                wlcalib.append(elem['aperture']['function']['coefficients'])
                if len(indices)==0:
                    indices.append(0)
                else:
                    indices.append(indices[-1])
            else:
                indices.append(indices[-1]+1)

        wlcalib_aux = np.asarray(wlcalib)
        final, wcsdata = self.resample_rss_flux(rss, wlcalib_aux, indices)

        hdu_f = fits.PrimaryHDU(final)
        hdu_f.writeto('resample_rss.fits', clobber=True)

        fibras.sort()

        suma = np.sum(final[fibras.astype(int),lambda1:lambda2],axis=1)
        sky = np.min(suma)
        sumaparcial = suma - sky

        neigh_info = neigh_info[np.argsort(neigh_info[:, 0])]

        centroid_x = np.multiply(sumaparcial,neigh_info[:,2])
        centroid_x = np.sum(centroid_x, axis=0)

        centroid_y = np.multiply(sumaparcial,neigh_info[:,3])
        centroid_y = np.sum(centroid_y, axis=0)

        sumatotal = np.sum(sumaparcial, axis=0)
        self.logger.info( "total sum: %s", sumatotal)

        second_order = []
        aux = np.sum(np.multiply(suma,(neigh_info[:,2] - np.mean(neigh_info[:,2]))**2),axis=0)
        second_order.append(np.divide(aux ,np.sum(suma, axis=0)))
        self.logger.info("Second order momentum X: %s", second_order[0])

        aux = np.sum(np.multiply(suma,(neigh_info[:,3] - np.mean(neigh_info[:,3]))**2),axis=0)
        second_order.append(np.divide(aux ,np.sum(suma, axis=0)))
        self.logger.info("Second order momentum Y: %s", second_order[1])

        aux = np.multiply(neigh_info[:,3] - np.mean(neigh_info[:,3]),neigh_info[:,2] - np.mean(neigh_info[:,2]))
        aux = np.sum(np.multiply(aux,suma))
        cov = np.divide(aux ,np.sum(suma, axis=0))
        self.logger.info("Cov X,Y: %s", cov)

        centroid_x = np.divide(centroid_x, sumatotal)
        self.logger.info( "centroid_x: %s", centroid_x)

        centroid_y = np.divide(centroid_y, sumatotal)
        self.logger.info("centroid_y: %s", centroid_y)

        centroid = [centroid_x, centroid_y]

        peak = np.sum(final[grid.get_fiber(centroid),lambda1:lambda2],axis=0)

        return centroid, sky, peak, second_order, cov

    def generate_solution(self, points, centroid, sky, fiber, peaks, second_order, cova):
        result = []
        for cont, value in enumerate(points):
            lista = (value[0], value[1], centroid[cont][0],centroid[cont][1], sky[cont], fiber[cont], peaks[cont], second_order[cont][0], second_order[cont][1], cova[cont])
            result.append(lista)
        return np.array(result, dtype=[('x_point','float'),('y_point','float'),('x_centroid','float'),('y_centroid','float'), ('sky','float'),('fiber','int'),('peak','float'),('x_second_order','float'), ('y_second_order','float'), ('covariance','float') ])

    def generateJSON(self, points, centroid, sky, fiber, peaks, second_order, cova):
        '''
        '''

        self.logger.info('start JSON generation')

        result = []
        for cont, value in enumerate(points):
            obj = {
                'points': value,
                'centroid': centroid[cont],
                'sky':sky[cont],
                'fiber': fiber[cont],
                'peak': peaks[cont],
                'second_order': second_order[cont],
                'covariance': cova[cont]
            }
            result.append(obj)

        self.logger.info('end JSON generation')

        return result

    def compute_dar(self, img):
        import numpy.polynomial.polynomial as pol
        import astropy.wcs

        fiberconf = self.datamodel.get_fiberconf(img)
        wlcalib = astropy.wcs.WCS(img[0].header)

        rssdata = img[0].data
        cut1 = 500
        cut2 = 3500
        colids = []
        x = []
        y = []
        for fiber in fiberconf.fibers.values():
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
            centroid = self.centroid(rssdata, fiberconf, c1, c2, point)
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
            centroid = self.centroid(rssdata, fiberconf, c1, c2, point)
            cols.append(c)
            xdar.append(centroid[0])
            ydar.append(centroid[1])
            point = centroid

        rr = [[col, 0] for col in cols]
        world = wlcalib.wcs_pix2world(rr, 0)

        if False:
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
            plt.title("Fiber map, %s %s" % (c1, c2))
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

        # fit something

        print('DAR, x:', pol.polyfit(world[:, 0], xdar, deg=3))
        print('DAR: y:', pol.polyfit(world[:, 0], ydar, deg=3))

    def centroid(self, rssdata, fiberconf, c1, c2, point):
        from scipy.spatial import KDTree

        self.logger.debug("LCB configuration is %s", fiberconf.conf_id)

        fibers = fiberconf.conected_fibers(valid_only=True)
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

        self.logger.info('Using %d nearest fibers', npoints)
        for diss, idxs, point in zip(dis_p, idx_p, points):
            # For each point
            self.logger.info('For point %s', point)
            colids = []
            coords = []
            for dis, idx in zip(diss, idxs):
                fiber = fibers[idx]
                colids.append(fiber.fibid - 1)
                coords.append((fiber.x, fiber.y))

            coords = np.asarray(coords)
            flux_per_cell = rssdata[colids, c1:c2].mean(axis=1)
            flux_per_cell_total = flux_per_cell.sum()
            flux_per_cell_norm = flux_per_cell / flux_per_cell_total
            # centroid
            scf = coords.T * flux_per_cell_norm
            centroid = scf.sum(axis=1)
            self.logger.info('centroid: %s', centroid)
            # central coords
            c_coords = coords - centroid
            scf0 = scf - centroid[:, np.newaxis] * flux_per_cell_norm
            mc2 = np.dot(scf0, c_coords)
            self.logger.info('2nd order moments, x2=%f, y2=%f, xy=%f', mc2[0,0], mc2[1,1], mc2[0,1])
            return centroid

    def extract_stars(self, final, position, npoints):
        from scipy.spatial import KDTree

        import matplotlib.pyplot as plt

        self.logger.info('extracting star')

        fiberconf = self.datamodel.get_fiberconf(final)
        self.logger.debug("Configuration UUID is %s", fiberconf.conf_id)
        rssdata = final[0].data
        pdata = final['wlmap'].data

        points = [position]
        fibers = fiberconf.conected_fibers(valid_only=True)
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

        self.logger.info('Using %d nearest fibers', npoints)
        totals = []
        for diss, idxs, point in zip(dis_p, idx_p, points):
            # For each point
            self.logger.info('For point %s', point)
            colids = []
            coords = []
            for dis, idx in zip(diss, idxs):
                fiber = fibers[idx]
                colids.append(fiber.fibid - 1)
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
            nz_max, = numpy.nonzero(numpy.diff(max_value_region))
            # Interval with at least 1 fiber
            nz_some, = numpy.nonzero(numpy.diff(some_value_region))

            # Collapse the flux in the optimal region
            perf = flux_fiber[:, nz_max[0]+ 1: nz_max[1] + 1].sum(axis=1)
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
            pack = flux_total_c, nz_max, nz_some
            # pack = flux_total_c
            totals.append(pack)

        return totals

    def generate_sensitivity(self, final, spectrum, star_interp, extinc_interp, cover1, cover2, sigma=20.0):

        from astropy.wcs import WCS
        import matplotlib.pyplot as plt
        from scipy.ndimage.filters import uniform_filter, gaussian_filter

        crpix, wlr0, delt = self.read_wcs(final[0].header)

        wcsl = WCS(final[0].header)

        r1 = numpy.arange(final[0].shape[1])
        r2 = r1 * 0.0
        lm = numpy.array([r1, r2])
        wavelen_ = wcsl.all_pix2world(lm.T, 0.0)
        wavelen = wavelen_[:, 0]

        airmass = final[0].header['AIRMASS']
        exptime = final[0].header['EXPTIME']

        response_0 = spectrum / exptime
        valid = response_0 > 0
        # In magAB
        # f(Jy) = 3631 * 10^-0.4 mAB

        response_1 = 3631 * numpy.power(10.0, -0.4 * (star_interp(wavelen) + extinc_interp(wavelen) * airmass))
        r0max = response_0.max()
        r1max = response_1.max()
        r0 = response_0 / r0max
        r1 = response_1 / r1max

        pixm1, pixm2 = cover1
        pixr1, pixr2 = cover2

        max_valid = numpy.zeros_like(valid)
        max_valid[pixm1:pixm2 + 1] = True

        partial_valid = numpy.zeros_like(valid)
        partial_valid[pixr1:pixr2 + 1] = True

        valid = numpy.ones_like(response_0)
        valid[pixm2:] = 0
        valid[:pixm1+1] = 0

        sampling = int(2.0 / delt)
        # print('sampling is', sampling)

        # with numpy.errstate(invalid='ignore', divide='ignore'):
        #    ratio = r0 / r1

        pixf1, pixf2 = int(math.floor(pixm1 +  2* sigma)), int(math.ceil(pixm2 - 2 * sigma))
        # calc1 = [[pixm1, pixm2, pixr1, pixr2, pixf1, pixf2], [0, 0, 0, 0, 0, 0]]

        flux_valid = numpy.zeros_like(valid, dtype='bool')
        flux_valid[pixf1:pixf2 + 1] = True

        # lm2 = numpy.array(calc1)
        # wavelen_ = wcsl.all_pix2world(lm2.T, 0.0)

        r0_ens = gaussian_filter(r0, sigma=sigma)
        # mf = uniform_filter(r0, size=resolution)
        #
        # plt.plot(wavelen, r0, 'b*-')
        # plt.plot(wavelen, r0_ens, 'r*-')
        # plt.axvspan(lm3[0], lm3[1], facecolor='g', alpha=0.2)
        # plt.axvspan(lm3[2], lm3[3], facecolor='r', alpha=0.2)
        # plt.axvspan(lm3[4], lm3[5], facecolor='b', alpha=0.2)
        # #plt.plot(wavelen, r1)
        # plt.show()

        ratio2 = r0_ens / r1
        s_response = ratio2 * (r0max / r1max)

        # FIXME: include history
        sens = fits.PrimaryHDU(s_response, header=final[0].header)
        sens.header['uuid'] = str(uuid.uuid1())
        sens.header['tunit'] = ('Jy', "Final units")
        sens.header['PIXLIMF1'] = (pixf1 + 1, "Start of valid flux calibration")
        sens.header['PIXLIMF2'] = (pixf2 + 1, "End of valid flux calibration")
        sens.header['PIXLIMR1'] = (pixr1 + 1, 'Start of region with at least one fiber')
        sens.header['PIXLIMR2'] = (pixr2 + 1, 'End of region with at least one fiber')
        sens.header['PIXLIMM1'] = (pixm1 + 1, 'Start of region with all fibers')
        sens.header['PIXLIMM2'] = (pixm2 + 1, 'End of region with all fibers')
        return sens
