#
# Copyright 2011-2022 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""LCB Standard Star Image Recipe for Megara"""


from astropy import constants as const
import astropy.io.fits as fits
import astropy.units as u
import astropy.wcs
from scipy.interpolate import interp1d

from numina.types.datatype import PlainPythonType
from numina.types.datatype import ListOfType
from numina.types.multitype import MultiType
from numina.array.numsplines import AdaptiveLSQUnivariateSpline
from numina.core import Result, Parameter
from numina.core.requirements import Requirement
from numina.core.validator import range_validator
from numina.exceptions import RecipeError
from numina.types.array import ArrayType

from megaradrp.instrument.focalplane import FocalPlaneConf
from megaradrp.ntypes import Point2D
from megaradrp.ntypes import ProcessedRSS, ProcessedFrame, ProcessedSpectrum
from megaradrp.ntypes import ReferenceSpectrumTable, ReferenceExtinctionTable
from megaradrp.ntypes import MasterSensitivity
from megaradrp.processing.extractobj import extract_star, generate_sensitivity
from megaradrp.processing.extractobj import mix_values, compute_broadening
from megaradrp.processing.centroid import calc_centroid_brightest
from megaradrp.recipes.scientific.base import ImageRecipe


class LCBStandardRecipe(ImageRecipe):
    """Process LCB Standard Star Recipe.

    This recipe processes a set of images
    obtained in **LCB Stardard Star image** mode and returns
    the total flux of the star.

    See Also
    --------
    megaradrp.recipes.calibration.mosstdstar.MOSStandardRecipe

    Notes
    -----
    Images provided by `obresult` are trimmed and corrected
    from overscan, bad pixel mask (if `master_bpm` is not None),
    bias, dark current (if `master_dark` is not None) and
    slit-flat (if `master_slitflat` is not None).

    Images thus corrected are then stacked using the median.
    The result of the combination is saved as an intermediate result, named
    'reduced_image.fits'. This combined image is also returned in the field
    `reduced_image` of the recipe result.

    The apertures in the 2D image are extracted, using the information in
    `master_apertures` and resampled according to the wavelength calibration in
    `master_wlcalib`. Then is divided by the `master_fiberflat`.
    The resulting RSS is saved as an intermediate
    result named 'reduced_rss.fits'. This RSS is also returned in the field
    `reduced_rss` of the recipe result.

    The sky is subtracted by combining the the fibers marked as `SKY`
    in the fibers configuration. The RSS with sky subtracted is returned ini the
    field `final_rss` of the recipe result.

    The flux of the star is computed by adding summing the fibers in `nrings` around
    the central spaxel containing the star and returned as `star_spectrum`.

    """
    position = Requirement(Point2D, "Position of the reference object", optional=True)
    nrings = Parameter(3, "Number of rings to extract the star",
                       validator=range_validator(minval=1))
    reference_spectrum = Requirement(ReferenceSpectrumTable, "Spectrum of reference star")
    reference_spectrum_velocity = Parameter(0.0, 'Radial velocity (km/s) of reference spectrum')
    reference_extinction = Requirement(ReferenceExtinctionTable, "Reference extinction")
    degrade_resolution_target = Parameter('object', 'Spectrum with higher resolution',
                                          choices=['object']
                                          )
    # TODO: Implement the possibility of the reference having higher resolution
    # degrade_resolution_target = Parameter('object', 'Spectrum with higher resolution',
    #                                       choices=['object', 'reference']
    #                                      )
    degrade_resolution_method = Parameter('fixed', 'Method to degrade the resolution',
                                          choices=['none', 'fixed', 'auto']
                                          )
    sigma_resolution = Parameter(20.0, 'sigma Gaussian filter to degrade resolution ')
    smoothing_knots = Requirement(
        MultiType(
            PlainPythonType(ref=3, validator=range_validator(minval=3)),
            ListOfType(PlainPythonType(ref=0.0), nmin=3)
        ),
        description='List of nodes or number of nodes for sensitivity smoothing',
        default=3,
        optional=True
    )

    reduced_image = Result(ProcessedFrame)
    final_rss = Result(ProcessedRSS)
    reduced_rss = Result(ProcessedRSS)
    sky_rss = Result(ProcessedRSS)
    star_spectrum = Result(ProcessedSpectrum)
    master_sensitivity = Result(MasterSensitivity)
    sensitivity_raw = Result(ProcessedSpectrum)
    fiber_ids = Result(ArrayType(fmt='%d'))
    sigma = Result(float)

    def set_base_headers(self, hdr):
        """Set metadata in FITS headers."""
        hdr = super(LCBStandardRecipe, self).set_base_headers(hdr)
        hdr['NUMTYPE'] = ('MASTER_SENSITIVITY', 'Product type')
        hdr['IMGTYPE'] = ('MASTER_SENSITIVITY', 'Product type')
        return hdr

    def run(self, rinput):

        self.logger.info('starting LCBStandardRecipe reduction')

        # Try to guard against receiving here something
        # that is not in magAB
        # TODO: implement this in ReferenceSpectrumTable
        maxm = max(rinput.reference_spectrum[:, 1])
        if maxm > 100:
            # If the column here has values greater than 100
            # this could not be a magnitude
            raise RecipeError("the maximum flux of 'reference_spectrum' is > 100, "
                              "check the flux unit (it has to be magAB)")

        # Create InstrumentModel
        # ins1 = rinput.obresult.configuration
        #
        reduced2d, rss_data = super(LCBStandardRecipe, self).base_run(rinput)
        # tags = rinput.obresult.tags
        ins2 = rinput.obresult.profile
        ins2.configure_with_image(rss_data)
        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(
            rss_data,
            sky_rss=rinput.sky_rss,
            ignored_sky_bundles=rinput.ignored_sky_bundles
        )
        self.logger.info('end sky subtraction')

        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring
        # 1 + 6 * Sum_i=0^n =  1 + 3 * n * (n +1)
        # Using three rings around central point

        # If position is None, find the brightest spaxel
        # and use the centroid
        if rinput.position is None:
            self.logger.info('finding centroid of brightest spaxel')
            extraction_region = [1000, 3000]
            nrings = rinput.nrings
            position = calc_centroid_brightest(final, extraction_region, nrings)
        else:
            position = rinput.position
        self.logger.info('central position is %s', position)

        self.logger.debug('adding %d nrings', rinput.nrings)
        npoints = 1 + 3 * rinput.nrings * (rinput.nrings + 1)
        self.logger.debug('adding %d fibers', npoints)

        fp_conf = FocalPlaneConf.from_img(final)
        spectra_pack = extract_star(final, position, npoints,
                                    fp_conf, logger=self.logger)

        spectrum, colids, wl_cover1, wl_cover2 = spectra_pack
        star_spectrum = fits.PrimaryHDU(spectrum, header=final[0].header)

        rad_vel = rinput.reference_spectrum_velocity * u.km / u.s
        factor = 1 + rad_vel / const.c

        star_interp = interp1d(rinput.reference_spectrum[:, 0] / factor,
                               rinput.reference_spectrum[:, 1])

        extinc_interp = interp1d(rinput.reference_extinction[:, 0],
                                 rinput.reference_extinction[:, 1])

        fiber_ids = [colid + 1 for colid in colids]

        wcsl = astropy.wcs.WCS(final[0].header)
        wl_aa, response_m, response_r = mix_values(wcsl, spectrum, star_interp)
        if rinput.degrade_resolution_method == 'none':
            sigma = 0
            self.logger.info('no broadening')
        elif rinput.degrade_resolution_method == 'fixed':
            sigma = rinput.sigma_resolution
            self.logger.info('fixed sigma=%3.0f', sigma)
        elif rinput.degrade_resolution_method == 'auto':
            self.logger.info('compute auto broadening')
            offset_broad, sigma_broad = compute_broadening(
                response_r.copy(), response_m.copy(), sigmalist=range(1, 101),
                remove_mean=False, frac_cosbell=0.10, zero_padding=50,
                fminmax=(0.003, 0.3), naround_zero=25, nfit_peak=21
            )
            sigma = sigma_broad
            self.logger.info('computed sigma=%3.0f', sigma)
        else:
            msg = f"'degrade_resolution_method' has value {rinput.degrade_resolution_method}"
            raise ValueError(msg)

        sens_raw = generate_sensitivity(final, spectrum, star_interp, extinc_interp, wl_cover1, wl_cover2, sigma)

        # Compute smoothed version
        self.logger.info('compute smoothed sensitivity')

        sens = sens_raw.copy()
        i_knots = rinput.smoothing_knots
        self.logger.debug(f'using adaptive spline with t={i_knots} interior knots')
        spl = AdaptiveLSQUnivariateSpline(x=wl_aa.value, y=sens_raw.data, t=i_knots)
        sens.data = spl(wl_aa.value)

        if self.intermediate_results:
            import matplotlib.pyplot as plt
            plt.plot(wl_aa, sens_raw.data, 'b')
            plt.plot(wl_aa, sens.data, 'r')
            plt.savefig('smoothed.png')
            plt.close()

        self.logger.info('end LCBStandardRecipe reduction')

        return self.create_result(
            reduced_image=reduced2d,
            final_rss=final,
            reduced_rss=origin,
            sky_rss=sky,
            star_spectrum=star_spectrum,
            master_sensitivity=sens,
            sensitivity_raw=sens_raw,
            fiber_ids=fiber_ids,
            sigma=sigma
        )
