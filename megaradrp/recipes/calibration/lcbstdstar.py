#
# Copyright 2011-2019 Universidad Complutense de Madrid
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
import matplotlib.pyplot as plt
import numpy
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import interp1d

from numina.array.numsplines import AdaptiveLSQUnivariateSpline
from numina.array.wavecalib.crosscorrelation import periodic_corr1d
from numina.core import Result, Parameter
from numina.core.requirements import Requirement
from numina.core.validator import range_validator
from numina.types.array import ArrayType

from megaradrp.processing.extractobj import extract_star, generate_sensitivity
from megaradrp.recipes.scientific.base import ImageRecipe
from megaradrp.types import ProcessedRSS, ProcessedFrame, ProcessedSpectrum
from megaradrp.types import ReferenceSpectrumTable, ReferenceExtinctionTable
from megaradrp.types import MasterSensitivity


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
    `master_traces` and resampled according to the wavelength calibration in
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
    position = Requirement(list, "Position of the reference object", default=(0, 0))
    nrings = Parameter(3, "Number of rings to extract the star",
                       validator=range_validator(minval=1))
    reference_spectrum = Requirement(ReferenceSpectrumTable, "Spectrum of reference star")
    reference_spectrum_velocity = Parameter(0.0, 'Radial velocity (km/s) of reference spectrum')
    reference_extinction = Requirement(ReferenceExtinctionTable, "Reference extinction")
    # sigma_auto = Parameter(False, 'sigma Gaussian is automatically computed')
    sigma_resolution = Parameter(20.0, 'sigma Gaussian filter to degrade resolution ')

    reduced_image = Result(ProcessedFrame)
    final_rss = Result(ProcessedRSS)
    reduced_rss = Result(ProcessedRSS)
    sky_rss = Result(ProcessedRSS)
    star_spectrum = Result(ProcessedSpectrum)
    master_sensitivity = Result(MasterSensitivity)
    sensitivity_raw = Result(MasterSensitivity)
    fiber_ids = Result(ArrayType(fmt='%d'))
    sigma = Result(float)

    def run(self, rinput):

        self.logger.info('starting LCBStandardRecipe reduction')

        # Create InstrumentModel
        # ins1 = rinput.obresult.configuration
        #
        reduced2d, rss_data = super(LCBStandardRecipe, self).base_run(rinput)
        # tags = rinput.obresult.tags
        ins2 = rinput.obresult.profile
        ins2.configure_with_image(rss_data)
        self.logger.info('start sky subtraction')
        final, origin, sky = self.run_sky_subtraction(rss_data, rinput.ignored_sky_bundles)
        self.logger.info('end sky subtraction')

        # 1 + 6  for first ring
        # 1 + 6  + 12  for second ring
        # 1 + 6  + 12  + 18 for third ring
        # 1 + 6 * Sum_i=0^n =  1 + 3 * n * (n +1)
        # Using three rings around central point
        self.logger.debug('adding %d nrings', rinput.nrings)
        npoints = 1 + 3 * rinput.nrings * (rinput.nrings +1)
        self.logger.debug('adding %d fibers', npoints)

        fiberconf = self.datamodel.get_fiberconf(final)
        spectra_pack = extract_star(final, rinput.position, npoints,
                                    fiberconf, logger=self.logger)

        spectrum, colids, wl_cover1, wl_cover2 = spectra_pack
        star_spectrum = fits.PrimaryHDU(spectrum, header=final[0].header)

        rad_vel = rinput.reference_spectrum_velocity * u.km / u.s
        factor = 1 + rad_vel / const.c
        star_interp = interp1d(rinput.reference_spectrum[:,0] / factor,
                               rinput.reference_spectrum[:,1])

        extinc_interp = interp1d(rinput.reference_extinction[:, 0],
                               rinput.reference_extinction[:, 1])

        fiber_ids = [colid + 1 for colid in colids]

        wcsl = astropy.wcs.WCS(final[0].header)
        wl_aa, response_m, response_r = mix_values(wcsl, spectrum, star_interp)

        if rinput.sigma_resolution < 0:
            self.logger.info('compute auto broadening')
            offset_broad, sigma_broad = compute_broadening(
                response_r, response_m, sigmalist=range(1, 101),
                remove_mean=False, frac_cosbell=0.10, zero_padding=50,
                fminmax=(0.003, 0.3), naround_zero=25, nfit_peak=21
            )
            sigma = sigma_broad
            self.logger.info('computed sigma=%3.0f', sigma)
        else:
            sigma = rinput.sigma_resolution
            self.logger.info('sigma=%3.0f', sigma)

        sens_raw = generate_sensitivity(final, spectrum, star_interp, extinc_interp, wl_cover1, wl_cover2, sigma)

        # Compute smoothed version
        self.logger.info('compute smoothed sensitivity')

        sens = sens_raw.copy()
        i_knots = 3
        self.logger.debug('using sdaptive spline with t=%d interior knots', i_knots)
        spl = AdaptiveLSQUnivariateSpline(x=wl_aa.value, y=sens_raw.data, t=i_knots)
        sens.data = spl(wl_aa.value)

        if self.intermediate_results:
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
            sensitivity_raw=sens,
            fiber_ids=fiber_ids,
            sigma=sigma
        )


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
