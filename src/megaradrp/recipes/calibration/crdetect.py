#
# Copyright 2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#
import contextlib

from numina.core import Result, Parameter
import numpy as np

from megaradrp.core.recipe import MegaraBaseRecipe
from megaradrp.ntypes import CRMasks
import megaradrp.requirements as reqs
from numina.array.crmasks.apply_crmasks import apply_crmasks
from numina.array.crmasks.compute_crmasks import compute_crmasks


class CRMasksRecipe(MegaraBaseRecipe):
    """Recipe for generation of cosmic ray mask"""

    # Requirements
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()

    crmethod = Parameter('lacosmic', description='Cosmic ray detection method')
    use_lamedian = Parameter(False, description='Use L.A. Median in cosmic ray removal')
    flux_factor = Parameter('none', description='Flux factor for cosmic ray detection')
    flux_factor_regions = Parameter('none', description='Regions to compute flux factor (FITS format)')
    apply_flux_factor_to = Parameter('simulated',
                                     description='Apply flux factor to simulated images or to the original data')
    rnoise = Parameter(3.4, description='Readout noise in electrons')
    interactive = Parameter(True, description='Interactive mode for cosmic ray detection')
    dilation = Parameter(1, description='Dilation factor for cosmic ray masks')
    regions_to_be_skipped = Parameter('none',
                                      description='Regions to be skipped for cosmic ray detection (FITS format)')
    pixels_to_be_flagged_as_cr = Parameter(
        'none',
        description='List of (X,Y) coordinates of pixels to be masked (FITS criterium)'
    )
    pixels_to_be_ignored_as_cr = Parameter(
        'none',
        description='List of (X,Y) coordinates of pixels to be excluded from masking (FITS criterium)'
    )
    pixels_to_be_replaced_by_local_median = Parameter(
        'none',
        description='List of (X,Y,Xwidth,Ywidth) values with coordinates and median window shape'
    )
    dtype = Parameter('float32', description='Data type for cosmic ray masks')
    verify_cr = Parameter(False, description='Verify cosmic ray detection')
    semiwindow = Parameter(15, description='Semi-window size to display suspected cosmic rays')
    color_scale = Parameter('minmax', description='Color scale for plots')
    maxplots = Parameter(-1, description='Maximum number of plots with suspected cosmic rays to generate')
    save_preprocessed = Parameter(False, description='Save preprocessed images for cosmic ray detection')
    debug = Parameter(False, description='Debug mode')

    # L.A. Cosmic parameters
    la_gain_apply = Parameter(False, description='Apply gain in L.A. Cosmic method')
    la_sigclip = Parameter(5.0, description='Sigma clipping for L.A. Cosmic method')
    la_sigfrac = Parameter(0.3, description='Fractional detection limit for L.A. Cosmic method')
    la_objlim = Parameter(5.0, description='Minimum contrast between Laplacian image and the fine structure image '
                          'for L.A. Cosmic method')
    la_satlevel = Parameter(65535.0, description='Saturation level for L.A. Cosmic method')
    la_niter = Parameter(4, description='Number of iterations for L.A. Cosmic method')
    la_sepmed = Parameter(True, description='Use the separable median filter instead of the full median filter '
                          'for L.A. Cosmic method')
    la_cleantype = Parameter('meanmask', description='Cleaning type for L.A. Cosmic method')
    la_fsmode = Parameter('convolve', description='Filter mode for L.A. Cosmic method')
    la_psfmodel = Parameter('gaussxy', description='PSF model for L.A. Cosmic method')
    la_psffwhm_x = Parameter(2.5, description='PSF FWHM in X direction (pixels) for L.A. Cosmic method')
    la_psffwhm_y = Parameter(2.5, description='PSF FWHM in Y direction (pixels) for L.A. Cosmic method')
    la_psfsize = Parameter(11, description='PSF size (pixels) for L.A. Cosmic method')
    la_psfbeta = Parameter(4.765, description='Beta parameter for Moffat PSF model in L.A. Cosmic method')
    la_verbose = Parameter(False, description='Verbose mode for L.A. Cosmic method')

    # Median-Minimum parameters
    mm_xy_offsets = Parameter('none', description='List of (X,Y) offsets for image alignment')
    mm_crosscorr_region = Parameter('none', description='Region for cross-correlation alignment (FITS format)')
    mm_boundary_fit = Parameter('spline', description='Type of fit to CR detection boundary')
    mm_knots_splfit = Parameter(3, description='Number of knots for spline fit to CR detection boundary')
    mm_fixed_points_in_boundary = Parameter('none', description='List of fixed points in boundary for CR detection')
    mm_nsimulations = Parameter(10, description='Number of simulations for cosmic ray detection')
    mm_niter_boundary_extension = Parameter(3, description='Iterations for boundary extension')
    mm_weight_boundary_extension = Parameter(10, description='Weight for boundary extension')
    mm_threshold = Parameter(0.0, description='Threshold for cosmic ray detection')
    mm_minimum_max2d_rnoise = Parameter(5.0, description='Minimum max2d readout noise')
    mm_seed = Parameter(1234, description='Random seed for cosmic ray detection')

    # Results
    crmasks = Result(CRMasks)

    def run(self, rinput):
        """Execute the recipe."""
        self.logger.info('start MegaraCrDetection recipe')

        flow1 = self.init_filters(rinput, rinput.obresult.configuration)

        # The reduction flows are split in two parts:
        # 1. Overscan and trimming
        # 2. Bias, dark, and gain. Note that since we have not
        #    included the flatfield as a requirement of this recipe,
        #    the second part will not apply flat-fielding.
        if len(flow1) != 2:
            raise ValueError('Invalid reduction flow for cosmic ray detection')
        reduction_flow_ot, reduction_flow_1im = flow1
        frames = rinput.obresult.frames

        with contextlib.ExitStack() as stack:
            hduls = [stack.enter_context(dframe.open()) for dframe in frames]

            # Apply the first part of the reduction flow
            hdul_ot = [reduction_flow_ot(hdul) for hdul in hduls]

            # Apply the second part of the reduction flow
            hdul_1im = [reduction_flow_1im(hdul) for hdul in hdul_ot]

            # Save individual preprocessed images (and median) if requested
            if rinput.save_preprocessed:
                for i, hdul in enumerate(hdul_1im):
                    self.save_intermediate_img(hdul, f'preprocessed_{i+1}.fits')
                median_data = np.median([hdul[0].data for hdul in hdul_1im], axis=0)
                self.save_intermediate_array(median_data.astype(rinput.dtype), 'preprocessed_median.fits')

            arrays = [hdul[0].data for hdul in hdul_1im]
            self.logger.info(f'{len(arrays)} images to generate CR masks')

            # Generate the cosmic ray masks
            hdul_masks = compute_crmasks(
                list_arrays=arrays,
                gain=1.0,  # arrays are already in electrons
                rnoise=rinput.rnoise,  # readout noise in electrons
                bias=0.0,  # bias level was already removed
                crmethod=rinput.crmethod,
                use_lamedian=rinput.use_lamedian,
                flux_factor=rinput.flux_factor,
                flux_factor_regions=rinput.flux_factor_regions,
                apply_flux_factor_to=rinput.apply_flux_factor_to,
                interactive=rinput.interactive,
                dilation=rinput.dilation,
                regions_to_be_skipped=rinput.regions_to_be_skipped,
                pixels_to_be_flagged_as_cr=rinput.pixels_to_be_flagged_as_cr,
                pixels_to_be_ignored_as_cr=rinput.pixels_to_be_ignored_as_cr,
                pixels_to_be_replaced_by_local_median=rinput.pixels_to_be_replaced_by_local_median,
                dtype=rinput.dtype,
                verify_cr=rinput.verify_cr,
                semiwindow=rinput.semiwindow,
                color_scale=rinput.color_scale,
                maxplots=rinput.maxplots,
                debug=rinput.debug,
                la_gain_apply=rinput.la_gain_apply,
                la_sigclip=rinput.la_sigclip,
                la_sigfrac=rinput.la_sigfrac,
                la_objlim=rinput.la_objlim,
                la_satlevel=rinput.la_satlevel,
                la_niter=rinput.la_niter,
                la_sepmed=rinput.la_sepmed,
                la_cleantype=rinput.la_cleantype,
                la_fsmode=rinput.la_fsmode,
                la_psfmodel=rinput.la_psfmodel,
                la_psffwhm_x=rinput.la_psffwhm_x,
                la_psffwhm_y=rinput.la_psffwhm_y,
                la_psfsize=rinput.la_psfsize,
                la_psfbeta=rinput.la_psfbeta,
                la_verbose=rinput.la_verbose,
                mm_xy_offsets=rinput.mm_xy_offsets,
                mm_crosscorr_region=rinput.mm_crosscorr_region,
                mm_boundary_fit=rinput.mm_boundary_fit,
                mm_knots_splfit=rinput.mm_knots_splfit,
                mm_fixed_points_in_boundary=rinput.mm_fixed_points_in_boundary,
                mm_nsimulations=rinput.mm_nsimulations,
                mm_niter_boundary_extension=rinput.mm_niter_boundary_extension,
                mm_weight_boundary_extension=rinput.mm_weight_boundary_extension,
                mm_threshold=rinput.mm_threshold,
                mm_minimum_max2d_rnoise=rinput.mm_minimum_max2d_rnoise,
                mm_seed=rinput.mm_seed
            )

            # Update header (basic information)
            hdr = hdul_masks[0].header
            self.set_base_headers(hdr)

            # Update header (UUID of the individual images))
            for hdul in hdul_1im:
                hdr.add_history('---')
                hdr.add_history(f'Image {hdul[0].header["UUID"]}')
                history_entries = hdul[0].header.get('history', None)
                for entry in history_entries:
                    hdr.add_history(entry)
            hdr.add_history('---')

            # Update header with cosmic ray mask information
            for extname in ['MEDIANCR', 'MEANCRT'] + [f'CRMASK{i+1}' for i in range(len(arrays))]:
                hdr.add_history(f'Extension {extname}: {np.sum(hdul_masks[extname].data)} masked pixels')

            # Save the corrected preprocessed images if requested
            if rinput.save_preprocessed:
                for combination in ['mediancr', 'meancrt', 'meancr']:
                    self.logger.info(f'Apply {combination} to preprocessed images')
                    combined2d, _, _ = apply_crmasks(
                        list_arrays=arrays,
                        hdul_masks=hdul_masks,
                        combination=combination,
                        use_lamedian=rinput.use_lamedian,
                        dtype=rinput.dtype,
                        apply_flux_factor=True,
                        bias=None
                    )
                    self.save_intermediate_array(combined2d, f'preprocessed_{combination}.fits')

        result = self.create_result(crmasks=hdul_masks)
        self.logger.info('end MegaraCrDetection recipe')
        return result
