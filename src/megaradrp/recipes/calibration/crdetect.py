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
from numina.array.crmasks.valid_parameters import VALID_COMBINATIONS


class CRMasksRecipe(MegaraBaseRecipe):
    """Recipe for generation of cosmic ray mask"""

    # Requirements
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()

    crmethod = Parameter('mm_pycosmic', description='Cosmic ray detection method')
    use_auxmedian = Parameter(False, description='Use auxiliary median in cosmic ray removal')
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
    save_postprocessed = Parameter(False, description='Save post-processed images after applying cosmic ray masks')
    debug = Parameter(False, description='Debug mode')

    # L.A. Cosmic parameters
    la_gain_apply = Parameter(False, description='Apply gain in L.A. Cosmic method')
    la_sigclip = Parameter([5.0,3.0], description='Sigma clipping for L.A. Cosmic method')
    la_sigfrac = Parameter([0.3,0.3], description='Fractional detection limit for L.A. Cosmic method')
    la_objlim = Parameter([5.0,5.0], description='Minimum contrast between Laplacian image and the fine structure image '
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
    la_padwidth = Parameter(10, description='Padding width for L.A. Cosmic method')

    # PyCosmic parameters
    pc_sigma_det = Parameter([5.0,3.0], description='Detection limit above the noise for PyCosmic method')
    pc_rlim = Parameter([1.2, 1.2], description="Detection threshold for PyCosmic method")
    pc_iterations = Parameter (4, description='Number of iterations for PyCosmic method')
    pc_fwhm_gauss_x = Parameter(2.5, description='FWHM of Gaussian smoothing kernel in X direction (pixels) for PyCosmic method')
    pc_fwhm_gauss_y = Parameter(2.5, description='FWHM of Gaussian smoothing kernel in Y direction (pixels) for PyCosmic method')
    pc_replace_box_x = Parameter(5, description='Median Box size in X direction for replacing masked pixels in PyCosmic method')
    pc_replace_box_y = Parameter(5, description='Median Box size in Y direction for replacing masked pixels in PyCosmic method')
    pc_replace_error = Parameter(1e6, description='Error value for bad pixels in PyCosmic method')
    pc_increase_radius = Parameter(0, description='Radius to increase the cosmic ray masks in PyCosmic method')
    pc_verbose = Parameter(False, description='Verbose mode for PyCosmic method')

    # deepCR parameters
    dc_mask = Parameter("ACS_WFC", description='Predefined mask for deepCR method')
    dc_threshold = Parameter(0.5, description='Threshold for deepCR method')
    dc_verbose = Parameter(False, description='Verbose mode for deepCR method')

    # Cosmic-CoNN parameters
    nn_model = Parameter("ground_imaging", description='Model for Cosmic-CoNN method')
    nn_threshold = Parameter(0.5, description='Threshold for Cosmic-CoNN method')
    nn_verbose = Parameter(False, description='Verbose mode for Cosmic-CoNN method')

    # Median-Minimum parameters
    mm_synthetic = Parameter('median', description='Type of synthetic image for median-minimum method')
    mm_hist2d_min_neighbors = Parameter(0, description='Minimum neighbors in 2D histogram for median-minimum method')
    mm_ydiag_max = Parameter(0, description='Maximum Y value for 2D histogram in median-minimum method (0=auto)')
    mm_dilation = Parameter(0, description='Dilation for median-minimum method')
    mm_xy_offsets = Parameter('none', description='List of (X,Y) offsets for image alignment')
    mm_crosscorr_region = Parameter('none', description='Region for cross-correlation alignment (FITS format)')
    mm_boundary_fit = Parameter('spline', description='Type of fit to CR detection boundary')
    mm_knots_splfit = Parameter(3, description='Number of knots for spline fit to CR detection boundary')
    mm_fixed_points_in_boundary = Parameter('none', description='List of fixed points in boundary for CR detection')
    mm_nsimulations = Parameter(10, description='Number of simulations for cosmic ray detection')
    mm_niter_boundary_extension = Parameter(3, description='Iterations for boundary extension')
    mm_weight_boundary_extension = Parameter(10, description='Weight for boundary extension')
    mm_threshold_rnoise = Parameter(0.0, description='Threshold for cosmic ray detection (in units of readout noise)')
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
                ldum = len(str(len(hdul_1im)))
                for i, hdul in enumerate(hdul_1im):
                    self.save_intermediate_img(hdul, f'preprocessed_{i+1:0{ldum}d}.fits')
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
                use_auxmedian=rinput.use_auxmedian,
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
                la_padwidth=rinput.la_padwidth,
                pc_sigma_det=rinput.pc_sigma_det,
                pc_rlim=rinput.pc_rlim,
                pc_iterations=rinput.pc_iterations,
                pc_fwhm_gauss_x=rinput.pc_fwhm_gauss_x,
                pc_fwhm_gauss_y=rinput.pc_fwhm_gauss_y,
                pc_replace_box_x=rinput.pc_replace_box_x,
                pc_replace_box_y=rinput.pc_replace_box_y,
                pc_replace_error=rinput.pc_replace_error,
                pc_increase_radius=rinput.pc_increase_radius,
                pc_verbose=rinput.pc_verbose,
                dc_mask=rinput.dc_mask,
                dc_threshold=rinput.dc_threshold,
                dc_verbose=rinput.dc_verbose,
                nn_model=rinput.nn_model,
                nn_threshold=rinput.nn_threshold,
                nn_verbose=rinput.nn_verbose,
                mm_synthetic=rinput.mm_synthetic,
                mm_hist2d_min_neighbors=rinput.mm_hist2d_min_neighbors,
                mm_ydiag_max=rinput.mm_ydiag_max,
                mm_dilation=rinput.mm_dilation,
                mm_xy_offsets=rinput.mm_xy_offsets,
                mm_crosscorr_region=rinput.mm_crosscorr_region,
                mm_boundary_fit=rinput.mm_boundary_fit,
                mm_knots_splfit=rinput.mm_knots_splfit,
                mm_fixed_points_in_boundary=rinput.mm_fixed_points_in_boundary,
                mm_nsimulations=rinput.mm_nsimulations,
                mm_niter_boundary_extension=rinput.mm_niter_boundary_extension,
                mm_weight_boundary_extension=rinput.mm_weight_boundary_extension,
                mm_threshold_rnoise=rinput.mm_threshold_rnoise,
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
            if rinput.save_postprocessed:
                for combination in VALID_COMBINATIONS:
                    self.logger.info("-" * 73)
                    self.logger.info(f'Apply {combination} to preprocessed images')
                    combined2d, _, _ = apply_crmasks(
                        list_arrays=arrays,
                        hdul_masks=hdul_masks,
                        combination=combination,
                        use_auxmedian=rinput.use_auxmedian,
                        dtype=rinput.dtype,
                        apply_flux_factor=True,
                        bias=None
                    )
                    self.save_intermediate_array(combined2d, f'preprocessed_{combination}.fits')

        result = self.create_result(crmasks=hdul_masks)
        self.logger.info('end MegaraCrDetection recipe')
        return result
