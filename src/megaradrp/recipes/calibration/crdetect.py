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
from numina.array.crmasks import compute_crmasks, apply_crmasks


class CRMasksRecipe(MegaraBaseRecipe):
    """Recipe for generation of cosmic ray mask"""

    # Requirements
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()

    crmethod = Parameter('lacosmic', description='Cosmic ray detection method')
    flux_factor = Parameter('none', description='Flux factor for cosmic ray detection')
    rnoise = Parameter(3.4, description='Readout noise in electrons')
    interactive = Parameter(True, description='Interactive mode for cosmic ray detection')
    dilation = Parameter(1, description='Dilation factor for cosmic ray masks')
    pixels_to_be_masked = Parameter(
        'none',
        description='List of (X,Y) coordinates of pixels to be masked (FITS criterium)'
    )
    pixels_to_be_excluded = Parameter(
        'none',
        description='List of (X,Y) coordinates of pixels to be excluded from masking (FITS criterium)'
    )
    dtype = Parameter('float32', description='Data type for cosmic ray masks')
    verify_cr = Parameter(False, description='Verify cosmic ray detection')
    semiwindow = Parameter(15, description='Semi-window size to display suspected cosmic rays')
    color_scale = Parameter('minmax', description='Color scale for plots')
    maxplots = Parameter(-1, description='Maximum number of plots with suspected cosmic rays to generate')
    save_preprocessed = Parameter(False, description='Save preprocessed images for cosmic ray detection')
    debug = Parameter(False, description='Debug mode')

    # L.A. Cosmic parameters
    la_sigclip = Parameter(5.0, description='Sigma clipping for L.A. Cosmic method')
    la_psffwhm_x = Parameter(2.5, description='PSF FWHM in X direction (pixels) for L.A. Cosmic method')
    la_psffwhm_y = Parameter(2.5, description='PSF FWHM in Y direction (pixels) for L.A. Cosmic method')
    la_fsmode = Parameter('convolve', description='Filter mode for L.A. Cosmic method')
    la_psfmodel = Parameter('gaussxy', description='PSF model for L.A. Cosmic method')
    la_psfsize = Parameter(11, description='PSF size (pixels) for L.A. Cosmic method')

    # Simulated boundary parameters
    sb_crosscorr_region = Parameter('none', description='Region for cross-correlation alignment (FITS format)')
    sb_boundary_fit = Parameter('spline', description='Type of fit to CR detection boundary')
    sb_knots_splfit = Parameter(3, description='Number of knots for spline fit to CR detection boundary')
    sb_fixed_points_in_boundary = Parameter('none', description='List of fixed points in boundary for CR detection')
    sb_nsimulations = Parameter(10, description='Number of simulations for cosmic ray detection')
    sb_niter_boundary_extension = Parameter(3, description='Iterations for boundary extension')
    sb_weight_boundary_extension = Parameter(10, description='Weight for boundary extension')
    sb_threshold = Parameter(0.0, description='Threshold for cosmic ray detection')
    sb_minimum_max2d_rnoise = Parameter(5.0, description='Minimum max2d readout noise')
    sb_seed = Parameter(1234, description='Random seed for cosmic ray detection')

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
                flux_factor=rinput.flux_factor,
                interactive=rinput.interactive,
                dilation=rinput.dilation,
                pixels_to_be_masked=rinput.pixels_to_be_masked,
                pixels_to_be_excluded=rinput.pixels_to_be_excluded,
                dtype=rinput.dtype,
                verify_cr=rinput.verify_cr,
                semiwindow=rinput.semiwindow,
                color_scale=rinput.color_scale,
                maxplots=rinput.maxplots,
                debug=rinput.debug,
                la_sigclip=rinput.la_sigclip,
                la_psffwhm_x=rinput.la_psffwhm_x,
                la_psffwhm_y=rinput.la_psffwhm_y,
                la_fsmode=rinput.la_fsmode,
                la_psfmodel=rinput.la_psfmodel,
                la_psfsize=rinput.la_psfsize,
                sb_crosscorr_region=rinput.sb_crosscorr_region,
                sb_boundary_fit=rinput.sb_boundary_fit,
                sb_knots_splfit=rinput.sb_knots_splfit,
                sb_fixed_points_in_boundary=rinput.sb_fixed_points_in_boundary,
                sb_nsimulations=rinput.sb_nsimulations,
                sb_niter_boundary_extension=rinput.sb_niter_boundary_extension,
                sb_weight_boundary_extension=rinput.sb_weight_boundary_extension,
                sb_threshold=rinput.sb_threshold,
                sb_minimum_max2d_rnoise=rinput.sb_minimum_max2d_rnoise,
                sb_seed=rinput.sb_seed
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
                        dtype=rinput.dtype,
                        apply_flux_factor=True,
                        bias=None
                    )
                    self.save_intermediate_array(combined2d, f'preprocessed_{combination}.fits')

        result = self.create_result(crmasks=hdul_masks)
        self.logger.info('end MegaraCrDetection recipe')
        return result
