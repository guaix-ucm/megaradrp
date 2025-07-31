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
from numina.array.crmasks import compute_crmasks


class CRMasksRecipe(MegaraBaseRecipe):
    """Recipe for generation of cosmic ray mask"""

    # Requirements
    master_bias = reqs.MasterBiasRequirement()
    master_dark = reqs.MasterDarkRequirement()
    master_bpm = reqs.MasterBPMRequirement()

    flux_factor = Parameter('auto', description='Flux factor for cosmic ray detection')
    knots_splfit = Parameter(3, description='Number of knots for spline fit')
    nsimulations = Parameter(10, description='Number of simulations for cosmic ray detection')
    niter_boundary_extension = Parameter(3, description='Iterations for boundary extension')
    weight_boundary_extension = Parameter(10, description='Weight for boundary extension')
    threshold = Parameter(0.0, description='Threshold for cosmic ray detection')
    minimum_max2d_rnoise = Parameter(5.0, description='Minimum max2d readout noise')
    interactive = Parameter(True, description='Interactive mode for cosmic ray detection')
    dilation = Parameter(1, description='Dilation factor for cosmic ray masks')
    dtype = Parameter('float32', description='Data type for cosmic ray masks')
    seed = Parameter(1234, description='Random seed for cosmic ray detection')
    plots = Parameter(False, description='Generate diagnostic plots')
    semiwindow = Parameter(15, description='Semi-window size to display suspected cosmic rays')
    color_scale = Parameter('minmax', description='Color scale for plots')
    maxplots = Parameter(10, description='Maximum number of plots with suspected cosmic rays to generate')

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

            for i, hdul in enumerate(hdul_1im):
                hdul.writeto(f'xxx_{i}.fits', overwrite=True)

            arrays = [hdul[0].data for hdul in hdul_1im]
            self.logger.info(f'{len(arrays)} images to generate CR masks')

            # Generate the cosmic ray masks
            hdul_masks = compute_crmasks(
                arrays,
                gain=1.0,    # arrays are already in electrons
                rnoise=3.4,  # readout noise in electrons
                bias=0.0,    # bias level was already removed
                flux_factor=rinput.flux_factor,
                knots_splfit=rinput.knots_splfit,
                nsimulations=rinput.nsimulations,
                niter_boundary_extension=rinput.niter_boundary_extension,
                weight_boundary_extension=rinput.weight_boundary_extension,
                threshold=rinput.threshold,
                minimum_max2d_rnoise=rinput.minimum_max2d_rnoise,
                interactive=rinput.interactive,
                dilation=rinput.dilation,
                dtype=rinput.dtype,
                seed=rinput.seed,
                plots=rinput.plots,
                semiwindow=rinput.semiwindow,
                color_scale=rinput.color_scale,
                maxplots=rinput.maxplots
            )

            # Update header
            hdr = hdul_masks[0].header
            self.set_base_headers(hdr)
            for hdul in hdul_1im:
                hdr.add_history('---')
                hdr.add_history(f'Image {hdul[0].header["UUID"]}')
                history_entries = hdul[0].header.get('history', None)
                for entry in history_entries:
                    hdr.add_history(entry)
            hdr.add_history('---')
            for extname in ['MEDIANCR', 'MEANCRT'] + [f'CRMASK{i+1}' for i in range(len(arrays))]:
                hdr.add_history(f'Extension {extname}: {np.sum(hdul_masks[extname].data)} masked pixels')

        result = self.create_result(crmasks=hdul_masks)
        self.logger.info('end MegaraCrDetection recipe')
        return result
