#
# Copyright 2016-2025 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Combination routines"""

import contextlib
import logging

from numina.array import combine
from numina.processing.combine import combine_imgs
import numpy as np


def basic_processing_with_combination(
        rinput, reduction_flows,
        method=combine.mean, method_kwargs=None,
        errors=True, prolog=None):

    if method.__name__ == 'mediancr':
        # Special case for combination using a cosmic ray mask
        return basic_processing_with_combination_frames_mediancr(
            rinput.obresult.frames, rinput.crmasks, reduction_flows,
            method=method, method_kwargs=method_kwargs
        )
    else:
        # General case for other combination methods
        return basic_processing_with_combination_frames(
            rinput.obresult.frames, reduction_flows,
            method=method, method_kwargs=method_kwargs,
            errors=errors, prolog=prolog
        )


def basic_processing_with_combination_frames(
        frames, reduction_flows,
        method=combine.mean, method_kwargs=None,
        errors=True, prolog=None):
    """Perform basic reduction on set of DataFrames.

    The reduction_flows are split in two parts.
    First part is performed in individual images (overscan and trimming)
    Then images are combined according to method and method_kwargs
    The resulting image is then processed with the
    second flow (bias, dark, gain and flat-fielding)
    """
    reduction_flow_ot, reduction_flow_1im = reduction_flows

    with contextlib.ExitStack() as stack:
        hduls = [stack.enter_context(dframe.open()) for dframe in frames]

        hdul_ot = [reduction_flow_ot(hdul) for hdul in hduls]

        hdu_combined = combine_imgs(hdul_ot, method=method, method_kwargs=method_kwargs,
                                    errors=errors, prolog=prolog)

        result = reduction_flow_1im(hdu_combined)

    return result


def basic_processing_with_combination_frames_mediancr(
        frames, crmasks, reduction_flows, method, method_kwargs,
        errors=True, prolog=None
):
    """Perform basic reduction on set of DataFrames using a cosmic ray mask.

    The reduction_flows are split in two parts:
    1: overscan and trimming;
    2: bias, dark, gain and flat-fielding

    In this case the two parts are applied to the individual images.
    After that, the images are combined using the requested cosmic ray mask.
    First part is performed in individual images (overscan and trimming)
    Then images are combined according to method and method_kwargs
    The resulting image is then processed with the
    second flow (bias, dark, gain and flat-fielding)
    """

    reduction_flow_ot, reduction_flow_1im = reduction_flows

    crmasks = crmasks.open()

    with contextlib.ExitStack() as stack:

        hduls = [stack.enter_context(dframe.open()) for dframe in frames]

        # apply overscan and trimming to individual images
        hdul_ot = [reduction_flow_ot(hdul) for hdul in hduls]

        # apply bias, dark, gain correction and flat-fielding to individual images
        hdul_otbg = [reduction_flow_1im(hdul) for hdul in hdul_ot]

        method_kwargs = method_kwargs or {}
        if 'dtype' not in method_kwargs:
            method_kwargs['dtype'] = 'float32'

        result = combine_imgs(hdul_otbg, method=method, method_kwargs=method_kwargs,
                              errors=errors, prolog=prolog, crmasks=crmasks)

    return result


def generate_crmask(
        rinput, reduction_flows,
        method=combine.maskscr, method_kwargs=None):
    """Generate a cosmic ray mask

    The reduction_flows are split in two parts. The two parts are
    applied to the individual images:
    - First part: overscan and trimming
    - Second part: bias, dark and gain correction. Important: flat-fielding
      is not applied here, as the calling recipe does not include it as a
      requirement.
    """

    _logger = logging.getLogger(__name__)

    reduction_flow_ot, reduction_flow_1im = reduction_flows
    frames = rinput.obresult.frames

    with contextlib.ExitStack() as stack:
        hduls = [stack.enter_context(dframe.open()) for dframe in frames]

        # apply overscan and trimming to invidual images
        hdul_ot = [reduction_flow_ot(hdul) for hdul in hduls]

        # apply bias, dark and gain correction to individual images
        hdul_otbg = [reduction_flow_1im(hdul) for hdul in hdul_ot]

        arrays = [hdul[0].data for hdul in hdul_otbg]
        _logger.info(f'{len(arrays)} images to generate CR masks using {method.__name__} method')

        method_kwargs = method_kwargs or {}
        if 'dtype' not in method_kwargs:
            method_kwargs['dtype'] = 'float32'

        # generate the cosmic ray masks
        hdul_masks = method(arrays, **method_kwargs)

        # update header
        for hdul in hdul_otbg:
            hdul_masks[0].header.add_history('---')
            hdul_masks[0].header.add_history(f'Image {hdul[0].header["UUID"]}')
            history_entries = hdul[0].header.get('history', None)
            for entry in history_entries:
                hdul_masks[0].header.add_history(entry)
        hdul_masks[0].header.add_history('---')
        for extname in ['MEDIANCR', 'MEANCR'] + [f'CRMASK{i+1}' for i in range(len(arrays))]:
            hdul_masks[0].header.add_history(f'Extension {extname}: {np.sum(hdul_masks[extname].data)} masked pixels')

    return hdul_masks
