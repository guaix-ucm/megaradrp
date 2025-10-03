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

from numina.array import combine
from numina.processing.combine import combine_imgs


def basic_processing_with_combination(
        rinput, reduction_flows,
        method=combine.mean, method_kwargs=None,
        errors=True, prolog=None):

    if method.__name__ in ['mediancr', 'meancrt', 'meancr']:
        # Special case for combination using a cosmic ray mask
        return basic_processing_with_combination_frames_crmasks(
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


def basic_processing_with_combination_frames_crmasks(
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

    if crmasks is None:
        raise ValueError(f'Cosmic ray masks are required for {method.__name__}')
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

        result[0].header.add_history(f'Masks uuid:{crmasks[0].header["UUID"]}')

    return result
