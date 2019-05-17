#
# Copyright 2016-2019 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Combination routines"""



try:
    import contextlib2 as contextlib  # python 2.7 compatibility
except ImportError:
    import contextlib

from numina.array import combine
from numina.processing.combine import combine_imgs


def basic_processing_with_combination(
        rinput, reduction_flows,
        method=combine.mean, method_kwargs=None,
        errors=True, prolog=None):

    return basic_processing_with_combination_frames(
        rinput.obresult.frames, reduction_flows,
        method=method, method_kwargs=method_kwargs,
        errors=errors, prolog=prolog
    )


def basic_processing_with_combination_frames(
        frames, reduction_flows,
        method=combine.mean, method_kwargs=None,
        errors=True, prolog=None):
    """Perform basic reduction on set of DataFrames

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