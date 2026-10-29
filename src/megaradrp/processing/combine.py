#
# Copyright 2016-2026 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""Combination routines

The combination of frames with several reduction flows is in
numina.processing.combine since numina 0.41. The functions of this
module are deprecated and will be removed in megaradrp 0.24.
"""

import warnings

import numina.processing.combine as _combine

# Functions moved to numina.processing.combine
_MOVED = ("basic_processing_with_combination", "basic_processing_with_combination_frames")


def basic_processing_with_combination_frames_crmasks(
    frames, crmasks, reduction_flows, method, method_kwargs, errors=True, prolog=None
):
    """Perform basic reduction on set of DataFrames using a cosmic ray mask.

    .. deprecated:: 0.23
        Use numina.processing.combine.basic_processing_with_combination_frames
        with the argument crmasks.
    """
    warnings.warn(
        "basic_processing_with_combination_frames_crmasks is deprecated, it will be removed in 0.24, "
        "use numina.processing.combine.basic_processing_with_combination_frames with crmasks",
        DeprecationWarning,
        stacklevel=2,
    )
    return _combine.basic_processing_with_combination_frames(
        frames,
        reduction_flows,
        method=method,
        method_kwargs=method_kwargs,
        errors=errors,
        prolog=prolog,
        crmasks=crmasks,
    )


def __getattr__(name):
    if name in _MOVED:
        warnings.warn(
            f"megaradrp.processing.combine.{name} is deprecated, it will be removed in 0.24, "
            f"use numina.processing.combine.{name}",
            DeprecationWarning,
            stacklevel=2,
        )
        return getattr(_combine, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
