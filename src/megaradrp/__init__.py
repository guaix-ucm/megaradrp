#
# Copyright 2011-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""The MEGARA Data Reduction Pipeline."""

import logging

from megaradrp._version import __version__  # noqa: F401


# Top level NullHandler
logging.getLogger("megaradrp").addHandler(logging.NullHandler())
