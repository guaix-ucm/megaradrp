#
# Copyright 2011-2023 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import matplotlib

pytest_plugins = ["numina.testing.plugins"]

# The recipes draw plots, use a backend without display
matplotlib.use("Agg")
