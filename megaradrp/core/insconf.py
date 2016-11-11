#
# Copyright 2016 Universidad Complutense de Madrid
#
# This file is part of Megara DRP
#
# Megara DRP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Megara DRP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Megara DRP.  If not, see <http://www.gnu.org/licenses/>.
#


from numina.core.pipeline import InstrumentConfiguration
import megaradrp.instrument.loader

class MegaraInstrumentConfiguration(InstrumentConfiguration):
    """The instrument configuration of MEGARA"""
    def __init__(self, name, values):
        super(MegaraInstrumentConfiguration, self).__init__(name, values)


def load_megara_instrument_configuration(key, value):
    return megaradrp.instrument.loader.build_instrument_config(value)
