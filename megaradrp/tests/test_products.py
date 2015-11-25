#
# Copyright 2015 Universidad Complutense de Madrid
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

from megaradrp.products import TraceMap
from tempfile import NamedTemporaryFile
import pytest

@pytest.mark.xfail
def test_traceMap(benchmark=None):
    my_obj = TraceMap()
    my_obj._datatype_load('')

def test_traceMap(benchmark=None):
    my_obj = TraceMap()
    my_file =NamedTemporaryFile()
    print (my_file.name)
    my_open_file = my_obj._datatype_load(my_file.name)
    print (my_open_file)

if __name__ == "__main__":
    test_traceMap()