from contextlib import nullcontext as does_not_raise

import pytest

from numina.core.dataholders import Requirement
import numina.core
import numina.dal.stored
from numina.core.validator import range_validator
from numina.types.datatype import PlainPythonType

from numina.exceptions import ValidationError
from numina.types.datatype import ListOfType
from numina.types.multitype import MultiType


class Dal(object):
    def search_parameter(self, name, type_, obsres, options=None):
        if name == 'req_null':
            return numina.dal.stored.StoredParameter(content=None)
        elif name == 'req_int':
            return numina.dal.stored.StoredParameter(content=2)
        elif name == 'smoothing_knots':
            return numina.dal.stored.StoredParameter(content=[1.0,2.0,3.0])
        else:
            raise numina.exceptions.NoResultFound('value not found')


@pytest.mark.parametrize(
    "example, expectation",
    [
        (1, does_not_raise()),
        ([1.0, 2.0, 3.0], does_not_raise()),
        pytest.param(
            [1.0], pytest.raises(ValidationError),
            marks=pytest.mark.xfail(reason="numina error")
        ),
    ]
)
def test_len_or_array_validate(example, expectation):
    """Test a type for positive integer or an array of floats"""

    len_or_array = MultiType(
        PlainPythonType(ref=3, validator=range_validator(minval=3)),
        ListOfType(PlainPythonType(ref=0.0), nmin=2)
    )

    with expectation:
        assert len_or_array.validate(example)


@pytest.mark.xfail(reason="numina error")
def test_len_or_array_validate2():
    """Test a type for positive integer or an array of floats"""

    len_or_array = MultiType(
        PlainPythonType(ref=3, validator=range_validator(minval=3)),
        ListOfType(PlainPythonType(ref=0.0), nmin=2)
    )

    with pytest.raises(ValidationError):
        len_or_array.validate(2)


def test_requirement():
    len_or_array = MultiType(
        PlainPythonType(ref=3, validator=range_validator(minval=3)),
        ListOfType(PlainPythonType(ref=0.0), nmin=3)
    )
    req = Requirement(len_or_array, 'List  of nodes or number of nodes',
                      default=3, optional=True, destination='smoothing_knots')


    dal = Dal()
    obsres = numina.core.ObservationResult()

    value = req.query(dal, obsres)
    assert value == [1.0,2.0,3.0]
