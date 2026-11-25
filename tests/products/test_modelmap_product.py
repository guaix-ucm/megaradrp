"""Serialization of the model map"""

import copy

from megaradrp.products.modelmap import GeometricModel

STATE = {
    "fibid": 1,
    "boxid": 0,
    "start": 1,
    "stop": 4096,
    "model": {
        "model_name": "gaussbox",
        "params": {
            "mean": {
                "function": "spline1d",
                "params": [
                    [100.0, 100.0, 100.0, 100.0, 4000.0, 4000.0, 4000.0, 4000.0],
                    [86.48, 84.15, 90.11, 104.86, 0.0, 0.0, 0.0, 0.0],
                    3,
                ],
            },
        },
    },
}


def test_getstate_does_not_modify_the_model():
    model = GeometricModel.__new__(GeometricModel)
    model.__setstate__(copy.deepcopy(STATE))
    assert callable(model.model["params"]["mean"])

    state1 = model.__getstate__()
    state2 = model.__getstate__()

    # the parameters of the model are still functions
    assert callable(model.model["params"]["mean"])
    assert state1 == state2
    assert state1["model"]["params"]["mean"]["function"] == "spline1d"


def test_flags():
    """The flags of the quality control are saved; models saved before them have none"""
    model = GeometricModel.__new__(GeometricModel)
    model.__setstate__(copy.deepcopy(STATE))
    assert model.flags == []
    model.add_flag("not_parallel")
    model.add_flag("not_parallel")
    state = model.__getstate__()
    assert state["flags"] == ["not_parallel"]
    other = GeometricModel.__new__(GeometricModel)
    other.__setstate__(copy.deepcopy(state))
    assert other.flags == ["not_parallel"]
    # the flags do not make the model invalid
    assert other.valid
