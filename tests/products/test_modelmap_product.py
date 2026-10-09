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
