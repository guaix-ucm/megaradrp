import numina.types
from numina.array.wavecalib.solutionarc import (
    WavecalFeature,
    CrLinear,
    SolutionArcCalibration,
)

from megaradrp.products import wavecalibration as wcal


def create_solution(orig):

    # Create Features
    features = []
    for feature in orig["features"]:
        features.append(WavecalFeature(**feature))

    cr_linear = CrLinear(**orig["cr_linear"])
    mm = SolutionArcCalibration(
        features, orig["coeff"], orig["residual_std"], cr_linear
    )
    return mm


orig = {
    "cr_linear": {
        "cdelt": 3.0,
        "crmax": 4600,
        "crmin": 2300,
        "crpix": 1200,
        "crval": 12,
    },
    "features": [
        {
            "category": "A",
            "xpos": 100,
            "reference": 3210,
            "wavelength": 1,
            "line_ok": True,
            "ypos": 0,
            "peak": 0,
            "funcost": 12.0,
            "fwhm": 0,
            "lineid": 1,
        },
        {
            "category": "A",
            "xpos": 150,
            "reference": 3310,
            "wavelength": 2,
            "line_ok": True,
            "ypos": 0,
            "peak": 0,
            "funcost": 12.0,
            "fwhm": 0,
            "lineid": 2,
        },
        {
            "category": "C",
            "xpos": 250,
            "reference": 3410,
            "wavelength": 3,
            "line_ok": True,
            "ypos": 0,
            "peak": 0,
            "funcost": 13.0,
            "fwhm": 0,
            "lineid": 3,
        },
    ],
    "coeff": [1.0, 0.1],
    "residual_std": 1.0,
}


def create_test_wavecalib():
    instrument = "TEST1"
    tags = {"insmode": "LCB", "vph": "LR-I"}
    uuid = "123456789"
    data = wcal.WavelengthCalibration(instrument=instrument)
    data.tags = tags
    data.uuid = uuid
    data.error_fitting = []
    data.missing_fibers = []
    data.total_fibers = 623
    meta_info = wcal.WavelengthCalibration.create_meta_info()
    meta_info["instrument_name"] = instrument
    meta_info["creation_date"] = data.meta_info["creation_date"]
    contents = []

    for key in range(10):
        fibid = key + 1
        solution = create_solution(orig)
        fibersolution = wcal.FiberSolutionArcCalibration(fibid, solution)
        data.contents.append(fibersolution)
        contents.append(data.contents[-1].__getstate__())

    for key in range(100, 101):
        fibid = key + 1
        solution = create_solution(orig)
        fibersolution = wcal.FiberSolutionArcCalibration(fibid, solution)
        data.contents.append(fibersolution)
        contents.append(data.contents[-1].__getstate__())

    state = dict(
        instrument=instrument,
        tags=tags,
        uuid=uuid,
        total_fibers=623,
        missing_fibers=[],
        error_fitting=[],
        meta_info=meta_info,
        type=data.name(),
        contents=contents,
        global_offset=[0.0],
        type_fqn="megaradrp.products.wavecalibration.WavelengthCalibration",
        quality_control=numina.types.qc.QC.UNKNOWN,
    )
    return data, state
