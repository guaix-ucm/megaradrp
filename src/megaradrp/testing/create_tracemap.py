import numina.types

from megaradrp.products import tracemap as tm


def create_test_tracemap():
    instrument = "TEST1"
    tags = {}
    uuid = "123456789"
    data = tm.TraceMap(instrument=instrument)
    data.tags = tags
    data.uuid = uuid
    data.total_fibers = 623
    data.expected_range = [2, 4092]
    data.ref_column = 2001
    meta_info = tm.TraceMap.create_meta_info()
    meta_info["instrument_name"] = instrument
    meta_info["creation_date"] = data.meta_info["creation_date"]
    state = dict(
        instrument=instrument,
        tags=tags,
        uuid=uuid,
        error_fitting=[],
        missing_fibers=[],
        total_fibers=623,
        meta_info=meta_info,
        contents=[],
        type_fqn="megaradrp.products.tracemap.TraceMap",
        boxes_positions=[],
        type=data.name(),
        ref_column=2001,
        expected_range=[2, 4092],
        global_offset=[0.0],
        quality_control=numina.types.qc.QC.UNKNOWN,
    )

    return data, state


def create_test_tracemap2():
    instrument = "MEGARA"
    tags = {"insmode": "LCB", "vph": "LR-R"}
    uuid = "123456789"
    data = tm.TraceMap(instrument=instrument)
    data.tags = tags
    data.uuid = uuid
    data.total_fibers = 623
    data.expected_range = [2, 4092]
    data.ref_column = 2001
    data.contents = [
        tm.GeometricTrace(fibid, 1, 4, 4090, fitparms=[200 + fibid * 3.5, 0.0])
        for fibid in range(1, 624)
    ]

    return data
