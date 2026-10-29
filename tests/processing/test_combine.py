"""The combination functions moved to numina.processing.combine"""

import pytest

import numina.processing.combine as numina_combine

import megaradrp.processing.combine as combine


@pytest.mark.parametrize("name", ["basic_processing_with_combination", "basic_processing_with_combination_frames"])
def test_moved_functions_deprecated(name):
    with pytest.warns(DeprecationWarning, match=f"use numina.processing.combine.{name}"):
        function = getattr(combine, name)
    assert function is getattr(numina_combine, name)


def test_crmasks_deprecated(monkeypatch):
    calls = []

    def record(*args, **kwargs):
        calls.append((args, kwargs))
        return "result"

    monkeypatch.setattr(numina_combine, "basic_processing_with_combination_frames", record)
    with pytest.warns(DeprecationWarning, match="with crmasks"):
        result = combine.basic_processing_with_combination_frames_crmasks(
            ["frame"], "masks", ["flow1", "flow2"], method=min, method_kwargs={}
        )
    assert result == "result"
    assert calls == [
        (
            (["frame"], ["flow1", "flow2"]),
            {"method": min, "method_kwargs": {}, "errors": True, "prolog": None, "crmasks": "masks"},
        )
    ]


def test_other_names():
    with pytest.raises(AttributeError):
        combine.other
