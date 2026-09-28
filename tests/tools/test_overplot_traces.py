import sys

import pytest

from megaradrp.tools.overplot_traces import main


@pytest.mark.xfail(
    sys.version_info >= (3, 14),
    reason="Issue #334: argparse.FileType PendingDeprecationWarning on Python 3.14",
    strict=True,
)
def test_overplot_traces_help(capsys):
    """Check that the program runs with --help"""
    try:
        main(["--help"])
    except SystemExit:
        pass

    out, err = capsys.readouterr()
    assert True
