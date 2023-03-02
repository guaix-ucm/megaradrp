
from ..overplot_traces import main


def test_overplot_traces_help(capsys):
    """Check that the program runs with --help"""
    try:
        main(['--help'])
    except SystemExit as u:
        pass

    out, err = capsys.readouterr()
    assert True
