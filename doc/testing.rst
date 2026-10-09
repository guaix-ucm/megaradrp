########
Testing
########

This section describes the testing framework and options for testing MEGARA DRP

**************
Running tests
**************

MEGARA DRP uses `pytest <http://pytest.org>`_ as its testing framework.
We require also the package `pytest-remotedata` to control access to online
resources during testing. The optional dependencies for testing are installed
with::

    pip install -e ".[test]"

The tests are in the directory ``tests``, and are run from the root of the
source tree::

    pytest

Some of the tests rely on data downloaded from a server. These tests are
skipped by default. To enable them run instead::

    pytest --remote-data=any

The reduction recipes are tested with remote data. Each recipe is run in
a directory created under the default ``$TMPDIR``, which is based on
the user temporal directory. The base of the created directories can be changed
with the option ``--basetemp=dir``::

    pytest --basetemp=/home/spr/test100 --remote-data=any

The tests can be run in parallel with
`pytest-xdist <https://pytest-xdist.readthedocs.io>`_, that is included in the
optional dependencies for testing. The tests of a file must be run in the same
process, as they share the downloaded data and can use the results of previous
tests, so use the option ``--dist loadfile``::

    pytest -n auto --dist loadfile --remote-data=any
