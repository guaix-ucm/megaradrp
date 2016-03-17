################
Testing the DFP
################

This section describes the testing framework and options for testing MEGARA DFP

*******************
Running unit tests
*******************

MEGARA DFP uses `py.test <http://pytest.org>`_ as its testing framework.

As MEGARA DFP does not contain C/Cython extensions, the tests can be run
directly in the source code, as::

    cd megaradrp-0.4.0
    cd src
    py.test megaradrp
    
Some of the tests rely on data downloaded from a server. These tests are
skipped by default. To enable them run instead::

    py.test --run-remote megaradrp

The reduction recipes are tested with remote data. Each recipe is run in
a directory created under the default ``$TMPDIR``, which is based on
the user temporal directory. The base of the created directories can be changed
with the option ``--basetemp=dir``::

    py.test --basetemp=/home/spr/test100 --run-remote megaradrp
