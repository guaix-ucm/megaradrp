############
Installation
############

The easiest way to install megaradrp is using `pip`, the default Python
package manager, in a virtual environment or in a conda environment.


megaradrp works with Python >= 3.11, and requires numina >= 0.41.


******************
Using PyPI and pip
******************

To install with pip, simply run:::

   pip install megaradrp

The latest stable version of MEGARA DRP can be downloaded from
https://pypi.python.org/pypi/megaradrp

With pip, we recommend to work in a virtual environment, see :ref:`deploy_venv`.


************
Using Conda
************

Create and activate a conda environment, and install megaradrp in it with
pip:::

   conda create --name megara python=3.12
   conda activate megara
   pip install megaradrp

The packages of megaradrp and numina in the `conda-forge
<https://conda-forge.org>`_ channel are no longer updated: they are older
versions, and some dependencies of numina, such as teareduce, are not in
conda-forge. Install megaradrp and numina with pip, also in conda
environments.


Building from source
====================

Obtaining the source
--------------------

You can obtain the source code from https://github.com/guaix-ucm/megaradrp

Or, if you happen to have git configured, you can clone the repository::

    git clone https://github.com/guaix-ucm/megaradrp.git


Installing from source
----------------------

You can install directly from the repository with::

    pip install git+https://github.com/guaix-ucm/megaradrp

If you have already downloaded the source code, you can install with::

    pip install .


.. note::
    There is a . in the command. Do not remove it.

This command will also download and install the dependencies.


Dependencies
------------

The MEGARA Pipeline package requires the following packages installed in order to
be able to be installed and work properly:

 - `python <https://www.python.org>`_ >= 3.10
 - `numina <https://pypi.org/project/numina/>`_ >= 0.40
 - `numpy <https://numpy.org/>`_
 - `scipy <https://scipy.org/>`_ >= 1.12
 - `astropy <https://www.astropy.org/>`_ >= 6.1.7
 - `scikit-image <https://scikit-image.org/>`_
 - `jsonschema <https://python-jsonschema.readthedocs.io/en/stable/>`_
 - `matplotlib <https://matplotlib.org/>`_
 - `pooch <https://www.fatiando.org/pooch/>`_

The complete list is in the file ``pyproject.toml`` of the source code.
If you install with pip, the dependencies will be installed automatically.


Additional packages are optionally required:

 - `pytest <https://pytest.org>`_ to run the tests (``pip install ".[test]"``)
 - `sphinx`_ to build the documentation (``pip install ".[docs]"``)

Checking the installation
-------------------------
Once the installation is finished, you can check
by listing the installed recipes with the command line interface tool ``numina``::

  $ numina show-instruments
  INFO: Numina simple recipe runner version 0.40.0
  Instrument: MEGARA
   version is '0.22'
   has profile 'Configuration at LICA' uuid=9a86b2b2-3f7d-48ec-8f4f-3780ec967c90
   has profile 'MEGARA@GTC after November 2025' uuid=9f3fecbf-376c-47ae-a188-313b7d829104
   has profile 'Configuration for the Simulator' uuid=66f2283e-3049-4d4b-8ef1-14d62fcb611d
   has profile 'MEGARA@GTC between 2017 and 2024' uuid=ca3558e3-e50d-4bbc-86bd-da50a0998a48
   has profile 'Test component' uuid=4fd05b24-2ed9-457b-b563-a3c618bb1d4c
   has datamodel 'megaradrp.datamodel.MegaraDataModel'
   has pipeline 'default', version 1



Building the documentation
--------------------------
The MEGARA DRP documentation is based on `sphinx`_.

Additional packages required to create the documentation can be installed with::

    pip install ".[docs]"

With these packages
installed, the html documentation can be built from the `doc` directory::

  $ cd doc
  $ make html

The documentation will be copied to a directory under `_build/html`.

The documentation can be built in different formats. The complete list will appear
if you type `make`


.. _deploy_venv:

Deployment in a virtual environment
===================================

`Virtualenv`_ is a tool to build isolated Python environments.

It's a great way to quickly test new libraries without cluttering your
global site-packages or run multiple projects on the same machine which
depend on a particular library but not the same version of the library.

Since Python version 3.3, there is also a module in the standard library
called `venv` with roughly the same functionality.

Create virtual environment
--------------------------
In order to create a virtual environment called e.g. megara using `venv`, run::

  $ python3 -m venv megara

Activate the environment
------------------------
Once the environment is created, you need to activate it. Just change
directory into it and source the script `bin/activate`.

With bash::

  $ cd megara
  $ . bin/activate
  (megara) $

With csh/tcsh::

  $ cd megara
  $ source bin/activate
  (megara) $

Notice that the prompt changes once you are activate the environment. To
deactivate it just type deactivate::

  (megara) $ deactivate
  $

After you have created the environment, you can install megaradrp in it with
the pip command::

    (megara) $ pip install megaradrp



.. _virtualenv: https://virtualenv.pypa.io/
.. _sphinx: https://www.sphinx-doc.org
