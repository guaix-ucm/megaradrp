############
Installation
############

The easiest way to install megaradrp is using `pip`, the default Python
package manager. We support also conda.


megaradrp works with Python >= 3.6.


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

`megaradrp` can be installed with conda using a custom channel.

From the shell, run:::

 conda install -c conda-forge megaradrp


Building from source
====================

Obtaining the source
--------------------

You can obtain the source code from https://github.com/guaix-ucm/megaradrp

Or, if you happen to have git configured, you can clone the repository::

    git clone git://github.com/guaix-ucm/megaradrp.git


Installing from source
----------------------

You can install directly from the repository with::

    pip install git+https://github.com/guaix-ucm/megaradrp

If you have already downloaded the source code, you can install with::

    pip install .


.. note::
    There is a . in the command. Do not remove it.

This command will also download and install the dependencies.


You can also install with::

    python setup.py install
    
The `install` command provides options to change the target directory. By 
default installation requires administrative privileges. The different 
installation options can be checked with::

   python setup.py install --help


Dependencies
------------

The MEGARA Pipeline package requires the following packages installed in order to
be able to be installed and work properly:

 - `python <https://www.python.org>`_ either 2.7 or >= 3.5
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numpy <https://www.numpy.org/>`_ >= 1.7
 - `scipy <https://www.scipy.org/>`_
 - `astropy <https://www.astropy.org/>`_ >= 2.0
 - `numina <https://pypi.python.org/pypi/numina/>`_ >= 0.17
 - `scikit-image <https://scikit-image.org/>`_
 - `jsonschema <https://python-jsonschema.readthedocs.io/en/stable/>`_

If you install with pip, the depencies will be installed automatically.


Additional packages are optionally required:

 - `py.test <http://pytest.org>`_ >= 2.5 to run the tests
 - `sphinx`_ to build the documentation

Checking the installation
-------------------------
Once the installation is finished, you can check
by listing the installed recipes with the command line interface tool ``numina``::

  $ numina show-instruments
  INFO: Numina simple recipe runner version 0.22.0
  Instrument: MEGARA
   has configuration 'default'
   has pipeline 'default', version 1
   has pipeline 'experimental', version 1



Building the documentation
--------------------------
The MEGARA DRP documentation is based on `sphinx`_.

Additional packages required to create the documentation can be installed with::

    pip install .[docs]

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

After you have created the environmet, you can install megaradrp in it with
the pip command::

    (megara) $ pip install megaradrp



.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org
