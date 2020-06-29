############
Installation
############
      
************
Requirements
************

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

Additional packages are optionally required:

 - `py.test <http://pytest.org>`_ >= 2.5 to run the tests
 - `sphinx`_ to build the documentation


*********************
Installing MEGARA DRP
*********************

Using Conda
===========

Using pip
=========
To install with pip, simply run:::

   pip install --no-deps megaradrp
   
.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to upgrade 
    your Numpy installation, which may not always be desired.


`megaradrp` can be installed with conda using a custom channel.

From the shell, execute:::

 conda install -c conda-forge megaradrp


Building from source
====================

The latest stable version of MEGARA DRP can be downloaded from
https://pypi.python.org/pypi/megaradrp

To install MEGARA DRP, use the standard installation procedure::

    $ tar zxvf megaradrp-X.Y.Z.tar.gz
    $ cd megaradrp-X.Y.Z
    $ python setup.py install
    
The `install` command provides options to change the target directory. By 
default installation requires administrative privileges. The different 
installation options can be checked with::

   $ python setup.py install --help


Checking the installation
=========================
Once the installation is finished, you can check
by listing the installed recipes with the command line interface tool ``numina``::

  (myenv) $ ./bin/numina show-instruments
  INFO: Numina simple recipe runner version 0.13.0
  Instrument: MEGARA
   has configuration 'default'
   has pipeline 'default', version 1
   has pipeline 'experimental', version 1


Development version
-------------------

The development version can be checked out with::

    $ git clone https://github.com/guaix-ucm/megaradrp.git

And then installed following the standard procedure::

    $ cd megaradrp
    $ python setup.py install

Building the documentation
--------------------------
The MEGARA DRP documentation is base on `sphinx`_. With the package 
installed, the html documentation can be built from the `doc` directory::

  $ cd doc
  $ make html
  
The documentation will be copied to a directory under `build/sphinx`.
  
The documentation can be built in different formats. The complete list will appear
if you type `make` 


Deployment with Virtualenv
==========================

`Virtualenv`_ is a tool to build isolated Python environments.

It's a great way to quickly test new libraries without cluttering your 
global site-packages or run multiple projects on the same machine which 
depend on a particular library but not the same version of the library.

Since Python version 3.3, there is also a module in the standard library 
called `venv` with roughly the same functionality.

Create virtual environment
--------------------------
In order to create a virtual environment called e.g. megara using `venv`:

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

.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org

