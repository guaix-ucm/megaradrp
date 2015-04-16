#####################
Installation
#####################
      
***********************
Requirements
***********************

The MEGARA Pipeline package requires the following packages installed in order to
be able to be installed and work properly:

 - `python 2.7 <https://www.python.org>`_
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numpy <http://www.numpy.org/>`_
 - `scipy <http://www.scipy.org/>`_
 - `astropy <http://www.astropy.org/>`_ >= 0.4
 - `numina <https://pypi.python.org/pypi/numina/>`_ >= 0.13

Additional packages are optionally required:

 - `py.test <http://pytest.org>`_ >= 2.5 to run the tests
 - `sphinx`_ to build the documentation


***********************
Installing MEGARA DRP
***********************

Using pip
---------
To install with pip, simply run:::

   pip install --no-deps megaradrp
   
.. note::

    The ``--no-deps`` flag is optional, but highly recommended if you already
    have Numpy installed, since otherwise pip will sometimes try to upgrade 
    your Numpy installation, which may not always be desired.

MEGARA DRP is registered in the Python Package Index. That means (among 
other things) that can be installed inside the environment with one command.


  (myenv) $ pip install megaradrp
  
The requirements of megaradrp will be downloaded and installed inside
the virtual environment. Once the installation is finished, you can check
by listing the installed recipes with the command line interface tool ``numina``::

  (myenv) $ ./bin/numina show-instruments
  INFO: Numina simple recipe runner version 0.13.0
  Instrument: MEGARA
   has configuration 'default'
   has pipeline 'default', version 1
   has pipeline 'experimental', version 1

***********************
Building from source
***********************


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
   
Development version
-------------------

The development version can be checked out with::

    $ hg clone http://guaix.fis.ucm.es/hg/megaradrp/

And then installed following the standard procedure::

    $ cd megaradrp
    $ python setup.py install

Building the documentation
---------------------------
The MEGARA DRP documentation is base on `sphinx`_. With the package 
installed, the html documentation can be built from the `doc` directory::

  $ cd doc
  $ make html
  
The documentation will be copied to a directory under `build/sphinx`.
  
The documentation can be built in different formats. The complete list will appear
if you type `make` 

***************************************
Deployment with Virtualenv
***************************************

`Virtualenv`_ is a tool to build isolated Python environments.

It's a great way to quickly test new libraries without cluttering your 
global site-packages or run multiple projects on the same machine which 
depend on a particular library but not the same version of the library.

Install virtualenv
------------------
I install it with the package system of my OS, so that it ends in my
global site-packages.

With Fedora/EL is just::

  $ sudo yum install python-virtualenv


Create virtual environment
--------------------------
Create the virtual environment enabling the packages already installed
in the global site-packages via the OS package system. Some requirements
(in particullar numpy and scipy) are difficult to build: they require
compiling and external C and FORTRAN libraries to be installed.

So the command is::

  $ virtualenv --system-site-packages myenv

If you need to create the virtualenv without global packages, drop the
system-site-packages flag.

Activate the environment
-------------------------
Once the environment is created, you need to activate it. Just change
directory into it and load with your command line interpreter the 
script bin/activate.

With bash::

  $ cd myenv
  $ . bin/activate
  (myenv) $

With csh/tcsh::

  $ cd myenv
  $ source bin/activate
  (myenv) $

Notice that the prompt changes once you are activate the environment. To 
deactivate it just type deactivate::

  (myenv) $ deactivate
  $ 

*********************
Installing MEGARA DFP
*********************

This section described how to install the MEGARA Pipeline inside
the GTC Control system.

Please refer to :ref:`Numina manual <numina:solaris10>` to install Numina
and its dependences under Solaris 10.




.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org

