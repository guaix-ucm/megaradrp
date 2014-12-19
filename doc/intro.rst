#####################
MEGARA DRP User Guide
#####################
      
This guide is intended as an introductory overview of MEGARA DRP and
explains how to install and make use of the most important features of
MEGARA DRP . For detailed reference documentation of the functions and
classes contained in the package, see the :ref:`reference`.
    
.. warning::

   This "User Guide" is still a work in progress; some of the material
   is not organized, and several aspects of MEGARA DRP are not yet covered
   sufficient detail.

***********************
MEGARA DRP Installation
***********************

MEGARA DRP is distributed under GNU GPL, either version 3 of the License, 
or (at your option) any later version. See the file LICENSE.txt for details.

MEGARA DRP requires the following packages installed in order to
be able to be installed and work properly:

 
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numina <http://guaix.fis.ucm.es/hg/numina/>`_

Additional packages are optionally required:
 - `py.test <http://pytest.org>`_ to run the tests
 - `sphinx`_ to build the documentation

Webpage: https://guaix.fis.ucm.es/projects/megara/

Maintainer: sergiopr@fis.ucm.es

Stable version
--------------

The latest stable version of PyMegara can be downloaded from  
https://pypi.python.org/pypi/megaradrp

To install PyMegara, use the standard installation procedure:::

    $ tar zxvf megaradrp-X.Y.Z.tar.gz
    $ cd megaradrp-X.Y.Z
    $ python setup.py install
    
The `install` command provides options to change the target directory. By 
default installation requires administrative privileges. The different 
installation options can be checked with::: 

   $ python setup.py install --help
   
Development version
-------------------

The development version can be checked out with:::

    $ hg clone http://guaix.fis.ucm.es/hg/megaradrp/

And then installed following the standard procedure:::

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
  
.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org

