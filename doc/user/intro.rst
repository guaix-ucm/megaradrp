
**********************
PyMegara Installation
**********************

PyMegara is distributed under GNU GPL, either version 3 of the License, 
or (at your option) any later version. See the file COPYING for details.

PyMegara requires the following packages installed in order to
be able to be installed and work properly:

 
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numina <http://guaix.fis.ucm.es/hg/numina/>`_

Additional packages are optionally required:
 - `nose <http://somethingaboutorange.com/mrl/projects/nose>`_ to run the tests
 - `sphinx`_ to build the documentation

Webpage: http://guaix.fis.ucm.es/hg/megara-drp/

Maintainer: sergiopr@fis.ucm.es

Stable version
--------------

The latest stable version of PyMegara can be downloaded from  
https://pypi.python.org/pypi/pymegara

To install PyMegara, use the standard installation procedure:::

    $ tar zxvf pymegara-X.Y.Z.tar.gz
    $ cd pymegara-X.Y.Z
    $ python setup.py install
    
The `install` command provides options to change the target directory. By 
default installation requires administrative privileges. The different 
installation options can be checked with::: 

   $ python setup.py install --help
   
Development version
-------------------

The development version can be checked out with:::

    $ hg clone http://guaix.fis.ucm.es/hg/megara-drp/

And then installed following the standard procedure:::

    $ cd megara-drp
    $ python setup.py install

Building the documentation
---------------------------
The PyMegara documentation is base on `sphinx`_. With the package installed, the 
html documentation can be built from the `doc` directory::

  $ cd doc
  $ make html
  
The documentation will be copied to a directory under `build/sphinx`.
  
The documentation can be built in different formats. The complete list will appear
if you type `make` 
  
.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org

