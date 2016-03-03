#####################
Installation
#####################
      
***********************
Requirements
***********************

The MEGARA Pipeline package requires the following packages installed 
in order to be able to be installed and work properly:

 - `python 2.7 <https://www.python.org>`_
 - `setuptools <http://peak.telecommunity.com/DevCenter/setuptools>`_
 - `numpy <http://www.numpy.org/>`_
 - `scipy <http://www.scipy.org/>`_
 - `astropy <http://www.astropy.org/>`_ >= 0.4
 - `numina <https://pypi.python.org/pypi/numina/>`_ >= 0.13

Additional packages are optionally required:

 - `py.test <http://pytest.org>`_ >= 2.5 to run the tests
 - `sphinx`_ to build the documentation


*********************
Installing MEGARA DFP
*********************

This section describes how to install the MEGARA Pipeline inside
the GTC Control system.

Please refer to :ref:`Numina manual <numina:solaris10>` to install Numina
and its dependences under Solaris 10.

Building from source
---------------------


The latest stable version of MEGARA DRP can be downloaded from
https://pypi.python.org/pypi/megaradrp

To install MEGARA DRP, use the standard installation procedure::

    $ tar zxvf megaradrp-X.Y.Z.tar.gz
    $ cd megaradrp-X.Y.Z
    $ python setup.py install --prefix GTC_PATH
    
The `install` command provides options to change the target directory. By
default installation requires administrative privileges. The different
installation options can be checked with::

   $ python setup.py install --help


.. _virtualenv: http://pypi.python.org/pypi/virtualenv
.. _sphinx: http://sphinx.pocoo.org

