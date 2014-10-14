
======================================
MEGARA DRP Deployment with Virtualenv
======================================

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

Install MEGARA DRP
-------------------

MEGARA DRP is registered in the Python Package Index. That means (among 
other things) that can be installed inside the environment with one command.


  (myenv) $ pip install megaradrp
  
The requirements of megaradrp will be downloaded and installed inside
the virtual environment. Once the installation is finished, you can check
by listing the installed recipes::

  (myenv) $ ./bin/numina show
  Bias Image: Recipe to process bias images
  Instrument: MEGARA
  Recipe: megara.drp.recipes.BiasRecipe
  Key: bias_image
  UUID: a7ea9c0c-76a6-4609-8413-0157102ec284
  --
  Dark current Image: Summary of Dark current Image
  Instrument: MEGARA
  Recipe: megara.drp.recipes.DarkRecipe
  Key: dark_image
  UUID: 5b15db30-7774-4823-9645-e8ca27a866af
  --
  ...


.. _virtualenv: http://pypi.python.org/pypi/virtualenv
