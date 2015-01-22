
#####################
Running the pipeline
#####################

The MEGARA DRP is run through a command line interface
provided by numina.

:ref:`numina:creation`

The run mode of numina requires:
 
  * A observation result file in YAML_ format
  * A requirements file in YAML_ format 
  * The raw images obtained in the observing block
  * The calibrations required by the recipe
 
 
---------------------------------
 Format of the observation result
---------------------------------

The contents of the file is a serialized dictionary with the
following keys:

*id*: not required, integer, defaults to 1
    Unique identifier of the observing block

*instrument*: required, string
    Name of the instrument, as it is returned by ``numina show-instruments``

*mode*: required, string
    Name of the observing mode, as returned by ``numina show-modes``

*children*: not required, list of integers, defaults to empty list
    Identifications of nested observing blocks

*frames*: required, list of strings
    List of images names

This is an example of the observation result file

.. code-block:: yaml

   id: 21
   instrument: MEGARA
   mode: dark_image
   frames:
   - r0121.fits
   - r0122.fits
   - r0123.fits
   - r0124.fits
   - r0125.fits
   - r0126.fits
   - r0127.fits
   - r0128.fits
   - r0129.fits
   - r0130.fits
   - r0131.fits
   - r0132.fits
   
---------------------------------
 Format of the requirements file
---------------------------------

This file contains configuration parameters for the recipes that
are not related to the particular instrument used.

The contents of the file are serialized as a dictionary with the
following keys:

requirements: required, dictionary
    A dictionary of parameter names and values.

products: optional, dictionary
    A dictionary with names for the products

logger: optional, dictionary
    A dictionary used to configure the custom file logger

Example requirements file

.. code-block:: yaml

   requirements:
     master_bias: master_bias-1.fits
   products:
     darkframe: 'master_dark.fits'
   logger:
     logfile: processing.log
     format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
     enabled: true
 
Generating template requirement files
-------------------------------------
Template requirement files can be generated by :program:`numina show-recipes`.

The flag generates templates for the named recipe or for all the available
recipes if no name is passed. 

For example::

  $ numina show-recipes -t megaradrp.recipes.DarkRecipe
  # This is a numina 0.13dev template file
  # for recipe 'megaradrp.recipes.DarkRecipe'
  #
  requirements: {master_bias: master_bias.fits}
  #products:
  # qc: QC.UNKNOWN
  # darkframe: darkframe.fits
  #logger:
  # logfile: processing.log
  # format: "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
  # enabled: true
  
The # character is a comment, so every line starting with it can safely 
removed. The names of FITS files in requirements must be edited to point
to existing files.

---------------------------------
Running the pipeline 
---------------------------------

Numina copies the images (calibrations and raw data) from directory ``datadir``
to directory ``workdir``, where the processing happens. The result is stored in
directory ``resultsdir``. The default values are ``_data``, ``_work`` and ``_results``.
All these directories can be defined in the command line.

See :ref:`numina:cli` for a full description of the command line interface.

Following the example, we create a directory ``_data`` in our current directory and copy
there the raw frames from ``r0121.fits`` to ``r0132.fits``and the master bias ``master_bias-1.fits``.

The we run::

  $ numina run obsresult.yaml -r requirements.yaml
  INFO: Numina simple recipe runner version 0.13dev
  INFO: Loading observation result from 'obsrun.yaml'
  INFO: Identifier of the observation result: 1
  INFO: instrument name: MEGARA
  ...
  numina.recipes.megara INFO stacking 4 images using median
  numina.recipes.megara INFO bias reduction ended
  INFO: result: BiasRecipeResult(qc=Product(type=QualityControlProduct(), dest='qc'), biasframe=Product(type=MasterBias(), dest='biasframe'))
  INFO: storing result

We get information of what's going on through logging messages. In the the end, the result and log files are stored in ``_results``.
The working directory ``_work`` can be inspected too. 


.. _YAML: http://www.yaml.org