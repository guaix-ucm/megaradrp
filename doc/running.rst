####################
Running the pipeline
####################

The MEGARA DRP is run through a command line interface
provided by :program:`numina`.

The run mode of numina requires:
 
  * A observation result file in YAML_ format
  * A requirements file in YAML format 
  * The raw images obtained in the observing block
  * The calibrations required by the recipe
 
The observation result file and the requirements file are created by the user,
the format is described in the following sections.
 
********************************
Format of the observation result
********************************

The contents of the file is a serialized dictionary with the
following keys:

*id*: not required, string, defaults to 1
    Unique identifier of the observing block

*instrument*: required, string
    Name of the instrument, as it is returned by ``numina show-instruments``

*mode*: required, string
    Name of the observing mode, as returned by ``numina show-modes``

*frames*: required, list of strings
    List of images names

*children*: not required, list of integers, defaults to empty list
    Identifications of nested observing blocks

This is an example of the observation result file

.. code-block:: yaml

    id: dark-test-21
    instrument: MEGARA
    mode: MegaraDarkImage
    images:
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
   
*******************************
Format of the requirements file
*******************************

This file contains calibrations obtained by running recipes (called **products**)
and other parameters (numeric or otherwise) required by the recipes (named **requirements**). The file
is serialized using YAML_


Example requirements file:

.. code-block:: yaml

    version: 1                                             (1)
    products:                                              (2)
      EMIR:
       - {id: 1, content: 'file1.fits', type: 'MasterFlat', tags: {'filter': 'J'}, ob: 200}     (3)
       - {id: 4, content: 'file4.fits', type: 'MasterBias', tags: {'readmode': 'cds'}, ob: 400} (3)
      MEGARA:
       - {id: 1, content: 'file1.fits', type: 'MasterFiberFlat', tags: {'vph': 'LR-U'}, ob: 1200} (3)
       - {id: 2, content: 'file2.yml', type: 'TraceMap', tags: {'vph': 'LR2', 'readmode': 'fast'}, ob: 1203} (3)
    requirements: (4)
      MEGARA:
         default:
           MegaraArcImage:  (5)
              polynomial_degree: 5 (6)
              nlines: [5, 5]       (6)


1. Mandatory entry, ``version`` must be 1
2. Products of other recipes are list, by instrument
3. The products of the reduction recipes are listed. Each result must contain:
    * A ``type``, one of the types of the products of the DRP in string format
    * A ``tags`` field, used to select the correct calibration based on the keywords of
      the input.
    * A ``content`` field, a pointer to the serialized version of the calibration.
    * A ``id`` field, unique integer
    * A ``ob`` field, optional integer, used to store the observation id of the images that
      created the calibration.
4. Numerical parameters of the recipes are stored in ``requirements``, with different sections
   per instrument.
5. The name of the observing mode.
6. Different parameters for the recipe corresponding to the observing mode in (5)


********************
Running the pipeline 
********************

:program:`numina` copies the images (calibrations and raw data) from directory 
``datadir`` to directory ``workdir``, where the processing happens. 
The result is stored in directory ``resultsdir``. 
The default values are for each directory are ``data``, ``obsid<id_of_obs>_work`` and ``obsid<id_of_obs>_results``.
All these directories can be defined in the command line using flags::

  $ numina run --workdir /tmp/test1 --datadir /scrat/obs/run12222 obs.yaml -r requires.yaml

See :ref:`numina:cli` for a full description of the command line interface.

Following the example, we create a directory ``data`` in our current directory and copy
there the raw frames from ``r0121.fits`` to ``r0132.fits`` and the master bias ``master_bias-1.fits``.

The we run::

  $ numina run obsresult.yaml -r requirements.yaml
  INFO: Numina simple recipe runner version 0.15
  INFO: Loading observation result from 'obsrun.yaml'
  INFO: Identifier of the observation result: 1
  INFO: instrument name: MEGARA
  ...
  numina.recipes.megara INFO stacking 4 images using median
  numina.recipes.megara INFO bias reduction ended
  INFO: result: BiasRecipeResult(qc=Product(type=QualityControlProduct(), dest='qc'), biasframe=Product(type=MasterBias(), dest='biasframe'))
  INFO: storing result

We get information of what's going on through logging messages. In the end, the result and log files are stored in ``obsid<id_of_obs>_results``.
The working directory ``obsid<id_of_obs>_work`` can be inspected too. Intermediate results will be saved here.


On the other hand, in the following we attach a short code to run megaradrp
by using a Python script. This is useful to use the Python debugger.

.. code-block:: python

    from numina.user.cli import main
    from megaradrp.loader import load_drp

    def run_recipe():
        main(['run', 'obsresult.yaml', '-r', 'requirements.yaml'])

    if __name__ == "__main__":
        run_recipe()


Pipeline's Flow Example
-----------------------
In this subsection, we detail an example about how to generate a called Master
Fiber Flat Image. To achieve our goal, a schematic flow can be seen in the next
Figure:

.. graphviz::

    digraph G {
        rankdir=LR;
        subgraph cluster_0 {
            style=filled;
            color=lightgrey;
            node [style=filled,color=white];
            edge[style=invis]
            a0 -> a5,a1 -> a4,a2 -> a3;
            #label = "Observing\nModes";
        }

        a0 -> a1 [rank=same];
        a1 -> a2 [rank=same];
        a1 -> a4 [rank=same];
        a2 -> a3 [rank=same];
        a4 -> a3 [rank=same];
        a5 -> a4 [rank=same];

        a0 [label="MegaraBiasImage"];
        a1 [label="MegaraTraceMap"];
        a2 [label="MegaraArcCalibration"];
        a3 [label="MegaraModelMap"];
        a4 [label="MegaraFiberFlat"];
        a5 [label="MegaraSlitFlat"];

    }

It is important to emphasize the fact that each time a Recipe is run, the results must be
renamed and copied to the ``data`` directory in order to be the input of the
next Recipe if it is needed. Taking this in mind, the content of the
``requirements.yaml`` file might well be and is common to all Recipes:

.. code-block:: yaml

    version: 1
    products:
      MEGARA:
      - {id: 1, type: 'LinesCatalog', tags: {}, content: 'ThAr_arc_LR-U.txt'}
      - {id: 2, type: 'MasterBias', tags: {}, content: 'master_bias.fits'}
      - {id: 3, type: 'TraceMap', tags: {}, content: 'master_traces.json'}
      - {id: 4, type: 'MasterFiberFlat', tags: {}, content: 'master_fiberflat.fits'}
      - {id: 5, type: 'WavelengthCalibration', tags: {}, content: 'master_wlcalib.json'}
      - {id: 6, type: 'MasterFiberFlatFrame', tags: {}, content: 'fiberflat_frame.fits'}
      - {id: 7, type: 'ModelMap', tags: {}, content: 'master_model.json'}
      - {id: 8, type: 'MasterSlitFlat', tags: {}, content: 'master_slitflat.fits'}
    requirements: {}

In order to run the next example, the user should execute the next command
at least 6 times taking into account that the file ``obsresult-%step.yaml`` should
change with each execution::

    $ numina run obsresult-1.yaml -r requirements.yaml
    $ numina run obsresult-2.yaml -r requirements.yaml
    ...
    $ numina run obsresult-6.yaml -r requirements.yaml

MegaraBiasImage file, obsresult-1.yaml:

.. code-block:: yaml

    id: 1
    instrument: MEGARA
    mode: MegaraBiasImage
    images:
      - bias1.fits
      - bias2.fits
      - bias3.fits
      - bias4.fits
      - bias5.fits

MegaraTraceMap, obsresult-2.yaml:

.. code-block:: yaml

    id: 2
    instrument: MEGARA
    mode: MegaraTraceMap
    images:
      - flat1.fits
      - flat2.fits
      - flat3.fits
      - flat4.fits
      - flat5.fits

MegaraArcCalibration, obsresult-3.yaml:

.. code-block:: yaml

    id: 3
    instrument: MEGARA
    mode: MegaraArcCalibration
    images:
      - arc1.fits
      - arc2.fits
      - arc3.fits
      - arc4.fits
      - arc5.fits

MegaraSlitFlat, obsresult-4.yaml:

.. code-block:: yaml

    id: 4
    instrument: MEGARA
    mode: MegaraSlitFlat
    images:
      - flat1.fits
      - flat2.fits
      - flat3.fits
      - flat4.fits
      - flat5.fits

MegaraModel, obsresult-5.yaml:

.. code-block:: yaml

    id: 5
    instrument: MEGARA
    mode: MegaraModelMap
    images:
      - flat1.fits
      - flat2.fits
      - flat3.fits
      - flat4.fits
      - flat5.fits

MegaraFiberFlat, obsresult-6.yaml:

.. code-block:: yaml

    id: 6
    instrument: MEGARA
    mode: MegaraFiberFlat
    images:
      - flat1.fits
      - flat2.fits
      - flat3.fits
      - flat4.fits
      - flat5.fits

Notice that if you would want to execute this example automatically, you could
code a script (following the same skeleton as shown above) with a loop flow to
read the .yaml files and the outputs that each recipe generates.

.. _YAML: http://www.yaml.org
