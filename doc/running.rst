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

*id*: required, string or integer
    Unique identifier of the observing block

*instrument*: required, string
    Name of the instrument, as it is returned by ``numina show-instruments``

*mode*: required, string
    Name of the observing mode, as returned by ``numina show-modes``

*frames* (or *images*): list of strings
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

    version: 1                                                       # (1)
    products:                                                        # (2)
      MEGARA:
        4fd05b24-2ed9-457b-b563-a3c618bb1d4c:                        # (3)
          - {id: 1, content: 'file1.fits', type: 'MasterFiberFlat', tags: {'vph': 'LR-U'}}  # (4)
          - {id: 2, content: 'file2.json', type: 'TraceMap', tags: {'vph': 'LR-R', 'insmode': 'LCB'}}  # (4)
    requirements:                                                    # (5)
      MEGARA:
        4fd05b24-2ed9-457b-b563-a3c618bb1d4c:                        # (3)
          default:                                                   # (6)
            MegaraArcCalibration:                                    # (7)
              - {name: polynomial_degree, tags: {}, content: 5}      # (8)
              - {name: nlines, tags: {'vph': 'LR-U'}, content: [5, 5]}  # (8)


1. Mandatory entry, ``version`` must be 1
2. Products of other recipes, by instrument and instrument profile
3. The instrument profile, the UUID of the instrument configuration
   (keyword ``INSCONF`` of the images)
4. The products of the reduction recipes. Each product contains:
    * A ``type``, one of the types of the products of the DRP in string format
    * A ``tags`` field, used to select the correct calibration based on the keywords of
      the input.
    * A ``content`` field, the path of the calibration, relative to the data directory.
    * A ``id`` field, unique integer
5. Numerical parameters of the recipes are stored in ``requirements``, by instrument,
   instrument profile, pipeline and observing mode.
6. The name of the pipeline.
7. The name of the observing mode.
8. The parameters of the recipe of the observing mode in (7). Each parameter contains
   its ``name``, its value (``content``) and the ``tags`` of the observations where it
   is used.

The file is validated with the schema ``control-schema.json`` of numina
(:mod:`numina.schemas`).


********************
Running the pipeline
********************

:program:`numina` links the images (calibrations and raw data) from directory
``datadir`` in directory ``workdir``, where the processing happens (with
``--copy-files``, the images are copied instead).
The result is stored in directory ``resultsdir``.
The default values of these directories are ``data``, ``obsid<id_of_obs>_work`` and ``obsid<id_of_obs>_results``,
inside ``basedir``, that is the current directory by default.
The names of the work and result directories are defined in the configuration of numina
(``workdir_tmpl`` and ``resultdir_tmpl`` in the section ``[tool.run]``).
The base and data directories can be defined in the command line using flags::

  $ numina run --basedir /tmp/test1 --datadir /scrat/obs/run12222 obs.yaml -r requires.yaml

See :ref:`numina:cli` for a full description of the command line interface.

Following the example, we create a directory ``data`` in our current directory and copy
there the raw frames from ``r0121.fits`` to ``r0132.fits`` and the master bias ``master_bias-1.fits``.

Then we run::

  $ numina run obsresult.yaml -r requirements.yaml
  INFO: Numina simple recipe runner version 0.40.0
  INFO: observation results from 'obsresult.yaml'
  INFO: reading control from requirements.yaml
  INFO: control format version 1
  INFO: procesing OB with id=dark-test-21
  INFO: running recipe
  megaradrp.core.correctors INFO loading master_bias
  ...
  megaradrp.processing.trimover DEBUG trimming image 2017-11-09T11:11:00.0
  INFO: stacking 12 images using 'median'
  INFO: storing result uuid=7d495e2e-c7ab-11f1-963b-a0369f681a26, quality=QC.GOOD
  INFO: storing task id=dark-test-21

The requirements file of this example is:

.. code-block:: yaml

    version: 1
    products:
      MEGARA:
        4fd05b24-2ed9-457b-b563-a3c618bb1d4c:
          - {id: 1, content: 'master_bias-1.fits', type: 'MasterBias', tags: {}}

We get information of what's going on through logging messages. In the end, the result and log files are stored in ``obsid<id_of_obs>_results``
(``master_dark.fits``, ``result.json``, ``task.json`` and ``processing.log``).
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
        4fd05b24-2ed9-457b-b563-a3c618bb1d4c:
          - {id: 1, type: 'LinesCatalog', tags: {}, content: 'ThAr_arc_LR-U.txt'}
          - {id: 2, type: 'MasterBias', tags: {}, content: 'master_bias.fits'}
          - {id: 3, type: 'TraceMap', tags: {}, content: 'master_traces.json'}
          - {id: 4, type: 'MasterFiberFlat', tags: {}, content: 'master_fiberflat.fits'}
          - {id: 5, type: 'WavelengthCalibration', tags: {}, content: 'master_wlcalib.json'}
          - {id: 6, type: 'ModelMap', tags: {}, content: 'master_model.json'}
          - {id: 7, type: 'MasterSlitFlat', tags: {}, content: 'master_slitflat.fits'}
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
