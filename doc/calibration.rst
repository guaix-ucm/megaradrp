Calibration Modes
=================

Arc
---

:Mode: Arc
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.ArcRecipe`

Requirements
++++++++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'obresult'``           | Product       | NA         |      Observation Result       |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``        | Product       | NA         |      Master Dark frame        |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         |      Master Bias frame        |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+
| ``'tracemap'``           | Product       | NA         |      TraceMap                 |
+--------------------------+---------------+------------+-------------------------------+
| ``'lines_catalog'``      | Product       | NA         |      Lines Catalog            |
+--------------------------+---------------+------------+-------------------------------+
| ``'polynomial_degree'``  | Product       | NA         |      Polynomial Degree        |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Arc.
This mode includes the required actions to translate the geometrical position
of each point in the detector into physical units of wavelength. The wavelength
calibration generated is used in other stages of the data processing.

Products
++++++++

+-------------------------+-------------------------------------------------+
| Name                    | Type                                            |
+=========================+=================================================+
| ``'arc_image'``         | :class:`~megaradrp.dataproducts.DataFrameType`  |
+-------------------------+-------------------------------------------------+
| ``'arc_rss'``           | :class:`~megaradrp.dataproducts.DataFrameType`  |
+-------------------------+-------------------------------------------------+
| ``'master_wlcalib'``    | :class:`~megaradrp.dataproducts.ArrayType`      |
+-------------------------+-------------------------------------------------+

A data structure containing information about wavelength calibrations
(the format is TBD), a QA flag, a text log file of the processing and a
structured text file containing information about the processing.

Bad-pixels mask
---------------

:Mode: Bad-pixels mask
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.BadPixelsMaskRecipe`


Requirements
++++++++++++

The observation result, that is: the series of images taken during the observing block and the
description of the observing block, including the parameters of the observing mode:

Procedure
+++++++++

The "User" processes an observing block obtained in the observing mode
Bad-pixels mask. This mode includes the required actions to obtain a bad-pixel
mask. The master bad pixel mask generated is used in other stages of the data
processing.

Products
++++++++

+-------------------+---------------------------------------------------------+
| Name              | Type                                                    |
+===================+=========================================================+
| ``'master_bpm'``  | :class:`~megaradrp.dataproducts.MasterBPM`              |
+-------------------+---------------------------------------------------------+

A bidimensional mask of bad pixels, a QA flag, a text log file of the
processing and a structured text file with information about the processing.





Bias
----

:Mode: Bias
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.BiasRecipe`

The actions to calibrate the zero (pedestal) level of the detector
plus associated control electronic by taken images with null
integration time.

Requirements
++++++++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The frames in the observed block are stacked together using the median of them as the final result.
The variance of the result frame is computed using two different methods.
The first method computes the variance across the pixels in the different frames stacked.
The second method computes the variance in each channel in the result frame.

Products
++++++++

+-------------------+---------------------------------------------------------+
| Name              | Type                                                    |
+===================+=========================================================+
| ``'master_bias'`` | :class:`~megaradrp.dataproducts.MasterBias`             |
+-------------------+---------------------------------------------------------+
| ``'stats'``       | :class:`~megaradrp.dataproducts.ChannelLevelStatistics` |
+-------------------+---------------------------------------------------------+

A bidimensional bias image, QA flag, a text log file of the processing and a
structured text file containing information about the processing.

Dark
----

:Mode: Dark
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.DarkRecipe`

Requirements
++++++++++++
+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'master_bias'``        | Product       | NA         |      Master Bias frame        |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Dark.
This mode includes the required actions to obtain a master dark frame. The
master dark generated is used in other stages of the data processing.

Products
++++++++
+------------------------------+-----------------------------------------------+
| Name                         | Type                                          |
+==============================+===============================================+
| ``'master_dark'``          | :class:`~megaradrp.dataproducts.MasterDark`     |
+------------------------------+-----------------------------------------------+

A bidimensional dark image, QA flag, a text log file of the processing and a
structured text file containing information about the processing.


Fiber-flat
----------

:Mode: Fiber-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.FiberFlatRecipe`

Requirements
++++++++++++
+---------------------------+---------------+------------+-------------------------------+
| Name                      | Type          | Default    | Meaning                       |
+===========================+===============+============+===============================+
| ``'obresult'``            | Product       | NA         |      Observation Result       |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``         | Product       | NA         |      Master Bias frame        |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``         | Product       | NA         |      Master Dark frame        |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_bpm'``          | Product       | NA         |      Master BPM frame         |
+---------------------------+---------------+------------+-------------------------------+
| ``'tracemap'``            | Product       | NA         |      TraceMap                 |
+---------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode
Fiber-flat. This mode includes the required actions to obtain a master
fiber-flat field. The master fiber-flat field generated is used in other stages
of the data processing.

Products
++++++++
+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
| ``'master_fiberflat_frame'`` | :class:`~megaradrp.dataproducts.MasterFiberFlatFrame` |
+------------------------------+-------------------------------------------------------+
| ``'master_fiberflat'``       | :class:`~megaradrp.dataproducts.MasterFiberFlat`      |
+------------------------------+-------------------------------------------------------+

A RSS master flat field, a QA flag, a text log file of the processing and a structured text file
containing information about the processing.


Slit-flat
---------

:Mode: Slit-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.SlitFlatRecipe`

Requeriments
++++++++++++
+---------------------------+---------------+------------+-------------------------------+
| Name                      | Type          | Default    | Meaning                       |
+===========================+===============+============+===============================+
| ``'obresult'``            | Product       | NA         |      Observation Result       |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``         | Product       | NA         |      Master Bias frame        |
+---------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``         | Product       | NA         |      Master Dark frame        |
+---------------------------+---------------+------------+-------------------------------+
| ``'window_length_x'``     | Product       | NA         |      Savitzky-Golay length    |
+---------------------------+---------------+------------+-------------------------------+
| ``'window_length_y'``     | Product       | NA         |      Savitzky-Golay length    |
+---------------------------+---------------+------------+-------------------------------+
| ``'polyorder'``           | Product       | NA         |      Savitzky-Golay order     |
+---------------------------+---------------+------------+-------------------------------+
| ``'median_window_length'``| Product       | NA         |      Median window width      |
+---------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode
Slit-flat. This mode includes the required actions to obtain a master slit-flat
field. The master slit-flat field generated is used in other stages of the data
processing.

Products
++++++++
+------------------------------+------------------------------------------------+
| Name                         | Type                                           |
+==============================+================================================+
| ``'master_slitflat'``        | :class:`~megaradrp.dataproducts.MasterSlitFlat`|
+------------------------------+------------------------------------------------+

A bidimensional master slit flat field, QA flag, a text log file of the
processing and a structured text file containing information about the
processing.


Trace
-----

:Mode: Trace
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.TraceMapRecipe`

Requirements
++++++++++++
+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'obresult'``           | Product       | NA         |      Observation Result       |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_dark'``        | Product       | NA         |      Master Dark frame        |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         |      Master Bias frame        |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bpm'``         | Product       | NA         |      Master BPM frame         |
+--------------------------+---------------+------------+-------------------------------+


Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Trace.
This mode includes the required actions to obtain a mapping of the trace of the
fibers. The master trace map generated is used in other stages of the data
processing.

Products
++++++++

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
| ``'master_fiberflat_frame'`` | :class:`~megaradrp.dataproducts.MasterFiberFlatFrame` |
+------------------------------+-------------------------------------------------------+
| ``'master_traces'``          | :class:`~megaradrp.dataproducts.TraceMap`             |
+------------------------------+-------------------------------------------------------+



Twilight fiber-flat
-------------------

:Mode: Twilight fiber-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.TwiligthFiberFlatRecipe`

Requeriments
++++++++++++

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'obresult'``           | Product       | NA         |      Observation Result       |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         |      Master Bias frame        |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++

The "User" processes an observing block obtained in the observing mode Twilight
Fiber Flat. This mode includes the required actions to obtain a master
illumination flat field. The master illumination flat field generated is used
in other stages of the data processing.

Products
++++++++

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
| ``'fiberflat_frame'``        | :class:`~megaradrp.dataproducts.MasterFiberFlatFrame` |
+------------------------------+-------------------------------------------------------+
| ``'fiberflat_rss'``          | :class:`~megaradrp.dataproducts.MasterFiberFlat`      |
+------------------------------+-------------------------------------------------------+
| ``'traces'``                 | :class:`~megaradrp.dataproducts.ArrayType`            |
+------------------------------+-------------------------------------------------------+

A RSS master illumination flat field, QA flag, a text log file of the
processing and a structured text file containing information about the
processing.














Standard star with the LCB IFU
------------------------------

:Mode: Standard start with the LCB IFU
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.LCB_IFU_StdStarRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++



Standard star with the Fiber MOS
--------------------------------

:Mode: Standard start with the FIBER MOS
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.FiberMOS_StdStarRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++



Linearity tests
---------------

:Mode: Linearity tests
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.LinearityTestRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++



