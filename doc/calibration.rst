Calibration Modes
===================

Bias
-----

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
The second method computes the variance en each channel in the result frame.

Products
++++++++

+-------------------+---------------------------------------------------------+
| Name              | Type                                                    |
+===================+=========================================================+
| ``'biasframe'``   | :class:`~megaradrp.dataproducts.MasterBias`             |
+-------------------+---------------------------------------------------------+
| ``'stats'``       | :class:`~megaradrp.dataproducts.ChannelLevelStatistics` |
+-------------------+---------------------------------------------------------+



Dark
-----

:Mode: Dark
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.DarkRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++


Slit-flat
------------

:Mode: Slit-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.SlitFlatRecipe`

Requeriments
++++++++++++

Procedure
+++++++++

Products
++++++++

.. _ff-mode-label:

Fiber-flat
------------

:Mode: Fiber-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.FiberFlatRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++

Twilight fiber-flat
---------------------

:Mode: Twilight fiber-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.TwiligthFiberFlatRecipe`

Requeriments
++++++++++++

Procedure
+++++++++

Products
++++++++

Arc
------------

:Mode: Arc
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.ArcRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++

Standard star with the LCB IFU
---------------------------------

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
----------------------------------

:Mode: Standard start with the FIBER MOS
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.FiberMOS_StdStarRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++

Bad-pixels mask
----------------------------------

:Mode: Bad-pixels mask
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.BadPixelsMaskRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++

Linearity tests
----------------------------------

:Mode: Linearity tests
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.LinearityTestRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++


Trace
----------------------------------

:Mode: Trace
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.TraceMapRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++
