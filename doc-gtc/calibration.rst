Calibration Modes
===================

Bias
-----

:Mode: Bias
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
:Recipe class: :class:`~megaradrp.recipes.calibration.TraceMapRecipe`

Requirements
++++++++++++

Procedure
+++++++++

Products
++++++++
