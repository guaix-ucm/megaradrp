Auxiliary Recipes
=================

Bias Image
----------

:Mode: Bias Image
:Recipe class: :class:`~megara.drp.recipes.BiasRecipe`
:Input class: :class:`~megara.drp.recipes.BiasRecipeInput`
:Result class: :class:`~megara.drp.recipes.BiasRecipeResult`

The actions to calibrate the zero (pedestal) level of the detector
plus associated control electronic by taken images with null
integration time.

Requeriments
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
| ``'biasframe'``   | :class:`~megara.drp.dataproducts.MasterBias`            |
+-------------------+---------------------------------------------------------+
| ``'stats'``       | :class:`~megara.drp.dataproducts.ChannelLevelStatistics`|
+-------------------+---------------------------------------------------------+

.. _ff-recipe-label:
