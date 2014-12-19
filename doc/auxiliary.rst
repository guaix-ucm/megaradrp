Auxiliary Recipes
=================

Bias Image
----------

:Mode: Bias Image
:Recipe class: :class:`~megaradrp.recipes.BiasRecipe`
:Input class: :class:`~megaradrp.recipes.BiasRecipeInput`
:Result class: :class:`~megaradrp.recipes.BiasRecipeResult`

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
| ``'biasframe'``   | :class:`~megaradrp.dataproducts.MasterBias`            |
+-------------------+---------------------------------------------------------+
| ``'stats'``       | :class:`~megaradrp.dataproducts.ChannelLevelStatistics`|
+-------------------+---------------------------------------------------------+

.. _ff-recipe-label:
