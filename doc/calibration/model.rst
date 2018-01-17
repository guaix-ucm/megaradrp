=========
Model Map
=========

:Mode: ModelMap
:Usage: Offline
:Key: MegaraModelMap
:Product: :class:`~megaradrp.products.modelmap.ModelMap`.
:Recipe: :class:`~megaradrp.recipes.calibration.modelmap.ModelMapRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.modelmap.ModelMapRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.modelmap.ModelMapRecipeResult`


Although for the majority of the observing modes described elsewhere in this
document the MEGARA off-line pipeline will perform its own fiber spectra
extraction from the 2D CCD FITS frame, there are cases where an archival master
"model map" should be used instead. Note that a different "model map" should be
available for each pseudo-slit and VPH combination.

Requirements
++++++++++++
This offline observing mode will use the images obtained in the online observing mode
"Trace".

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Trace.
This mode includes the required actions to obtain a mapping of the profiles of the
fibers. The master model map generated is used in other stages of the data
processing.

Products
++++++++

Trace map image sets are to be obtained both as part of the activities related
to the verification of the instrument status and for processing data for
scientific exploitation. Note, however, that the use of this observing mode for
scientific exploitation should be limited as it could affect to the general
performance of the on-line quick-look software.

This mode produces the profile information required to perform advanced extraction of
the fibers.
The result is stored in an object named  ``master_model``
of type :class:`~megaradrp.products.modelmap.ModelMap`.


Recipe
------

.. autoclass:: megaradrp.recipes.calibration.modelmap.ModelMapRecipe
   :members:
