==============
Combined Modes
==============

The observing modes described in this section are those that require
the outputs of previous observing blocks.


Compute Sensitivity from Std Stars
----------------------------------

:Mode:
:Usage: Offline
:Key: MegaraSensitivityStar
:Product: :class:`~megaradrp.types.MasterSensitivity`
:Recipe: :class:`~megaradrp.recipes.combined.sensstar.Recipe`
:Recipe input: :class:`~megaradrp.recipes.combined.sensstar.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.combined.sensstar.RecipeResult`

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.combined.sensstar.Recipe
      :members:



Compute Extinction and Sensitivity from Std Stars
-------------------------------------------------

:Mode: Compute Extinction from Std Stars
:Usage: Offline
:Key: MegaraExtinctionStar
:Product: :class:`~megaradrp.types.MasterSensitivity`, :class:`~megaradrp.types.Extinction`
:Recipe: :class:`~megaradrp.recipes.combined.extinctionstar.Recipe`
:Recipe input: :class:`~megaradrp.recipes.combined.extinctionstar.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.combined.extinctionstar.RecipeResult`

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.combined.extinctionstar.Recipe
      :members:


