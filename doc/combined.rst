==============
Combined Modes
==============

The observing modes described in this section are those that require
the outputs of previous observing blocks.


Compute Extinction from Std Stars
---------------------------------

:Mode: Compute Extinction from Std Stars
:Usage: Offline
:Key: MEGARA_EXTINCTION_STAR
:Recipe: :class:`~megaradrp.recipes.combined.extinctionstar.Recipe`
:Recipe input: :class:`~megaradrp.recipes.combined.extinctionstar.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.combined.extinctionstar.RecipeResult`

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.combined.extinctionstar.Recipe
      :members:


Compute Extinction and Sensitivity from Std Stars
-------------------------------------------------


:Mode:
:Usage: Offline
:Key: MEGARA_SENSITIVITY_STAR
:Recipe: :class:`~megaradrp.recipes.combined.sensivitystar.Recipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.sensivitystar.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.sensivitystar.RecipeResult`

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.combined.sensivitystar.Recipe
      :members:
