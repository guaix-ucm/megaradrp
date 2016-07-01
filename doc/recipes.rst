Reduction Recipes
=================

Execution environment of the Recipes
------------------------------------

Recipes have different execution environments. Some recipes are designed to
process observing modes required while observing at the telescope. These modes
are related to visualization, acquisition and focusing. The corresponding
Recipes are integrated in the GTC environment. We call these recipes the **Data
Factory Pipeline**, (DFP).

Other group of recipes are devoted to scientific observing modes and auxiliary calibrations. 
These Recipes constitute the **Data Reduction Pipeline**, (DRP). The software is meant to be standalone,
users shall download the software and run it in their own computers, with
reduction parameters and calibrations provided by the instrument team.

Users of the DRP may use the simple Numina CLI (Command Line Interface). 
Users of the DFP shall interact with the software through the GTC Inspector. 

Recipe Parameters
-----------------
MEGARA Recipes based on Numina have a list of required parameters needed to configure the Recipe properly.
The Recipe announces the required parameters with the following syntax (the syntax is subject to changes).

.. code-block:: python

    class SomeRecipeInput(RecipeInput):
        master_dark = DataProductParameter(MasterDark, 'Master dark image') 
        some_numeric_value = Parameter(0.45, 'Some numeric value'),

    @define_input(SomeRecipeInput)
    class SomeRecipe(RecipeBase):        
        ...


When the Recipe is configured properly, it is executed with an observing block
data structure as input. When is run using Numina CLI, this data structure is
created from an user-provided text file. The recipe requirements values are 
either provided in a text file or have default values.

Recipe Products
--------------- 
Recipes based on Numina provide a list of products created by the recipe.
The Recipe announces the required parameters with the following syntax
(the syntax is subject to changes).

.. code-block:: python

    class SomeRecipeInput(RecipeInput):
        master_dark = DataProductParameter(MasterDark, 'Master dark image') 
        some_numeric_value = Parameter(0.45, 'Some numeric value'),
        
    class SomeRecipeResult(RecipeResult):
        master_flat = Product(MasterDark) 
        
    @define_input(SomeRecipeInput)
    @define_result(SomeRecipeResult)
    class SomeRecipe(RecipeBase):        
        ...

        
The data products of the MEGARA DRP are describe in :doc:`products`