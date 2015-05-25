Reduction Recipes
==================

Recipe Parameters
-----------------
MEGARA Recipes based on Numina have a list of required parameters needed
to properly configure the Recipe.
The Recipe announces the required parameters with the following
syntax (the syntax is subject to changes).

.. code-block:: python

    class SomeRecipeInput(RecipeInput):
        master_dark = DataProductParameter(MasterDark, 'Master dark image') 
        some_numeric_value = Parameter(0.45, 'Some numeric value'),

    @define_input(SomeRecipeInput)
    class SomeRecipe(RecipeBase):        
        ...


When the Recipe is properly configured, it is executed with RecipeInput 
data estrcuture as input.

Recipe Products
--------------- 
Recipes based on Numina provide a list of products created by the recipe.
The Recipe announces the provided outputs with the following syntax
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

        
The data products of the MEGARA DFP recipes are described in :doc:`products`

