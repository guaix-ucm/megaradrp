Reduction Recipes
==================

Execution environment of the Recipes
------------------------------------

Recipes have different execution environments. Some recipes are designed to
process observing modes required while observing at the telescope. These modes
are related to visualization, acquisition and focusing. The corresponding
Recipes are integrated in the GTC environment. We call these recipes the **Data
Factory Pipeline**, (DFP).

Other group of recipes are devoted to scientific observing modes: imaging, spectroscopy and auxiliary calibrations. These Recipes constitute the
**Data Reduction Pipeline**, (DRP). The software is meant to be standalone,
users shall download the software and run it in their own computers, with
reduction parameters and calibrations provided by the instrument team.

Users of the DRP may use the simple Numina CLI (Command Line Interface). 
Users of the DFP shall interact with the software through the GTC Inspector. 

Recipe Parameters
-----------------
MEGARA Recipes based on Numina have a list of required parameters needed to properly configure the Recipe.
The Recipe announces the required parameters with the following syntax (the syntax is subject to changes).

.. code-block:: python

    class SomeRecipeInput(RecipeInput):
        master_dark = DataProductParameter(MasterDark, 'Master dark image') 
        some_numeric_value = Parameter(0.45, 'Some numeric value'),

    @define_input(SomeRecipeInput)
    class SomeRecipe(RecipeBase):        
        ...

When the reduction is run from the command line using Numina CLI, the program
checks that the required values are provided or have default values. When the
reduction is automatically executed using Pontifex, the program searches the
operation database looking for the most appropriated data products (in this
case, a MasterDark frame).

When the Recipe is properly configured, it is executed with an observing block
data structure as input. When run using Numina CLI, this data structure is
created from a text file. When run with Pontifex, the observing block data
structure is created from the contents of the database.

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


In the following two sections, we list the Reduction Recipes for the DRP and
for the DFP. The format is: name of the Python class of the recipe, name of the
observing mode, required parameters and data products provided. From the fully
quallified name of the recipe we have removed the initial ``megara.drp.recipes.``.

The name of the parameters are prefixed with **Product** if the parameter is
the result provided by another Recipe. If not, the value is a **Parameter**,
or an **OptionalParameter** that will be ignored if not present.

.. raw:: pdf

    PageBreak

DFP Recipes Parameters
++++++++++++++++++++++

:class:  ``focus.TelescopeRoughFocusRecipe``  
:mode:  TS rough focus 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Parameter: objects 
    -  Parameter: focus_range         
:provides:  ``TelescopeFocus`` 
