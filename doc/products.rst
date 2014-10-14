
*********************
MEGARA Data Products
*********************

Each recipe of the MEGARA Pipeline produces a set of predefined results, known
as *data products*. In turn, the recipes may request diferent data products
as computing requeriments, efectively chaining the recipes.

For example, the requirements of the :ref:`ff-recipe-label` recipe include a
:class:`~megara.drp.dataproducts.MasterDark` object. This object is produced by 
the recipe :class:`~megara.drp.recipes.DarkRecipe`, which in turn requires a 
:class:`~megara.drp.dataproducts.MasterBias` object.

.. toctree::

   keywords
   rawimages
   images
