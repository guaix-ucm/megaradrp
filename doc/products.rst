
*********************
Data Products
*********************

Each recipe of the MEGARA Pipeline produces a set of predefined results, known
as *data products*. In turn, the recipes may request different data products
as computing requirements, effectively chaining the recipes.

For example, the requirements of the :ref:`ff-mode-label` recipe include a
:class:`~megara.drp.dataproducts.MasterDark` object. This object is produced by 
the recipe :class:`~megaradrp.recipes.calibration.DarkRecipe`, which in turn requires a 
:class:`~megara.drp.dataproducts.MasterBias` object.

.. toctree::

   keywords
   rawimages
   images
