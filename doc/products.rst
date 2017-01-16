
*************
Data Products
*************

Each recipe of the MEGARA Pipeline produces a set of predefined results, known
as *data products*. In turn, the recipes may request different data products
as computing requirements, effectively chaining the recipes.

For example, the requirements of :class:`~megaradrp.recipes.calibration.flat.FiberFlatRecipe`
include a :class:`~megaradrp.types.MasterDark` object. This object is produced by
the recipe :class:`~megaradrp.recipes.calibration.dark.DarkRecipe`, which in turn requires a
:class:`~megaradrp.types.MasterBias` object.

.. toctree::

   keywords
   images
