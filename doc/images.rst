
Image data products
===================

These data products are saved to disk as FITS files. MEGARA DRP makes use of the FITS headers
to record information about the data processing. This information may be recorded using other
methods as well, such as the GTC Database.

The following headers are included in all image data products and record information
about the version of Numina and the name and version of the recipe used.

  ::

   NUMXVER = '0.13.0   '           / Numina package version                         
   NUMRNAM = 'BiasRecipe'          / Numina recipe name                            
   NUMRVER = '0.1.0   '            / Numina recipe version                                     
   NUMTYP  = 'TARGET  '            / Data product type  

``HISTORY`` keywords may be used also, but the information in these keyword may not be easily indexed.

Master Bias frames
*******************

Bias frames are produced by the recipe :class:`~megaradrp.recipes.BiasRecipe`. Each bias frame is a 
multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The bias level
  ``VARIANCE``       Image                Variance of the bias level
  ``MAP``            Image                Number of pixels used to compute the bias level
  ===============    =======   ========   =======================

Master bias frames are represented by :class:`~megaradrp.dataproducts.MasterBias`.

Master Dark frames
******************

Master dark frames are produced by the recipe :class:`~megaradrp.recipes.DarkRecipe`. Each dark frame is a 
multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The dark level
  ``VARIANCE``       Image                Variance of the dark level
  ``MAP``            Image                Number of pixels used to compute the dark level
  ===============    =======   ========   =======================

Master dark frames are represented by :class:`~megaradrp.dataproducts.MasterDark`.


