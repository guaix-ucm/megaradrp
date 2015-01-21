Observing modes
==================

TBD

.. toctree::
   :maxdepth: 1

   auxiliary
   calibration
   scientific


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


