

MEGARA Reduction Recipes
=========================

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

Users of the DRP may use the simple Numina CLI (Command Line Interface) or the
higher level, database-driven Pontifex. Users of the DFP shall interact with
the software through the GTC Inspector. 

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

-----

:class:  ``focus.TelescopeFineFocusRecipe``  
:mode:  TS fine focus 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Parameter: objects 
:provides:  ``TelescopeFocus`` 

-----

:class:  ``focus.DTUFocusRecipe``  
:mode:  MEGARA focus control 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Parameter: objects 
    -  Parameter: msm_pattern 
    -  Parameter: dtu_focus_range 
:provides:  ``DTUFocus`` 

-----

:class:  ``acquisition.MaskCheckRecipe`` 
:mode:  Target acquisition 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
:provides: ``TelescopeOffset``

-----

:class:  ``acquisition.MaskImagingRecipe``       
:mode:  Mask image 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
:provides:  ``MSMPositions`` 

-----

:class:  ``acquisition.MaskCheckRecipe``  
:mode:  MSM and LSM check 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
:provides: ``TelescopeOffset``, ``MSMPositions`` 

.. raw:: pdf

    PageBreak

DRP Recipes Parameters
++++++++++++++++++++++

:class: ``auxiliary.BiasRecipe``
:mode: Bias image 
:requires:
:provides: ``MasterBias`` 

------------

:class: ``auxiliary.DarkRecipe``
:mode: Dark image 
:requires: Product: ``MasterBias``
:provides: ``MasterDark`` 

------------

:class: ``auxiliary.IntensityFlatRecipe``
:mode:  Intensity flat-field
:requires: 
        - Product: ``MasterBias``
        - Product: ``MasterDark``
        - Product: ``MasterBadPixelMask``
        - Product: ``NonLinearityCorrection``
:provides:  ``MasterIntensityFlat``

------------

:class:  ``auxiliary.SpectralFlatRecipe``  
:mode:  MSM spectral flat-field 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
:provides:  ``MasterSpectralFlat`` 

-----

:class:  ``auxiliary.SlitTransmissionRecipe``  
:mode:  Slit transmission calibration 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
:provides:  ``SlitTransmissionCalibration`` 

-----

:class:  ``auxiliary.WavelengthCalibrationRecipe``  
:mode:  Wavelength calibration 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Product: ``MasterSpectralFlatField``  
    -  Parameter: line_table (with wavelengths of arc lines)
:provides:  ``WavelengthCalibration`` 

-----

:class:  ``image.StareImageRecipe`` 
:mode:  Stare image 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  OptionalParameter: sources (list of sources coordinates)
:provides: ``Image``, ``SourcesCatalog``

-----

:class:  ``image.NBImageRecipe`` 
:mode:  Nodded/Beamswitched images 
:requires:
    - Product: ``MasterBias`` 
    - Product: ``MasterDark``  
    - Product: ``MasterBadPixelMask`` 
    - Product: ``NonLinearityCorrection`` 
    - Product: ``MasterIntensityFlatField``
    - Parameter: extinction (Mean atmospheric extinction)
    - Parameter: iterations 
    - Parameter: sky_images (Images used to estimate the background before and after current image)
    - Parameter: sky_images_sep_time (Maximum separation time between consecutive sky images in minutes)
    - Parameter: check_photometry_levels (Levels to check the flux of the objects)
    - Parameter: check_photometry_actions (Actions to take on images)
    - OptionalParameter: offsets (list of integer offsets between images)
:provides: ``Image``, ``SourcesCatalog``

-----

:class:  ``image.DitheredImageRecipe`` 
:mode:  Dithered images 
:requires:
    - Product: ``MasterBias`` 
    - Product: ``MasterDark``  
    - Product: ``MasterBadPixelMask`` 
    - Product: ``NonLinearityCorrection`` 
    - Product: ``MasterIntensityFlatField`` 
    - Parameter: extinction (Mean atmospheric extinction)
    - Parameter: iterations 
    - Parameter: sky_images (Images used to estimate the background before and after current image)
    - Parameter: sky_images_sep_time (Maximum separation time between consecutive sky images in minutes)
    - Parameter: check_photometry_levels (Levels to check the flux of the objects)
    - Parameter: check_photometry_actions (Actions to take on images)
:provides: ``Image``, ``SourcesCatalog``

-----

:class:  ``image.MicroditheredImageRecipe`` 
:mode:  Micro-dithered images 
:requires:
    -  *All the parameters of* ``image.DitheredImageRecipe``
    -  Parameter: subpixelization (number of subdivisions in each pixel side)
:provides: ``Image``, ``SourcesCatalog``

-----

:class:  ``image.MosaicRecipe`` 
:mode:  Mosaiced images 
:requires:  
:provides: ``Image``, ``SourcesCatalog``

-----

:class:  ``mos.StareSpectraRecipe`` 
:mode:  Stare spectra 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Product: ``MasterSpectralFlatField`` 
    -  Product: ``SlitTransmissionCalibration`` 
    -  Product: ``WavelengthCalibration``
    -  Parameter: lines (wavelength to measure)
:provides: ``Spectra``, ``LinesCatalog``

-----

:class:  ``mos.DNSpectraRecipe`` 
:mode:  Dithered/Nodded spectra along the slit
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Product: ``MasterSpectralFlatField``
    -  Product: ``SlitTransmissionCalibration`` 
    -  Product: ``WavelengthCalibration`` 
    -  Parameter: lines (wavelegnth to measure)
    -  OptionalParameter: offsets (list of integer offsets between images)
:provides: ``Spectra``, ``LinesCatalog``

-----

:class:  ``mos.OffsetSpectraRecipe`` 
:mode:  Offset spectra beyond the slit 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Product: ``MasterSpectralFlatField``
    -  Product: ``SlitTransmissionCalibration``
    -  Product: ``WavelengthCalibration``
    -  Parameter: lines (wavelegnth to measure)
    -  OptionalParameter: offsets (list of integer offsets between images)
:provides: ``Spectra``, ``LinesCatalog``

-----

:class:  ``mos.RasterSpectraRecipe`` 
:mode:  Raster spectra 
:requires:
    -  Product: ``MasterBias`` 
    -  Product: ``MasterDark``  
    -  Product: ``MasterBadPixelMask`` 
    -  Product: ``NonLinearityCorrection`` 
    -  Product: ``MasterIntensityFlatField`` 
    -  Product: ``MasterSpectralFlatField`` 
    -  Product: ``SlitTransmissionCalibration`` 
    -  Product: ``WavelengthCalibration`` 
    -  Parameter: lines (wavelegnth to measure)
:provides: ``DataCube``

-----

:class:  ``engineering.DTU_XY_CalibrationRecipe`` 
:mode:  DTU X_Y calibration 
:requires:
    -  Parameter: slit_pattern 
    -  Parameter: dtu_range 
:provides: ``DTU_XY_Calibration``

-----

:class:  ``engineering.DTU_Z_CalibrationRecipe`` 
:mode:  DTU Z calibration 
:requires:  Parameter: dtu_range 
:provides: ``DTU_Z_Calibration``

-----

:class: ``engineering.DTUFlexureRecipe`` 
:mode:  DTU Flexure compensation 
:requires:  
:provides: ``DTUFlexureCalibration``

-----

:class:  ``engineering.CSU2DetectorRecipe`` 
:mode:  CSU2Detector calibration 
:requires:  Parameter: dtu_range 
:provides: ``DTU_XY_Calibration``

-----

:class:  ``engineering.FocalPlaneCalibrationRecipe`` 
:mode:  Lateral colour 
:requires:  
:provides: ``PointingOriginCalibration``

-----

:class:  ``engineering.SpectralCharacterizationRecipe`` 
:mode:  Spectral characterization 
:requires:  
:provides: ``WavelengthCalibration``

-----

:class:  ``engineering.RotationCenterRecipe`` 
:mode:  Centre of rotation 
:requires:  
:provides: ``PointingOriginCalibration``

-----

:class:  ``engineering.AstrometricCalibrationRecipe`` 
:mode:  Astrometric calibration 
:requires:  
:provides: ``Image``

-----

:class:  ``engineering.PhotometricCalibrationRecipe`` 
:mode:  Photometric calibration 
:requires:  Parameter: phot 
:provides: ``PhotometricCalibration``

-----

:class:  ``engineering.SpectroPhotometricCalibrationRecipe`` 
:mode:  Spectrophotometric calibration 
:requires:  Parameter: sphot 
:provides: ``SpectroPhotometricCalibration``

