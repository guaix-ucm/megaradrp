=============
Data products
=============

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

*************
Generic types
*************


Processed Frame
***************

Processed Frame is the type of any image produced by the pipeline that represents a view of the detector.
Its size may be 4096x4112 for trimmed images or 4196x4212 for unprocessed, raw images.

Processed Frame is represented by :class:`~megaradrp.types.ProcessedFrame`.

Processed RSS
*************
Processed Row Stacked Spectra is the type of any image produced by the pipeline that represents
a view of the focal plane, using extracted fibers. It will have 623 rows in LCB mode and 644 rows
in MOS mode, each row representing the extracted spectrum of one fiber. The number of columns
will be 4112 or larger, depending on the stage of the reduction.

Procesed RSS images will tipically have a `FIBERS` extension.

Processed RSS is represented by :class:`~megaradrp.types.ProcessedRSS`.


Processed Spectrum
******************
Processed spectrum is the type of any image produced by the pipeline that represents
the spectrum of one object. Its data content will be a 1D array.

Processed Spectrum is represented by :class:`~megaradrp.types.ProcessedSpectrum`.


************
Calibrations
************

Master Bias frames
******************

Bias frames are produced by the recipe :class:`~megaradrp.recipes.calibration.bias.BiasRecipe`. Each bias frame is a
multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The bias level
  ``VARIANCE``       Image                Variance of the bias level
  ``MAP``            Image                Number of pixels used to compute the bias level
  ===============    =======   ========   =======================

Master bias frames are represented by :class:`~megaradrp.types.MasterBias`.

Master Dark frames
******************

Master dark frames are produced by the recipe :class:`~megaradrp.recipes.calibration.dark.DarkRecipe`. Each dark frame is a
multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The dark level
  ``VARIANCE``       Image                Variance of the dark level
  ``MAP``            Image                Number of pixels used to compute the dark level
  ===============    =======   ========   =======================

Master dark frames are represented by :class:`~megaradrp.types.MasterDark`.


Master Bad Pixel Mask
*********************

Master Bad Pixel Mask is produced by the recipe :class:`~megaradrp.recipes.calibration.bpm.BadPixelsMaskRecipe`.
Each bad pixel mask frame is a multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The Bad Pixel Mask level
  ===============    =======   ========   =======================

Master bad pixel mask frames are represented by :class:`~megaradrp.types.MasterBPM`.

Master Slit Flat
****************

Master Slit Flat is produced by the recipe :class:`~megaradrp.recipes.calibration.slitflat.SlitFlatRecipe`.
Each slit flat frame is a multiextension FITS file with the following extensions.

  ===============    =======   =======================
  Extension name     Type      Contents
  ===============    =======   =======================
  ``PRIMARY``        Primary   The Slit Flat level
  ===============    =======   =======================

Masterslit flat frames are represented by :class:`~megaradrp.types.MasterSlitFlat`.

Master Traces
*************

Master Fiber Flat is produced by the recipe :class:`~megaradrp.recipes.calibration.trace.TraceMapRecipe`.
The result is a JSON_ file where each one of the records belongs to a given fiber
in the RSS file. Moreover, each one of the records has the next information:


  ===============    =======   =======================
  Field              Type      Contents
  ===============    =======   =======================
  ``boxid``          Integer   Number of the box
  ``fibid``          Integer   Number of the fiber
  ``fitparms``       Primary   Polyfit algorithm result
  ``start``          Integer   X-Coordenate in the Flat image
  ``stop``           Integer   X-Coordenate in the Flat image
  ===============    =======   =======================


In the following, a real example of the fourth fiber which is in the first box  can be seen in the yaml format:

.. code-block:: yaml

  - boxid: 1
    fibid: 4
    fitparms: [2.6909627476636523e-18, -3.0949058966515047e-14, 1.872326137294402e-10,1.1602592442769502e-06, -0.0009443161994027746, 262.01840282676613]
    start: 4
    stop: 3594

Master Tracemap files are represented by :class:`~megaradrp.products.tracemap.TraceMap`.

..
  Master Weights
  **************

  Master weights file is produced by the recipe :class:`~megaradrp.recipes.calibration.weights.WeightsRecipe`.
  This is a .tar file which is made up of 4096 .npz files (one per fiber). These are
  ``numpy`` files where the ndarray are stored.

  This file is compulsary to generate the master fiber flat.

  Master weights files are represented by :class:`~megaradrp.products.MasterWeights`.


Master Wavelength Calibration
*****************************

Master wavelength calibration is produced by the recipe :class:`~megaradrp.recipes.calibration.arc.ArcCalibrationRecipe`.
The result is a JSON_ file where each one of the records belongs to a given fiber
in the RSS file. Moreover, each one of the records or ``apertures`` has the next fields:

  ============    ==========    =======================
  Field           Type          Contents
  ============    ==========    =======================
  ``features``    List          List with the arc's information
  ``function``    Dictionary    Number of pixels used to compute the dark level
  ``id``          Integer       Number the corresponding fiber
  ============    ==========    =======================

Additionally, each one of the elements that belongs to the ``features``
corresponds to each one of the arc lines that has been found in the RSS image.
The dictionary that each element has, contains the next information:

  ===============    =======     =======================
  Field              Type        Contents
  ===============    =======     =======================
  ``category``       String      Type of the arc
  ``flux``           Float       Flux of the arc
  ``fwhm``           Float       Full Width at Half Maximum of the arc
  ``reference``      Float       Line in the Catalog lines
  ``wavelength``     Float       Predicted line
  ``xpos``           Float       X-coordenate of the arc in the RSS image
  ``ypos``           Float       Y-coordenate of the arc in the RSS image
  ===============    =======     =======================

Finally, the ``function`` dictionary has three elements: ``coefficients``,
``method`` and ``order`` fields. Coefficients has the result of executing
the ``polynomial.polyfit`` numpy method. Method field has the name of the
algorithm used. Order field has the polynomial degree.

In the following, an example of the first fiber of a real JSON file with only
two arc lines can be seen:

.. code-block:: json

  {
    "aperture": {
      "features": [
        {
          "category": "E",
          "flux": 50212.563405324945,
          "fwhm": 3.438967092459162,
          "reference": 6013.2816999999995,
          "wavelength": 6013.2847301957181,
          "xpos": 33.267395825699928,
          "ypos": 251.10097403866305
        },
      ],
      "function": {
        "coefficients": [6001.573165443434,0.35298729563735487,-2.898410563853586e-05,1.858317850662985e-08,-8.411429549924489e-12,1.4341696725726076e-15],
        "method": "least squares",
        "order": 5
      },
      "id": 2
    }


Master Wavelength calibration file is represented by :class:`~megaradrp.products.wavecalibration.WavelengthCalibration`.


Master Fiber Flat
*****************

Master Fiber Flat is produced by the recipe :class:`~megaradrp.recipes.calibration.flat.FiberFlatRecipe`.
Each master fiber flat frame is a multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The Fiber Flat level
  ``FIBERS``         Image                Description of the focal plane
  ===============    =======   ========   =======================

Master fiber flats frames are represented by :class:`~megaradrp.types.MasterFiberFlat`.

Master Twilight Flat
********************

Master Twilight Flat is produced by the recipe :class:`~megaradrp.recipes.calibration.twilight.TwilightFiberFlatRecipe`.
Each twilight flat frame is a multiextension FITS file with the following extensions.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The Twilight Flat level
  ``FIBERS``         Image                Description of the focal plane
  ===============    =======   ========   =======================

Master twilight flat frames are represented by :class:`~megaradrp.types.MasterTwilightFlat`.


Master Sensitivity
******************

Master sensitivity star image is produced by the recipe :class:`~megaradrp.recipes.composed.sensstar.Recipe`.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The Sensitivity Star Image level
  ===============    =======   ========   =======================

Master sensitivity star image is represented by :class:`~megaradrp.types.MasterSensitivity`.


Master Extinction
*****************

Master extinction star image is produced by the recipe :class:`~megaradrp.recipes.composed.extinctionstar.Recipe`.

  ===============    =======   ========   =======================
  Extension name     Type      Version    Contents
  ===============    =======   ========   =======================
  ``PRIMARY``        Primary              The Extinction Star Image level
  ===============    =======   ========   =======================

Master extinction star image is represented by :class:`~megaradrp.types.Extinction`.


**********************
Reference calibrations
**********************

The following types represent types used for calibration, but that are not the
result of any recipe. Examples of this type are the spectra of flux standars
or the tables of spectral lines of calibration lamps.


Reference Spectrum
******************

A tabular representation of the spectral energy distribution of a standard star.
The first column contains wavelength (in Angstroms) and the second column
the flux in erg/s/cm^2/Angstrom

Reference spectrum is represented by :class:`~megaradrp.types.ReferenceSpectrum`.


.. _JSON: http://www.json.org/
