########################
Validation of the inputs
########################

Before running a recipe, the pipeline can check that its inputs are valid:
that the raw images have the headers that MEGARA images must have, and
that they are the right kind of image for the observing mode (for example,
that the images of a bias have no exposure time, or that the lamps of an
arc are on). Raw images from GTC sometimes arrive with incomplete headers,
and the validation stops the reduction with a clear error instead of
failing later or producing wrong results.

The validation is disabled by default.

*************************
Enabling the validation
*************************

Use the option ``--validate`` of ``numina run``::

    numina run --validate obsresult.yaml -r control.yaml

or enable it in the configuration of numina, in the file
``.numina.cfg`` of the working directory, or in
``~/.config/numina/numina.cfg``:

.. code-block:: ini

    [tool.run]
    validate: True

The same value is used by :class:`numina.user.session.Session`, unless the
arguments ``validate_inputs`` or ``validate_results`` of
:meth:`~numina.user.session.Session.run` are given.

With the validation enabled, the results of the recipes are also
validated, after the recipe runs.

***************
What is checked
***************

Raw images
==========

Each raw image of the observation result is checked as the type of raw
image of the observing mode:

.. list-table::
   :header-rows: 1

   * - Observing modes
     - Type of raw image
     - Additional checks
   * - MegaraBiasImage
     - ``IMAGE_BIAS``
     - one HDU, ``EXPTIME`` and ``DARKTIME`` are 0, ``OBSMODE`` is
       ``MegaraBiasImage``
   * - MegaraDarkImage
     - ``IMAGE_DARK``
     - one HDU, ``OBSMODE`` is ``MegaraDarkImage``
   * - MegaraTraceMap, MegaraModelMap, MegaraFiberFlatImage
     - ``IMAGE_FLAT``
     - two HDUs, some incandescent lamp is on, all the comparison lamps
       are off
   * - MegaraArcCalibration
     - ``IMAGE_COMP``
     - two HDUs, all the incandescent lamps are off, some comparison lamp
       is on
   * - MegaraLcbImage, MegaraMosImage, MegaraLcbStdStar,
       MegaraMosStdStar, MegaraExtinctionStar, MegaraSensitivityStar,
       MegaraCrDetection
     - ``IMAGE_TARGET``
     - two HDUs, all the lamps are off
   * - MegaraSlitFlat
     - ``IMAGE_SLITFLAT``
     -
   * - MegaraTwilightFlatImage
     - ``IMAGE_TWILIGHT``
     -
   * - MegaraSuccess, MegaraFail, MegaraBadPixelMask,
       MegaraFocusSpectrograph, MegaraFocusTelescope,
       MegaraLcbAcquisition, MegaraMosAcquisition
     - ``IMAGE_TEST``
     -

In every case, the headers of the image are validated with the JSON
schema ``baseimage.json``, in the package :mod:`megaradrp.schemas`. The
schema describes the keywords of the primary header (``OBSMODE``,
``DATE-OBS``, ``INSTRUME``, ``INSCONF``, ``EXPTIME`` and ``VPH`` are
required), the size of the raw images (4196 x 4212 pixels) and the
``FIBERS`` extension. The ``FIBERS`` extension must have the keywords of
every bundle and fiber of the instrument mode (``INSMODE``), that must be
the same as in the primary header.

Some images of the commissioning of MEGARA, in 2017, do not follow the
current headers: they lack ``UUID``, or have an integer ``COVER``. The
schema accepts them for now, see
`issue #340 <https://github.com/guaix-ucm/megaradrp/issues/340>`_.

Observation results
===================

After the raw images, some observing modes check the observation result as
a whole:

* MegaraArcCalibration: all the images have the same ``SPECLAMP``.
* MegaraFocusSpectrograph: the images have a numeric ``FOCUS``, with at
  least two different values.

Calibrations
============

The calibrations given as requirements are checked as their type: a master
bias must be a processed image of 4096 x 4112 pixels with
``NUMTYPE = MasterBias``, a master fiber flat must be an RSS with three
HDUs, and so on. The structured calibrations
(trace maps, model maps and wavelength calibrations) are validated with
the schema ``basestruct.json``.

*****************************
When the validation fails
*****************************

If an input is not valid, the reduction stops before running the recipe.
The error lists the invalid inputs and, for each raw image, the name of
the file, the HDU and the keyword:

.. code-block:: text

    numina.exceptions.ValidationError: invalid BiasRecipeInput: obresult:
    invalid raw images: bad_bias.fits: ValidationError: HDU 0: 'VPH' is
    a required property

Other examples of errors:

.. code-block:: text

    image has 1 HDUs, not the expected number of HDUs (2)
    HDU 0, keyword EXPTIME: 10.0 is greater than the maximum of 0
    HDU 1: keyword BUN093_P is missing
    some incandescent lamps are ON

**************************
Checking images in Python
**************************

An image can be checked without running a recipe, with the function that
MEGARA registers in numina:

.. code-block:: python

    import astropy.io.fits as fits
    import numina.core.config as cfg
    import numina.drps
    from megaradrp.datatype import MegaraDataType

    numina.drps.get_system_drps()  # loads the DRPs and their checkers

    with fits.open("0001310880-20170827-MEGARA-MegaraBiasImage.fits.gz") as hdulist:
        cfg.check("MEGARA", hdulist, astype=MegaraDataType.IMAGE_BIAS)

The function returns ``True`` if the image is valid, and raises
:class:`numina.exceptions.ValidationError` otherwise. Without ``astype``,
the type of the image is inferred from its headers and shape, see
:func:`megaradrp.datamodel.megara_inferr_datatype`.

The details of the checks are in :mod:`megaradrp.validators`.
