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

If an input is not valid, the reduction stops before running the recipe,
and ``numina run`` exits with status 1. The error lists the invalid inputs
and, for each raw image, the name of the file, the HDU and the keyword:

.. code-block:: text

    ERROR: invalid BiasRecipeInput: obresult: invalid raw images:
    bad_bias.fits: ValidationError: HDU 0: 'VPH' is a required property

With the option ``-d`` of numina (``numina -d run ...``), the traceback
of the error is also shown.

Other examples of errors:

.. code-block:: text

    image has 1 HDUs, not the expected number of HDUs (2)
    HDU 0, keyword EXPTIME: 10.0 is greater than the maximum of 0
    HDU 1: keyword BUN093_P is missing
    some incandescent lamps are ON

*********************************
Checking images before reducing
*********************************

The raw images can be checked without running a recipe, with
``numina verify``, for example to discard the images with incomplete
headers as soon as they arrive from GTC::

    $ numina verify 0001310880-20170827-MEGARA-MegaraBiasImage.fits.gz 0001310889-20170827-MEGARA-MegaraBiasImage.fits
    OK 0001310880-20170827-MEGARA-MegaraBiasImage.fits.gz
    INVALID 0001310889-20170827-MEGARA-MegaraBiasImage.fits: HDU 0: 'VPH' is a required property
    1 valid, 1 invalid, 0 not checked

The type of each image is inferred from its headers and shape (see
:func:`megaradrp.datamodel.megara_inferr_datatype`). With ``--mode``, the
images are checked as the raw images of an observing mode, as in the table
above::

    $ numina verify --mode MegaraArcCalibration 0001312246-20170831-MEGARA-MegaraSuccess.fits.gz
    INVALID 0001312246-20170831-MEGARA-MegaraSuccess.fits.gz: some incandescent lamps are ON
    0 valid, 1 invalid, 0 not checked

With ``--ob``, the files are observation results: their raw images, in the
directory given by ``--datadir``, are checked as the raw images of the
observing mode, and then the observation result with the checks of the
mode, as ``numina run --validate`` does::

    $ numina verify --ob --datadir data 3_M15_wavecalib.yaml
    OK data/0001312249-20170831-MEGARA-MegaraSuccess.fits.gz
    OK data/0001312250-20170831-MEGARA-MegaraSuccess.fits.gz
    OK data/0001312251-20170831-MEGARA-MegaraSuccess.fits.gz
    OK 3_M15_wavecalib.yaml (id 3_HR-R)
    4 valid, 0 invalid, 0 not checked

The exit status of ``numina verify`` is 1 if some file is invalid, so it
can be used in scripts. The structured calibrations (JSON files) can be
checked too. ``numina verify`` is available in the versions of numina
after 0.37.1.

**************************
Checking images in Python
**************************

In Python, an image can be checked with the function that MEGARA
registers in numina, the same used by ``numina verify``:

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
