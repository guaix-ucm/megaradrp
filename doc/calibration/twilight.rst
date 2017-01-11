Twilight fiber-flat
-------------------

:Mode: Twilight fiber-flat
:Usage: Offline
:Recipe class: :class:`~megaradrp.recipes.calibration.twilight.TwilightFiberFlatRecipe`

Depending on the final performance of the ICM (provided by the GTC) at F-C the
twilight fiber-flat mode (proposed in this section) might be offered as
optional to the observer or a must should a proper data reduction be required.
In any case this must be always available as an observing mode.

The twilight fiber-flat observing mode should include the actions required to
calibrate the low-frequency sensitivity variation in the spatial direction of
the detector. In principle, the lamp fiber-flat should suffice to correct the
change in sensitivity along both the spatial (fiber-to-fiber relative
transmission) and the spectral direction of the system. The latter only
combined with flux standard-star observations since the spectral shape of the
ICM lamps is not known with enough accuracy.

The twilight fiber-flat is based on the observation of the blank twilight sky.
This can safely assume to homogeneously illuminate the entire MEGARA field of
view (3.5 arcmin x 3.5 arcmin).

Requeriments
++++++++++++
The focal plane should be uniformly illuminated with twilight-sky light. As the
illumination conditions change during twilight, each image set has a different
exposure time. The purpose is to obtain a similar (linear) level of DUs at the
detector (counts) under different illumination conditions.

This mode requires having the focal-plane cover configured (at least one of the
sides should be open), the instrument shutter open, the telescope tracking, to
move the pseudo-slit to that of the instrument mode of choice, to configure the
VPH wheel mechanism in order to select the grating to be used, to move the
focusing mechanism to the position pre-defined for the specific VPH of choice
and to take a series of exposures with different exposure times and to readout
the detector for this series of exposures, being these series the twilight
image set, each with a different exposure time, but with similar level of
counts.

+--------------------------+---------------+------------+-------------------------------+
| Name                     | Type          | Default    | Meaning                       |
+==========================+===============+============+===============================+
| ``'obresult'``           | Product       | NA         |      Observation Result       |
+--------------------------+---------------+------------+-------------------------------+
| ``'master_bias'``        | Product       | NA         |      Master Bias frame        |
+--------------------------+---------------+------------+-------------------------------+

Procedure
+++++++++

The "User" processes an observing block obtained in the observing mode Twilight
Fiber Flat. This mode includes the required actions to obtain a master
illumination flat field. The master illumination flat field generated is used
in other stages of the data processing.

Products
++++++++
Twilight-sky fiber-flat image sets are expected to be obtained as part of the
routine calibration activities performed by the observer since are needed for
processing any scientific-valid data. Therefore, this observing mode should be
considered as part of the "Daily Calibration Modes".

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
| ``'reduced_frame'``          | :class:`~megaradrp.types.ProcessedFrame`              |
+------------------------------+-------------------------------------------------------+
| ``'reduced_rss'``            | :class:`~megaradrp.types.ProcessedRSS`                |
+------------------------------+-------------------------------------------------------+
| ``'master_twilight_flat'``   | :class:`~megaradrp.products.MasterTwilightFlat`       |
+------------------------------+-------------------------------------------------------+

A RSS master illumination flat field, QA flag, a text log file of the
processing and a structured text file containing information about the
processing.

.. autoclass:: megaradrp.recipes.calibration.twilight.TwilightFiberFlatRecipe
      :members: