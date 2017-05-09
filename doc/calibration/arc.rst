Arc
---

:Mode: Arc
:Usage: Offline, Online
:Key: MegaraArcCalibration
:Product: :class:`~megaradrp.products.wavecalibration.WavelengthCalibration`
:Recipe: :class:`~megaradrp.recipes.calibration.arc.ArcCalibrationRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.arc.ArcCalibrationRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.arc.ArcCalibrationRecipeResult`

This mode sequence includes the required actions to translate the geometrical
position of each point in the detector into physical units of wavelength. The
calibration is performed by means of using reference calibration lamps
(arc lamps) that should be part of the Instrument Calibration Module (ICM) at
F-C. Note that the optical distortions in the spectrograph will lead to
different wavelength calibrations from each individual fiber, therefore the
entire focal plane should be illuminated by the corresponding arc lamp of
choice. Given the relatively high spectral resolution and broad wavelength
coverage of MEGARA we anticipate that more than one arc lamp will be needed at
the ICM. The lamps at ICM have to deliver enough bright spectral lines for
calibrating the whole range of MEGARA spectral resolutions and wavelength
ranges (for HR modes only two VPHs shall be provided by MEGARA but more
gratings could come funded by GTC users). MEGARA has provided a whole review of
the possible illumination systems in the document TEC/MEG/151, but the
responsibility of the development of the ICM module is on the GTC side.

Requirements
++++++++++++
The entire focal plane should be illuminated with light from the ICM arc lamp
with the required  input focal ratio. This mode requires having the ICM turned
on, one of the arc lamps at the ICM also turned on, the focal-plane cover
configured (at least one of the sides should be open), the instrument shutter
open, to move the pseudo-slit to that of the instrument mode of choice, to
configure the VPH wheel mechanism in order to select the grating to be used, to
move the focusing mechanism to the position pre-defined for the specific VPH of
choice and to expose a certain time and to readout the detector a series of
exposures, being this series the arc image set.

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Arc.
This mode includes the required actions to translate the geometrical position
of each point in the detector into physical units of wavelength. The wavelength
calibration generated is used in other stages of the data processing.

Products
++++++++

Arc image sets are to be obtained both as part of the activities related to the
verification of the instrument status and for processing data for scientific
exploitation and are part of the "Daily Calibration Modes".

A data structure containing information about wavelength calibrations
(the format is TBD), a QA flag, a text log file of the processing and a
structured text file containing information about the processing.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.arc.ArcCalibrationRecipe
      :members:
