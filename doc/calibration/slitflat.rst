Slit-flat
---------

:Mode: Slit-flat
:Usage: Offline, Online
:Key: MegaraSlitFlat
:Product: :class:`~megaradrp.types.MasterSlitFlat`
:Recipe: :class:`~megaradrp.recipes.calibration.slitflat.SlitFlatRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.slitflat.SlitFlatRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.slitflat.SlitFlatRecipeResult`

In the case of fiber-fed spectrographs the correction for the detector
pixel-to-pixel variation of the sensibility is usually carried out using data
from laboratory, where the change in efficiency of the detector at different
wavelengths is computed and then used to correct for this effect for each
specific instrument configuration (VPH setup in the case of MEGARA).

Requeriments
++++++++++++
In the case of MEGARA we will offset the pseudo-slit from its optical focus
position to ensure that the gaps between fibers are also illuminated when a
continuum (halogen) lamp at the ICM is used. The NSC zemax model of the
spectrograph indicates that by offsetting 3mm the pseudo-slit we would already
obtain a homogenous illumination of the CCD. A series of images with different
count levels would be obtained.

The quality of present-day CCDs leads to a rather small impact of these
pixel-to-pixel variations in sensitivity on either the flux calibration and the
cosmetics of the scientific images, especially considering that not one but a
number of pixels along the spatial direction are extracted for each fiber and
at each wavelength. Therefore, we anticipate that this correction might not be
needed or that, as a maximum, a first-order correction based on laboratory data
might suffice. However, before the results of the analysis of the
pixel-to-pixel variations in sensitivity planned using our CCD230 e2V test CCD
are obtained we will consider this observing mode as TBC.

This mode requires having the ICM halogen lamp on, the instrument shutter open,
to move the pseudo-slit to the open position, to configure the VPH wheel
mechanism in order to select the grating to be used, to move the focusing
mechanism to the position pre-defined for the specific VPH of choice but offset
by 3mm and to expose a certain time and to readout the detector a series of
exposures, being this series the slit-flat image set.

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode
Slit-flat. This mode includes the required actions to obtain a master slit-flat
field. The master slit-flat field generated is used in other stages of the data
processing.

Products
++++++++
Slit-flat image sets are to be obtained both as part of the activities related
to the verification of the instrument status (such as for evaluating the status
of the MEGARA spectrograph) and also for processing data for scientific
exploitation (correction for the pixel-to-pixel variation in sensitivity). The
frequency at which these detector flat images should be acquired is TBC.
Although defined in this document as a mode to be considered part of the
"Daily Calibration Modes" if it is finally used only sporadic it should be
considered as part of the "System Calibration Modes" instead.


A bidimensional master slit flat field, QA flag, a text log file of the
processing and a structured text file containing information about the
processing.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.slitflat.SlitFlatRecipe
      :members:
