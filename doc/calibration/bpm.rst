Bad-pixels mask
---------------

:Mode: Bad-pixels mask
:Usage: Offline, Online
:Key: MegaraBadPixelMask
:Product: :class:`~megaradrp.types.MasterBPM`
:Recipe: :class:`~megaradrp.recipes.calibration.bpm.BadPixelsMaskRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.bpm.BadPixelsMaskRecipe.BadPixelsMaskRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.bpm.BadPixelsMaskRecipe.BadPixelsMaskRecipeResult`

Although science-grade CCD detectors show very few bad pixels / bad columns
there will be a number of pixels (among the ~17 Million pixels in the MEGARA
CCD) whose response could not be corrected by means of using calibration images
such as dark frames or flat-field images. These pixels, commonly called either
dead or hot pixels, should be identified and masked so their expected signal
could be derived using dithered images or, alternatively, locally interpolated.
While a bad-pixels mask will be generated as part of the AIV activities, an
increase in the number of such bad pixels with time is expected. Therefore, we
here define an observing mode that the observatory staff could use to generate
an updated version of the bad-pixels masks should the number of bad pixels
increase significantly.

In the case of fiber-fed spectrographs the fiber flats (either lamp or twilight
flats) are not optimal for generating bad-pixels masks as these leave many
regions in the CCD not exposed to light. The whole CCD should be illuminated at
different intensity levels in order to clearly identify both dead and hot
pixels.

Requirements
++++++++++++

In the case of MEGARA we will offset the pseudo-slit from its optical focus
position to ensure that the gaps between fibers are also illuminated when a
continuum (halogen) lamp at the ICM is used. The NSC zemax model of the
spectrograph indicates that by offsetting 3mm the pseudo-slit we would already
obtain a homogenous illumination of the CCD. A series of images with different
count levels would be obtained.

This mode requires having the ICM halogen lamp on, the instrument shutter open,
to move the pseudo-slit to the open position, to configure the VPH wheel
mechanism in order to select the grating to be used, to move the focusing
mechanism to the position pre-defined for the specific VPH of choice but offset
by 3mm and to expose a certain time and to readout the detector a series of
exposures, being this series the slit-flat image set. Note that only one
Bad-pixels mask will be used for all spectral setups. The specific choice for
the VPH will depend on the actual color of the ICM halogen lamp and on the
actual response of the VPHs. In principle, we should choose the VPH at the peak
of the lamp spectral energy distribution but we should also consider the fact
that the VPH should have the flattest spectral response possible. We call this
specific VPH the "BPM VPH". LR-R and LR-I are currently the best candidates for
finally being the BPM VPH.

Procedure
+++++++++

The "User" processes an observing block obtained in the observing mode
Bad-pixels mask. This mode includes the required actions to obtain a bad-pixel
mask. The master bad pixel mask generated is used in other stages of the data
processing.

Products
++++++++
This Bad-pixels mask observing mode will be used only sporadically as it is
considered part of the "System Calibration Modes".

A bidimensional mask of bad pixels, a QA flag, a text log file of the
processing and a structured text file with information about the processing.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.bpm.BadPixelsMaskRecipe
      :members:

.. autoclass:: megaradrp.recipes.calibration.bpm::BadPixelsMaskRecipe.BadPixelsMaskRecipeInput
      :members:

.. autoclass:: megaradrp.recipes.calibration.bpm::BadPixelsMaskRecipe.BadPixelsMaskRecipeResult
      :members:
