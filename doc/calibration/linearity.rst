Linearity test
--------------

:Mode: Linearity test
:Usage: Offline, Online
:Key: MegaraLinearityTest
:Product:
:Recipe: :class:`~megaradrp.recipes.calibration.linearity.LinearityTestRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.linearity.LinearityTestRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.linearity.LinearityTestRecipeResult`

Although the linearity of the MEGARA CCD are well characterized at the LICA lab
already, it might be advisable to generate linearity test frames both as part
of the AIV activities and after changes in the MEGARA DAS.

The MEGARA e2V 231-84 CCD offers a full-well capacity of 350,000 ke-. Linearity
tests carried out in instruments already using this type of CCD indicate a
linearity better than ±0.4% at 100 kpix/sec in the range between 140 to 40,000
e- (Reiss et al. 2009 for MUSE@VLT). Given these good linearity results (up to
40,000 e-) and considering the fact that at the spectral resolutions of MEGARA
we will rarely reach those signals from astronomical targets linearity can be
considered negligible. Despite these facts, it is advisable to carry out this
kind of tests both at the lab and at the telescope on the MEGARA CCD itself.

While Linearity tests will be generated as part of the characterization
activities at the lab, the use of the ICM would also allow carrying them out as
part of AIV activities and routinely as part of the "System Calibration Modes".
Therefore, we define here an observing mode that the observatory staff could
use to generate updated Linearity tests should these be needed.

In the case of fiber-fed spectrographs the fiber flats (either lamp or twilight
flats) are not optimal for carrying out Linearity tests as these leave many
regions in the CCD not exposed to light. The whole CCD should be illuminated at
different intensity levels in order for properly carrying out these tests.


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
mechanism to the position pre-defined for the specific VPH of choice and to
expose a certain time and to readout the detector a series of exposures, being
this series the slit-flat image set. Note that the Linearity test will be done
using only one spectral setup as this is independent of the VPH of use. The
specific choice for the VPH will depend on the actual color of the ICM halogen
lamp and on the actual response of the VPHs. In principle, we should choose the
VPH at the peak of the lamp spectral energy distribution but we


Procedure
+++++++++

Products
++++++++

This Linearity-test observing mode will be used only sporadically as it is
considered part of the "System Calibration Modes".


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.linearity.LinearityTestRecipe
      :members:
