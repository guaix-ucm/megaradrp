Standard star with the LCB IFU
------------------------------

:Mode: Standard star with the LCB IFU
:Usage: Offline, Online
:Key: MegaraLcbStdStar
:Recipe: :class:`~megaradrp.recipes.calibration.lcbstdstar.LCBStandardRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.lcbstdstar.LCBStandardRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.lcbstdstar.LCBStandardRecipeResult`

This observing mode includes the required actions to obtain those calibration
images needed to correct for the variation in the response of the system along
the spectral direction. This signature is manifested by a change in the
conversion factor between the energy surface density hitting the telescope
primary mirror and the DUs per CCD pixel with wavelength. Its effect is already
present in the original data but it could get modified during the reduction
process, e.g. after the fiber-flat correction is applied.

The flux calibration is performed by observing one or several
spectrophotometric stars with the same instrument configuration that for the
scientific observations. Depending on the number of standard stars observed and
on the weather conditions (mainly transparency) two different types of
calibration could be achieved:

* Absolute-flux calibration: The weather conditions during the night should be photometric and a number of
  spectrophotometric standard stars at different airmasses should be observed. This allows to fully correct from
  DUs per CCD pixel to energy surface density (typically in erg s-1 cm-2 Å-1) incident at the top of the atmosphere.
  If only one single standard star is observed (at the airmass of the science object) this correction allows deriving
  the energy surface density hitting the telescope primary mirror exclusively, unless an atmospheric extinction curve
  for the observatory and that particular night is assumed. In order to properly flux-calibrate scientific observations
  at all airmasses several stars should be observed during the night.

* Relative-flux calibration: If the weather conditions are not photometric this correction only allows normalizing the
  DUs per CCD pixel along the spectral direction so the conversion to incident energy at the top of the atmosphere is
  the same at all wavelengths. In order for this calibration to be valid the assumption that the effect of the
  atmosphere (including atmospheric cirrus and possibly thick clouds) on the wavelength dependence of this
  correction is that given by the atmospheric extinction curve adopted.

Since the observing sequence needed for both types of flux calibration is
identical only one observing mode (standard star) needs to be defined.

We will use this same observing mode also for the observation of either
telluric standards or radial-velocity standards. The former are needed to
correct for the presence of telluric absorptions mainly in the red part of the
spectrum and are achieved by means of observing A-type stars at the same
airmass and very close in time to the corresponding scientific observation.
The latter can used to determine a precise zero point velocity for the
instrument at a specific night and to verify its stability from night to night
and season to season.

Requirements
++++++++++++

This mode requires the entire flux of the spectrophotometric standard star to
be recovered (even if the star is a telluric or radial-velocity standard),
especially when an absolute-flux calibration is needed, so the LCB IFU bundle
must be used. The FOV of the LCB IFU is large enough for these observations to
be carried out with one of the sides of the focal-plane cover closed. When this
calibration is aimed for a set of Fiber-MOS scientific observations,
complementary observations of standard stars through the Fiber-MOS minibundles
might be also required. This allows verifying the quality and stability of the
calibration when two different pseudo-slits are used. Such observing mode is
described later.

This mode requires having the focal-plane cover configured (at least one of the
sides should be open), the instrument shutter open, to configure the VPH
mechanism to select the grating to be used, to set the instrument mode to LCB
IFU, to move the focusing mechanism to the position pre-defined for the
specific VPH of choice, and to expose a certain time and to readout the
detector in a series of exposures, being this series the image set containing
the spectral energy distribution of the spectrophotometric standard star.

In order to distribute the flux from the star across multiple spaxels in the
LCB IFU bundle (particularly important in the case of very bright
spectrophotometric standard stars) we might also need to apply a small drift
motion (typically of a few arcsec per second) to one of the telescope axes at
the start of the observation or, more likely, slightly defocus the telescope.

Products
++++++++
Standard star image sets are to be obtained only as part of the routine
calibration activities performed by the observer and that are needed for
processing data for scientific exploitation.


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.lcbstdstar.LCBStandardRecipe
      :members:
