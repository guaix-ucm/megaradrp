Scientific Modes
================

The observing modes described in this section are those intended for the
acquisition of scientific data by the observer. Here we describe all possible
scientific observations to be carried out either with one of the MEGARA IFUs or
the Fiber MOS.

LCB IFU scientific observation
------------------------------

:Mode: LCB IFU scientific observation
:Usage: Online, Offline
:Key: MEGARA_LCB_IMAGE
:Recipe: :class:`~megaradrp.recipes.scientific.lcb.LCBImageRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.lcb.LCBImageRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.lcb.LCBImageRecipeResult`

This mode sequence includes the required actions to obtain scientifically-valid
data with the LCB IFU instrument mode of MEGARA.

Requirements
++++++++++++
This mode requires having the focal-plane cover configured, the instrument
shutter open, to configure the VPH mechanism to select the grating to be used,
to set the instrument mode to LCB IFU, to move the focusing mechanism to the
position pre-defined for the specific VPH of choice, and to expose a certain
time and to readout the detector in a series of exposures, being this series
the LCB IFU image set.

A pre-requisite for this observing mode is to have previously executed a "Fine
acquisition with the LCB IFU" auxiliary observing mode (TBC depending on the
system absolute positioning precision and the observer requirements).

As part of the MEGARA on-line quick-look software the images to be obtained by
this observing mode should be processed and the spectra extracted so to produce
a view of the field at a selectable wavelength within the wavelength range
covered by the VPH of choice. As a number of Fiber MOS positioners are devoted
to the measure of the sky background simultaneously with the LCB IFU
observations the on-line quick-look software should be able to subtract the
spectrum of the sky from the spectra in each IFU spaxel and from the maps
generated above. This software should have information on the specific
focal-plane cover configuration being used.

In the case of observing relatively extended targets (comparable in size or
larger than the LCB IFU field of view) but a mapping observing mode is not
required to be used a blank-sky image set should be obtained using the same
instrumental configuration as for the science target.

Procedure
+++++++++

Products
++++++++
The observer will obtain LCB IFU image sets as part of the routine scientific
operation of the instrument. The observatory staff could also make use of this
observing mode to verify the status of the instrument using any source
different from a standard star. In the case of observing a standard star the
calibration mode standard star with the LCB IFU could be used instead.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.lcb.LCBImageRecipe
      :members:

LCB IFU fast mapping
--------------------

:Mode: LCB IFU fast mapping
:Usage: Online, Offline
:Key: MEGARA_LCB_FAST_MAPPING
:Recipe: :class:`~megaradrp.recipes.scientific.lcbfastmapping.LCBFastMappingRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.lcbfastmapping.LCBFastMappingRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.lcbfastmapping.LCBFastMappingRecipeResult`


This mode sequence includes the required actions to obtain scientifically valid
data with the LCB IFU instrument mode of MEGARA covering a relatively large
area of the sky (of an extension equivalent to several LCB FOV). As this mode
is intended for the observation of large areas of the sky the exposure times of
individual exposures are expected to be short. That implies that the use of the
A&G system is not mandatory. By skipping the use of the A&G system throughout
the mapping process the observing efficiency dramatically increases. A (slow)
mapping observation, where each position in the grid should make use of the A&G
system, would have to be defined by the means of individual LCB IFU scientific
observing modes. The current observing mode only contemplates the possibility
of using the A&G system for the starting position of the mapping grid.

Requirements
++++++++++++
This mode requires having the focal-plane cover configured, the instrument
shutter open, to configure the VPH mechanism to select the grating to be used,
to set the instrument mode to LCB IFU, to move the focusing mechanism to the
position pre-defined for the specific VPH of choice, to expose a certain time
and to readout the detector in a series of exposures, and to repeat the data
acquisition in a regular grid of positions on the sky, being this series the
LCB IFU fast-mapping image set.

The execution of a "Fine acquisition with the LCB IFU" auxiliary observing mode
is optional in this case and it would be only possible for the starting
position in the grid.

As part of the MEGARA on-line quick-look software the images to be obtained by
this observing mode should be processed and the spectra extracted so to produce
a view of the entire mapped field at a selectable wavelength within the
wavelength range covered by the VPH of choice. The user can provide a blank-sky
position away from the mapping area. As a number of Fiber MOS positioners are
devoted to the measure of the sky background simultaneously with the LCB IFU
observations, alternatively, the on-line quick-look software could also
subtract the spectrum of the sky from the spectra in each IFU spaxel and from
the maps generated above. In the latter case (i.e. when no blank-sky position
is defined) showing sky-subtracted maps is optional as the sky positioners
could be strongly contaminated by emission from the extended target. Note that,
in any case, the on-line quick-look software should have information on the
configuration of the focal-plane cover as well.


Procedure
+++++++++

Products
++++++++
The observer will obtain LCB IFU image sets as part of the routine scientific
operation of the instrument. The observatory staff could also make use of this
observing mode to verify the status of the instrument using any source
different from a standard star.



Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.lcbfastmapping.LCBFastMappingRecipe
      :members:

Fiber MOS scientific observation
--------------------------------

:Mode: Fiber MOS scientific observation
:Usage: Online, Offline
:Key: MEGARA_MOS_IMAGE
:Recipe: :class:`~megaradrp.recipes.scientific.mos.MOSImageRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.mos.MOSImageRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.mos.MOSImageRecipeResult`


This mode sequence includes the required actions to observe a list of targets
with known celestial coordinates with MEGARA using the Fiber MOS instrument
mode. The information on the assignment of targets and positioners is included
in the form of a set of input catalogues generated off-line by the MEGARA
Observing Preparation Software Suite (MOPSS). The reference position for each
positioner is the center of the central fiber of the 7-fiber minibundle. This
observing mode could be run with one of the sides of the focal-plane cover
closed in order to reduce the cross-talk between positioners that would be
placed adjacent in the pseudo-slit. Thus, the input catalogues should be
specific of the focal-plane cover configuration to be used and both the on-line
quick-look software and the off-line pipeline should include that information
as a parameter or set of parameters. Note that the Fiber MOS catalogues and
configuration files could be designed for their use under any focal-plane cover
configuration. Even in that case the data processing software should know under
which configuration a given image set was obtained.

Requirements
++++++++++++
This observing mode requires having the focal-plane cover configured, the
instrument shutter open, to configure the VPH mechanism to select the grating
to be used, to set the instrument mode to Fiber MOS, to move the focusing
mechanism to the position pre-defined for the specific VPH of choice, to move
all robotic positioners with a target associated in the input catalogues to the
position of the corresponding target (these include science targets, reference
stars for fine acquisition and positioners devoted to blank-sky measurements)
and to expose a certain time and to readout the detector in a series of
exposures, being this series the Fiber MOS image set.

A pre-requisite for running this observing mode is to have previously executed
a "Fine acquisition with the Fiber MOS" auxiliary observing mode on the same
field.

As part of the MEGARA on-line quick-look software, the image (or images)
obtained should be processed and the spectra extracted. The observer might
define a number of positioners to be placed on blank-sky regions of the field
in order to improve sky subtraction. Alternatively, the user can also define a
blank sky position. This is particularly important when observing individual
stars in a nearby (Local Group) galaxy, for example, where the emission from
the host galaxy is expected to contaminate even the outermost positioners.
Should that be the case, the on-line quick-look software should be able to
derive a sky spectrum from the blank-sky observation (if present) or the
spectra of these positioners (if defined and no blank-sky observation is
available) and subtract it from the spectra of the targets. The processed
spectra should then be visualized using the on-line quick-look software. If
neither a blank-sky observation nor blank-sky positioners are available no sky
subtraction will be performed.


Procedure
+++++++++

Products
++++++++
The observer will obtain Fiber MOS image sets as part of the routine scientific
operation of the instrument. The observatory staff could also make use of this
observing mode to verify the status of the instrument.


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.mos.MOSImageRecipe
      :members:

Standard star with the LCB IFU
------------------------------

:Mode: Standard start with the LCB IFU
:Usage: Offline
:Key: MEGARA_LCB_STD_STAR
:Recipe: :class:`~megaradrp.recipes.scientific.lcbstdstar.LCBStandardRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.lcbstdstar.LCBStandardRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.lcbstdstar.LCBStandardRecipeResult`

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

* Absolute-flux calibration: The weather conditions during the night should be photometric and a number of spectrophotometric standard stars at different airmasses should be observed. This allows to fully correct from DUs per CCD pixel to energy surface density (typically in erg s-1 cm-2 Å-1) incident at the top of the atmosphere. If only one single standard star is observed (at the airmass of the science object) this correction allows deriving the energy surface density hitting the telescope primary mirror exclusively, unless an atmospheric extinction curve for the observatory and that particular night is assumed. In order to properly flux-calibrate scientific observations at all airmasses several stars should be observed during the night.

* Relative-flux calibration: If the weather conditions are not photometric this correction only allows normalizing the DUs per CCD pixel along the spectral direction so the conversion to incident energy at the top of the atmosphere is the same at all wavelengths. In order for this calibration to be valid the assumption that the effect of the atmosphere (including atmospheric cirrus and possibly thick clouds) on the wavelength dependence of this correction is that given by the atmospheric extinction curve adopted.

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


Procedure
+++++++++

Products
++++++++
Standard star image sets are to be obtained only as part of the routine
calibration activities performed by the observer and that are needed for
processing data for scientific exploitation.


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.lcbstdstar.LCBStandardRecipe
      :members:

Standard star with the Fiber MOS
--------------------------------

:Mode: Standard start with the FIBER MOS
:Usage: Offline
:Key: MEGARA_MOS_STD_STAR
:Recipe: :class:`~megaradrp.recipes.scientific.mosstdstar.MOSStandardRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.mosstdstar.MOSStandardRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.mosstdstar.MOSStandardRecipeResult`


This observing mode includes the required actions to obtain those calibration
images needed to correct for the variation in the response of the system along
the spectral direction. The difference between this mode and the two precedent
observing modes is that in this case the spectrophotometric standard star is
observed through one of the robotic positioners of the Fiber-MOS subsystem.

As in Standard star with the IFUs observing modes, the calibration is performed
by observing one or several spectrophotometric stars with the same instrument
configuration that for the scientific observations. Depending on the number of
standard stars observed and the weather conditions two different types of
calibration could be achieved, absolute or relative. In the case of the former
calibration an aperture correction should be applied to take into account the
possible flux losses from the standard stars when observed through one of the
~1.6-arcsec-wide robotic positioners.


Requirements
++++++++++++
This mode requires having the focal-plane cover configured, the instrument
shutter open, to configure the VPH mechanism to select the grating to be used,
to set the instrument mode to Fiber MOS, to move the focusing mechanism to the
position pre-defined for the specific VPH of choice, to move one of the robotic
positioners to the position of the spectrophotometric standard star (other
positioners could be also moved if needed) and to expose a certain time and to
readout the detector in a series of exposures, being this series the image set
containing the spectral energy distribution of the spectrophotometric standard
star.

This observing mode could still be carried out with one of the sides of the
focal-plane cover closed. However, as the (commonly rather bright)
spectrophotometric standard star is the only object of interest in the field,
the other positioners would not be observing scientific targets, so the level
of cross-talk between these and the positioner devoted to the standard star
should be negligible. Thus, the use of the focal-plane cover, although
considered, is not recommended for this specific observing mode.

In order to place the robotic positioner(s) on the corresponding target(s) a
set of input catalogues previously generated by the observer using MOPSS
(MEGARA Observing Preparation Software Suite) are needed.


Procedure
+++++++++

Products
++++++++
Standard star image sets are to be obtained only as part of the routine
calibration activities performed by the observer that are needed for processing
data for scientific exploitation.


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.mosstdstar.MOSStandardRecipe
      :members:

MEGARA_FLUX_CALIBRATION
-----------------------

:Mode:
:Usage: Offline
:Key: MEGARA_FLUX_CALIBRATION
:Recipe: :class:`~megaradrp.recipes.calibration.fluxcal.PseudoFluxCalibrationRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.fluxcal.PseudoFluxCalibrationRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.fluxcal.PseudoFluxCalibrationRecipeResult`


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.fluxcal.PseudoFluxCalibrationRecipe
      :members:

MEGARA_EXTINCTION_STAR
----------------------

:Mode:
:Usage: Offline
:Key: MEGARA_MOS_STD_STAR
:Recipe: :class:`~megaradrp.recipes.scientific.extinctionstar.ExtinctionStarRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.extinctionstar.ExtinctionStarRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.extinctionstar.ExtinctionStarRecipeResult`

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.extinctionstar.ExtinctionStarRecipe
      :members:

MEGARA_SENSITIVITY_STAR
-----------------------

:Mode:
:Usage: Offline
:Key: MEGARA_FLUX_CALIBRATION
:Recipe: :class:`~megaradrp.recipes.scientific.sensivitystar.SensivityStarRecipe`
:Recipe input: :class:`~megaradrp.recipes.scientific.sensivitystar.SensivityStarRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.scientific.sensivitystar.SensivityStarRecipeResult`

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.scientific.sensivitystar.SensivityStarRecipe
      :members:
