Scientific Modes
================

The observing modes described in this section are those intended for the
acquisition of scientific data by the observer. Here we describe all possible
scientific observations to be carried out either with one of the MEGARA IFUs or
the Fiber MOS.

LCB IFU scientific observation
------------------------------
:Usage: Online and Offline

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

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
+------------------------------+-------------------------------------------------------+

Procedure
+++++++++

Products
++++++++
The observer will obtain LCB IFU image sets as part of the routine scientific
operation of the instrument. The observatory staff could also make use of this
observing mode to verify the status of the instrument using any source
different from a standard star. In the case of observing a standard star the
calibration mode standard star with the LCB IFU could be used instead.

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
+------------------------------+-------------------------------------------------------+



LCB IFU fast mapping
--------------------
:Usage: Online and Offline

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

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
+------------------------------+-------------------------------------------------------+

Procedure
+++++++++

Products
++++++++
The observer will obtain LCB IFU image sets as part of the routine scientific
operation of the instrument. The observatory staff could also make use of this
observing mode to verify the status of the instrument using any source
different from a standard star.

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
+------------------------------+-------------------------------------------------------+



Fiber MOS scientific observation
--------------------------------
:Usage: Online and Offline

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

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
+------------------------------+-------------------------------------------------------+

Procedure
+++++++++

Products
++++++++
The observer will obtain Fiber MOS image sets as part of the routine scientific
operation of the instrument. The observatory staff could also make use of this
observing mode to verify the status of the instrument.

+------------------------------+-------------------------------------------------------+
| Name                         | Type                                                  |
+==============================+=======================================================+
+------------------------------+-------------------------------------------------------+