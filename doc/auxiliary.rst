Auxiliary Modes
===============
The auxiliary modes refer to those tasks that are carried out in preparation
for other observing modes and that are needed to ensure the quality of the
MEGARA observations. In general, these observing modes are not useful per se,
neither for straight scientific exploitation nor even data-processing purposes.

These Auxiliary modes are intended to (1) allow the observatory staff to
prepare the instrument for its optimal exploitation or (2) to improve the
performance of the telescope procedures in terms of the target acquisition and
focus.

With regard to the former, the following modes are defined:

* Telescope focus

* Spectrograph focus

These observing modes might have to be run before any observing run (in case of
visitor mode observations) and certainly after a long period of inactivity of
the instrument. Below we provide a detailed description of each of these
observing modes. Regarding the auxiliary observing modes defined to improve the
default telescope procedures these are:

* (The two previous observing modes)

* Fine acquisition with the LCB IFU

* Fine acquisition with the Fiber MOS

These latter three observing modes should be run before any observing run (in
case of visitor mode observations) and after a long period of inactivity of the
instrument; at least in the instrument mode to be used (LCB or MOS).


Telescope focus
---------------

:Mode: Telescope Focus
:Usage: Online
:Key: MegaraFocusTelescope
:Recipe: :class:`~megaradrp.recipes.auxiliary.focustel.FocusTelescopeRecipe`
:Recipe input: :class:`~megaradrp.recipes.auxiliary.focustel.FocusTelescopeRecipe.FocusTelescopeRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.auxiliary.focustel.FocusTelescopeRecipe.FocusTelescopeRecipeResult`

This observing mode includes the required actions to focus GTC using MEGARA. A
bright point source should be identified for this purpose. This mode is an
alternative to obtain the best focus of the telescope using the A&G system.
Since this mode requires having precise information on the spatial distribution
of the flux coming from a point source in a region a few arcsec in diameter the
use of the LCB  IFU bundle.

Requirements
++++++++++++
This mode requires the observation of a bright point source on the sky
continuously while the observing mode is being run. Although photometric
conditions are not needed the transparency should allow to properly measure the
source FWHM on the sky in any exposure of the series. The FOV of the LCB IFU is
large enough for these observations to be carried out with one of the sides of
the focal-plane cover closed. However, as the default configuration of the
instrument with the focal-plane cover in the open position the alignment will
be more likely carried out with the cover in that position.

This mode requires having the focal-plane cover configured (at least one of the
sides should be open; the most likely configuration will be with the
focal-plane cover fully open), the instrument
shutter open, to configure the VPH mechanism to select the grating to be used,
to set the instrument mode to LCB IFU, to move the focusing
mechanism to the position pre-defined for the specific VPH of choice, to expose
a certain time and to readout the detector in a series of exposures, being this
series the telescope focus image set. A pause in between every exposure in the
series should be introduced in order to give time for the M2 to adjust each new
focus position in the series and for the M2 control system to inform about its
new position, which should then re-start the observing sequence.

Products
++++++++
The observatory staff should obtain the telescope focus image sets as part of
standard preparatory observations (according to GRANTECAN this is usually done
every night). Should the focus offset between the ASG and SFS arms (or an
imaging instrument in other focus) and that of the MEGARA FC be stable
overtime, the use of this mode would be limited to the early stages of
characterization of the optimal telescope configuration for MEGARA, or after
major changes in the instrument, or if any problem arises. The observatory
staff should obtain the telescope focus image sets as part of standard
preparatory observations obtained at twilight (although not necessarily every
night) or when problems with the telescope focus model are suspected. Also,
telescope focus should be revised periodically (e.g. monthly) to correct for
potential temperature effects. In general, Auxiliary modes will be typically
run once every observing run (e.g. the fine-acquisition ones) or, in the best
(most relaxed) case, after a long period of inactivity.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.auxiliary.focustel.FocusTelescopeRecipe
      :members:

Spectrograph focus
------------------

:Mode: Spectrograph Focus
:Usage: Online
:Key: MegaraFocusSpectrograph
:Recipe: :class:`~megaradrp.recipes.auxiliary.focusspec.FocusSpectrographRecipe`
:Recipe input: :class:`~megaradrp.recipes.auxiliary.focusspec.FocusSpectrographRecipe.FocusSpectrographRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.auxiliary.focusspec.FocusSpectrographRecipe.FocusSpectrographRecipeResult`


This mode sequence includes the required actions to focus the MEGARA
spectrograph. The arc lamps from the ICM will be used for this purpose. MEGARA
includes a focusing mechanism at the position of the pseudo-slits. The best
image quality for the spectrograph is achieved by using a different focus
position for each disperser element (VPH). The focus position is independent of
the instrument mode in use (TBC; see next paragraph).

As said above, the two pseudo-slits (LCB and MOS) have to be at the same
exact place at the spectrograph entrance to yield the same focus position. If
this is not the case MEGARA can correct from this effect by focusing with a
fixed offset for all the VPHs. This value has to be added to the nominal
focusing position of a VPH with the pseudo-slit used as reference and shall be
taken into account in the MEGARA Control System look-up table. Whether these
should be one look-up table (i.e. a different focus compromise) for each
focal-plane cover configuration is TBD.


Requirements
++++++++++++
This mode requires having the focal-plane cover configured (at least one of the
sides should be open), the instrument shutter
open, to configure the VPH mechanism to select the grating to be used, to set
the instrument mode to use, to move the focusing mechanism to the pre-defined
focus position for specific VPH of choice, to expose a certain time and to
readout the detector in a series of exposures, being this series the
spectrograph focus image-sets. The focusing mechanism should be able to change
its position in between every two exposures in the series. A pause in between
every exposure should be introduced in order to give time for the focusing
mechanism to adjust to the new focus position in the series and for the MEGARA
control system to inform about its new position, which should then re-start the
observing sequence (see below). Whether the focus positions (or range) for each
VPH are already pre-defined is TBD.

As part of the on-line quick-look software all images in the series should be
pre-processed and several spectra along the pseudo-slit should be extracted and
analyzed. This analysis should include the computation of the FWHM of a few
unresolved spectra lines at different wavelengths. This software should then
decide based on the FWHM values computed at different wavelengths and positions
along the pseudo-slit the best focus compromise. The best focus obtained for
the VPH of choice should then be stored and used to determine the best foci for
all spectral configurations (and instrument modes; TBC).

Procedure
+++++++++
Spectrograph focus image sets through; at least, one of the MEGARA VPHs should
be obtained at the beginning of every observing night by either the observer or
the staff of the observatory (TBD). Once a VPH is checked, the rest of the
values could be corrected relative to this one. It is expected that minor focus
corrections should be done as the temperature changes. This could be modeled in
further phases and checked at laboratory and/or at the telescope. The
observatory staff should obtain an entire sequence of spectrograph focus image
sets through all VPHs (and instrument modes; TBC) after major changes in the
instrument, long periods of inactivity or when the relative-focus prescriptions
(i.e. the spectrograph focus model) are suspected to be inaccurate.

The focus difference (obtained by measuring a particular VPH) will provide the
offset focus (due to temperature) and this value will be the same for all VPHs.
The Control System will be prepared to update the look-up table with this
offset focus value due to temperature.

Products
++++++++

The best focus, the goodness of the fit of the best focus, a table with the
FWHM of the spectral line corresponding to each focus, position along the slit
and wavelength, the collapsed PSFs, QA flag, a text log file of the processing
and a structured text file containing information about the processing.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.auxiliary.focusspec.FocusSpectrographRecipe
      :members: FocusSpectrographRecipeInput, FocusSpectrographRecipeResult


Fine acquisition with the LCB IFU
---------------------------------

:Mode: LCB Acquisition
:Usage: Online
:Key: MegaraLcbAcquisition
:Recipe: :class:`~megaradrp.recipes.auxiliary.acquisitionlcb.AcquireLCBRecipe`
:Recipe input: :class:`~megaradrp.recipes.auxiliary.acquisitionlcb.AcquireLCBRecipe.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.auxiliary.acquisitionlcb.AcquireLCBRecipe.RecipeResult`


This mode sequence includes the required actions to acquire a target with known
celestial coordinates and place it at a reference position inside the LCB IFU
instrument mode. The reference position for each mode is defined as the center
of the fibers (or its associated microlens) that is closest to the bundle
footprint geometrical center. In the case of the LCB the reference position
will depend on the focal-plane cover configuration. This mode is a refinement
of acquisition performed by the telescope or A&G systems.

Requirements
++++++++++++
This mode requires having the focal-plane cover configured, the instrument
shutter open, to configure the VPH mechanism to select the grating to be used,
to set the instrument mode to LCB, to move the focusing mechanism to the
position pre-defined for the specific VPH of choice, and to expose a certain
time and to readout the detector in a series of exposures, being this series
the fine acquisition image set.

As part of the MEGARA on-line quick-look software the image (or images)
obtained as part of this observing mode should be processed and the spectra
extracted so to determine the position of the centroid of the target in the
corresponding field of view. A view of the field should be also produced in
order to evaluate whether or not the angle of the Folded-Cass rotator matches
that specified by the observer.

Products
++++++++
Fine acquisition image sets should be obtained at the beginning of the
observing night by either the observer or the staff of the observatory (TBD) or
every time a problem with the telescope absolute pointing is suspected. Such
image sets should be also obtained when an absolute positioning precision of
the order of a fraction of the spaxel size is required, better than 0.62 arcsec
in this case for the LCB.

The observatory staff should decide whether or not the corrections derived
must be applied to the acquisition of other targets during the same observing
night or exclusively to the target currently being observed.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.auxiliary.acquisitionlcb.AcquireLCBRecipe
      :members:

Fine acquisition with the Fiber MOS
-----------------------------------

:Mode: MOS Acquisition
:Usage: Online
:Key: MegaraLcbAcquisition
:Recipe: :class:`~megaradrp.recipes.auxiliary.acquisitionmos.AcquireMOSRecipe`
:Recipe input: :class:`~megaradrp.recipes.auxiliary.acquisitionmos.AcquireMOSRecipe.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.auxiliary.acquisitionmos.AcquireMOSRecipe.RecipeResult`

The sequence for this observing mode includes the required actions to acquire a
list of targets with known celestial coordinates and place each target at the
center of a different robotic positioner. The information on the assignment of
targets and positioners is included in the form of a set of input catalogues
generated off-line by the MEGARA Observing Preparation Software Suite (MOPSS).
The reference position for each positioner is the center of the central fiber
of the 7-fiber minibundle. This mode is a refinement of the acquisition
performed by the telescope or A&G systems.

Requirements
++++++++++++
This mode requires having the focal-plane cover configured, the instrument
shutter open, to configure the VPH mechanism to select the grating to be used,
to set the instrument mode to Fiber MOS, to move the focusing mechanism to the
position pre-defined for the specific VPH of choice, to move all robotic
positioners with a target associated in the input catalogues to the position of
the corresponding target and to expose a certain time and to readout the
detector in a series of exposures, being this series the Fiber-MOS fine
acquisition image set.

As part of the MEGARA on-line quick-look software, the image (or images)
obtained should be processed and the spectra extracted so to determine the
position of the centroid of a number of reference targets included in the
corresponding field of view and identified as such in the set of input
catalogues used for this observing mode. A minimum of three reference sources
should be included in each Fiber MOS configuration block in order for this
observing mode to generate a solution. The quick-look software should compare
the expected and the actual positions of these reference sources in order to
determine the best-fitting set of offsets (both in X and Y) and rotation angle
to apply to the telescope and Folded-Cass rotator, respectively, to then
continue with one of the scientific observing modes described in next Section.

Products
++++++++
Fine acquisition image sets should be obtained by the observer at the beginning
of the observation of each field with the Fiber MOS. The observatory staff
should decide whether or not the corrections derived (telescope offset and
Folded-Cass rotator angle) must be applied to the acquisition of other fields
with the Fiber MOS during the same observing night or exclusively to the target
currently being observed.

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.auxiliary.acquisitionmos.AcquireMOSRecipe
      :members:
