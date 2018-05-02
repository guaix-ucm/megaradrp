Fiber-flat
----------

:Mode: Fiber-flat
:Usage: Offline, Online
:Key: MegaraFiberFlatImage
:Product: :class:`~megaradrp.types.MasterFiberFlat`
:Recipe: :class:`~megaradrp.recipes.calibration.flat.FiberFlatRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.flat.FiberFlatRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.flat.FiberFlatRecipeResult`

In fiber-fed spectrographs such as MEGARA each optical fiber behaves like a
different optical system, and therefore, its optical transmission is different
and individual, with different wavelength dependence. In the Preliminary Design
phase this mode was named "Lamp fiber flat".

Requirements
++++++++++++
This observing mode should include the actions to calibrate the low-frequency
variations in transmission in between fibers and as a function of wavelength in
MEGARA. A fiber-flat should be used to perform this correction and is the
result of illuminating the instrument focal plane with a flat source that can
be either a continuum (halogen) lamp that is part of the GTC Instrument
Calibration Module (ICM) or the twilight sky. The fiber-flat observing mode
discussed here assumes that the focal plane is illuminated with a halogen lamp
located at the ICM. The ICM beam has to have the same focal ratio arriving to
the first MEGARA optical element (the MEGARA telecentricity-correction lens in
this case) simulating as much as possible the real GTC mirrors beam at F-C.

These fiber-flat images are also used to trace the fiber spectra on the
detector for each specific spectral setup. Finally, they are also useful to
verify the status of the optical link between the F-C focal plane and the
platform where the spectrographs are located.

This mode requires having the ICM turned on, one of the halogen lamps at the
ICM also turned on, to configure the focal-plane cover (at least one of the
sides should be open), to have the instrument shutter open, to move the
pseudo-slit to that of the instrument mode of choice, to configure the VPH
wheel mechanism in order to select the grating to be used, to move the
focusing mechanism to the position pre-defined for the specific VPH of choice
and to expose a certain time and to readout the detector a series of exposures,
being this series the fiber-flat image set.


Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode
Fiber-flat. This mode includes the required actions to obtain a master
fiber-flat field. The master fiber-flat field generated is used in other stages
of the data processing.

Products
++++++++
Fiber-flat image sets are to be obtained both as part of the activities related
to the verification of the instrument status and for processing data for
scientific exploitation.

.. _ff-mode-label:

Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.flat.FiberFlatRecipe
   :members:
