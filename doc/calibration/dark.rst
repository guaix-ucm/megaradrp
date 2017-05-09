Dark Image
----------

:Mode: Dark
:Usage: Offline, Online
:Key: MegaraDarkImage
:Product: :class:`~megaradrp.types.MasterDark`
:Recipe: :class:`~megaradrp.recipes.calibration.dark.DarkRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.dark.DarkRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.dark.DarkRecipeResult`

The potential wells in CCD detectors spontaneously generate electron-ion pairs
at a rate that is a function of temperature. For very long exposures this
translates into a current that is associated with no light source and that is
commonly referred to as dark current.

Requirements
++++++++++++
While in imaging or low-resolution spectroscopy this is nowadays a negligible
effect thanks to the extremely low dark current levels of state-of-the-art CCDs
(typically < 1 e-/hour) when working at intermediate-to-high spectral
resolutions where the emission per pixel coming from the sky background and the
astronomical source can be very low this is worth considering.


The sequence for this observing mode should include the actions to measure the
variation of the intrinsic signal of the system by taking images under zero
illumination condition and long integration time. This mode requires that the
focal-plane cover is configured (it should be fully closed), the shutter is
closed and to expose a certain time and readout the detector a series of
exposures, being this series the dark image set.

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Dark.
This mode includes the required actions to obtain a master dark frame. The
master dark generated is used in other stages of the data processing.

Products
++++++++
Dark image sets are to be obtained both as part of the activities related to
the verification of the instrument status and for processing data for
scientific exploitation.

A bidimensional dark image, QA flag, a text log file of the processing and a
structured text file containing information about the processing.


Recipe, inputs and results
++++++++++++++++++++++++++

.. autoclass:: megaradrp.recipes.calibration.dark.DarkRecipe
      :members:
