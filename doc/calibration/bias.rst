**********
Bias Image
**********

:Mode: Bias
:Usage: Offline, Online
:Key: MegaraBIAS_IMAGE
:Product: :class:`~megaradrp.types.MasterBias`
:Recipe: :class:`~megaradrp.recipes.calibration.bias.BiasRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.bias.BiasRecipe.RecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.bias.BiasRecipe.RecipeResult`

Before the Analog-to-Digital conversion is performed a pedestal (electronic)
level is added to all images obtained with the MEGARA CCD. This is a standard
procedure in CCD imaging and spectroscopy applications for Astronomy and is
intended to minimize the ADC errors produced when very low analog values are
converted to DUs.


Requirements of the mode
++++++++++++++++++++++++
The sequence for this observing mode should include the actions to calibrate
the pedestal level of the detectors and associated control electronics by
taking images with null integration time. This mode requires having the shutter
closed and to readout the detector in a series of exposures with null
integration time, being this series the bias image set.

Procedure
+++++++++
The frames in the observed block are stacked together using the median of them
as the final result. The variance of the result frame is computed using two
different methods. The first method computes the variance across the pixels in
the different frames stacked. The second method computes the variance in each
channel in the result frame.

Products
++++++++

Bias image sets are to be obtained both as part of the activities related to
the verification of the instrument status and for processing data for
scientific exploitation.


Recipe, inputs and results
++++++++++++++++++++++++++


.. autoclass:: megaradrp.recipes.calibration.bias.BiasRecipe
      :members:


