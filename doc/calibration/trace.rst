=====
Trace
=====

:Mode: Trace
:Usage: Offline, Online
:Key: MegaraTraceMap
:Product: :class:`~megaradrp.products.tracemap.TraceMap`.
:Recipe: :class:`~megaradrp.recipes.calibration.trace.TraceMapRecipe`
:Recipe input: :class:`~megaradrp.recipes.calibration.trace.TraceMapRecipeInput`
:Recipe result: :class:`~megaradrp.recipes.calibration.trace.TraceMapRecipeResult`


Although for the majority of the observing modes described elsewhere in this
document the MEGARA off-line pipeline will perform its own fiber spectra
extraction from the 2D CCD FITS frame, there are cases where an archival master
"trace map" should be used instead. Note that a different "trace map" should be
available for each pseudo-slit and VPH combination.

Requirements
++++++++++++
This observing mode should include the actions needed to obtain a series of
Fiber-flats that should be combined to generate a master "trace map". This will
be done by means of illuminating the instrument focal plane with a continuum
(halogen) lamp that is part of the GTC Instrument Calibration Module (ICM). The
use of the twilight sky is not recommended in this case as the twilight sky can
present strong absorption lines that could lead to errors in the resulting
trace map at specific wavelengths.

This mode requires having the ICM turned on, one of the halogen lamps at the
ICM also turned on, to configure the focal-plane cover (at least one of the
sides should be open), to have the instrument shutter open, to move the
pseudo-slit to that of the instrument mode of choice, to configure the VPH
wheel mechanism in order to select the grating to be used, to move the focusing
mechanism to the position pre-defined for the specific VPH of choice and to
expose a certain time and to readout the detector a series of exposures, being
this series the trace map image set.

Procedure
+++++++++
The "User" processes an observing block obtained in the observing mode Trace.
This mode includes the required actions to obtain a mapping of the trace of the
fibers. The master trace map generated is used in other stages of the data
processing.

Products
++++++++

Trace map image sets are to be obtained both as part of the activities related
to the verification of the instrument status and for processing data for
scientific exploitation. Note, however, that the use of this observing mode for
scientific exploitation should be limited as it could affect to the general
performance of the on-line quick-look software.

This mode produces the tracing information required to extract the flux of the fibers.
The result is stored in an object named  ``master_traces``
of type :class:`~megaradrp.products.tracemap.TraceMap`.


Recipe
------

.. autoclass:: megaradrp.recipes.calibration.trace.TraceMapRecipe
   :members:
