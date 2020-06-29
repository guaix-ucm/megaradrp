Calibration Modes
=================

The calibration modes are those operating modes that are intended (1) to either
analyze the state of the instrument or (2) to be used for processing scientific
observations from raw to science-grade data.

With respect to the determination of the status of the instrument the following
calibration images should be acquired:

* Bias

* Dark

* Fiber-flat

* Arc

Regarding the processing of scientific data this basically implies obtaining
calibration images in number and quality required to remove the instrumental
signatures so to obtain a science-grade image. These images, which are taken as
part of routine scientific operations, include:

* Bias

* Dark

* Slit-flat

* Fiber-flat

* Twilight fiber-flat

* Arc

* Standard star

Except for the slit-flat, that might be taken only occasionally, the rest of
this latter set of observing modes will be taken routinely as part of either
daytime or nighttime operations. We will refer to these as
"Daily CalibrationModes". Besides these modes we have identified a series of
calibration modes (named "System Calibration Modes") that are also necessary
for processing MEGARA observations but that are only produced occasionally as
part of long-term calibrations of the instrument to be carried out by the
observatory staff.

Thus, the "System Calibration Modes" will be:

* Bad-pixels mask

* Linearity Test

* Slit-flat. Whether the slit-flat should be considered as a "Daily Calibration" or "System Calibration" mode is TBD and will depend on the stability of the pixel-to-pixel efficiency of the MEGARA CCD.

In this latter case, the difference between a "System Calibration Mode" and the
corresponding "Auxiliary Mode" described in Section 3 depends on the frequency
the observing mode has to be executed. Auxiliary modes are typically run once
every observing run (e.g. the fine-acquisition ones) or, in the best case,
after a long period of inactivity. System Calibration modes, on the other hand,
are expected to be run only after major changes in the telescope or the
instrument or if a degradation of any of the subsystems of the instrument is
suspected.

The need for obtaining all these sets of images drives the requirements and
characteristics of the Calibration modes described below as defined by the
MEGARA team.


.. toctree::
   :maxdepth: 1

   bias
   dark
   slitflat
   trace
   model
   arc
   fiberflat
   twilight
   bpm
   lcbstd
   mosstd
