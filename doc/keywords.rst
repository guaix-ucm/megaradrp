
FITS Keywords
=============
The FITS keywords used by MEGARA are described in full detail
elsewhere (document TEC/MEG/146).

In the following sections, we describe the keywords that are used by the pipeline.


Primary header
--------------
==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 L   SIMPLE     T                        Standard FITS format
 I   BITPIX     16                       One of -64,-32,8,16,32
 I   NAXIS      2                        # of axes in frame
 I   NAXIS1    2048                      # of pixels per row
 I   NAXIS2    2048                      # of rows
 S   ORIGIN    'GTC'                     FITS originator
 S   OBSERVAT  'ORM'                     Observatory
 S   TELESCOP  'GTC'                     The telescope
 S   INSTRUME  'MEGARA'                  The instrument
 S   OBJECT    'NGC 4594'                Target designation
 S   OBSERVER  'OBSERVER'                Who adquired the data
 S   DATE-OBS  '2012-09-20T12:00:11.50'  Date of the start of the observation
 S   DATE      '2012-09-20T12:14:12.78'  Date the file was written
==== ========  ========================  ============================================

Required by the pipeline
------------------------

==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 R   AIRMASS   1.1908                    Mean airmass of the observation
 R   MJD-OBS   72343.34324               Modified JD of the start of the observation
 S   IMAGETYP  'FLAT'                    Type of the image
 S   VPH       'LR-R'                    Type of VPH
 S   OBSTYPE   'SLITFLAT'                Type of observation
 R   EXPOSED                             Exposure time in seconds
 R   EXPTIME                             Exposure time in seconds (synonim)
 R   DARKTIME                            TBD
 S   OBSMODE   'SLITFLAT                 Identifier of the observing mode
==== ========  ========================  ============================================

FIBERS extension
----------------

The state of the focal plane of MEGARA is stored in a dedicated extension names `FIBERS`.
This extension contains only headers, the data part will be empty.

==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 I   NFIBERS   643                       Number of fibers
 I   NSPAXEL   644                       Number of spaxels
 I   NBUNDLES  92                        Number of fiber bundles
 S   INSMODE   LCB                       Name of active pseudo slit
 S   CONFID    'b7d35e7df0274fde..'      Unique identificator of the configuration
 I   BUNnnn_P  0                         Priority of the target in this bundle
 S   BUNnnn_I  'unknown '                Name of the target
 S   BUNnnn_T  'UNASSIGNED'              Type of target ('STAR', 'SKY', 'TARGET', 'UNASSIGNED'
 I   FIBmmm_B   nnn                      ID of the bundle
 F   FIBmmm_D   +3.34565                 Declination of the spaxel
 F   FIBmmm_R   12.342223                Right Ascension of the spaxel
 F   FIBmmm_O   0.0                      Position Angle of the Fiber
 L   FIBmmm_A   T                        Is fiber active?
 F   FIBmmm_X   -0.4646226291303512      X position of the fiber in the focal plane
 F   FIBmmm_Y   63.63025                 Y position of the fiber in the focal plane
==== ========  ========================  ============================================
