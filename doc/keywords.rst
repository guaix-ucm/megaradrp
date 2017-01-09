
FITS Keywords
=============
The FITS keywords used by MEGARA are described in full detail
elsewhere (document TEC/MEG/146).

In the following sections, we describe the keywords that are used by the pipeline.

TBD

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
 R   EXPOSED                             Photometric time?
 R   DARKTIME                            TBD
 R   EXPTIME                             TBD
 R   ELAPSED                             Time between resets?
 S   OBSMODE   'SLITFLAT                 Identifier of the observing mode
==== ========  ========================  ============================================
