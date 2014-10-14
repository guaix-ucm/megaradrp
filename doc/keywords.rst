
FITS Keywords
=============
The description of the keywords follows a convention found in other FITS 
keyword dictionaries, for example the list in http://fits.gsfc.nasa.gov/fits_dictionary.html. 
The keyword name is expressed, with the reference to the paper where it is 
included. Following the type of HDU where the keyword can appear. The value 
shows the kind of variable represented by the keyword. The comment is a 
example of the comment associated with the keyword and definition is a 
explanation in detail of the usage of the keyword.

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
 S   OBSTYPE   'FLATON'                  Type of observation
 L   READPROC   T                        The frame has been preprocessed after readout
 R   EXPOSED                             Photometric time?
 R   DARKTIME                            TBD
 R   EXPTIME                             TBD
 R   ELAPSED                             Time between resets?
 I   OBSID      567                      Identifier of the observing block
 S   OBSMODE   'BIAS'                    Identifier of the observing mode
==== ========  ========================  ============================================

Coordinate system
-----------------
The specifications of world coordinates are treated in a series of four papers.
By World coordinates we mean coordinates that serve to locate a measurement in 
some multidimensional parameter space. They include, for example, a measurable 
quantity such as frequency or wavelength associated with each point of the 
spectrum or a longitude and latitude in a conventional spherical coordinate 
system.

“Representation of world coordinates in FITS”, Greisen, E.W. & Clabretta, 
M.R. 2002 A&A 395, 1061 (hereafter Paper I) describes a very general method 
for specifying coordinates. A pixel-to-coordinate matrix PCj_i will replace 
CROTAj, units will be described with a new keyword CUNITj, and secondary 
sets of coordinate descriptions may be specified. A complete system of unit
specification is described and is expected to supplement the IAU standard 
system of units. Methods for describing the coordinates of matrices in binary 
tables are also described. 

“Representation of celestial coordinates in FITS”, Clabretta, M.R. & 
Greisen, E.W. 2002 A&A 395, 1077 (hereafter Paper II) applies the general 
rules of Paper I to the specific problem of specifying celestial coordinates 
in a two-dimensional projection of the sky. The coordinate system is 
specified with the new keyword RADESYS and a large number of projections are 
defined. Oblique projections are described and illustrated. Several examples 
of header interpretation and construction are given including one that 
specifies coordinates on a planetary body rather than the celestial sphere. 
The application to binary tables is described. 

“Representation of spectral coordinates in FITS”, Greisen, E.W. et al. 2004 (hereafter Paper III) is
still open to comments from the FITS community. It applies the general rules and practices developed
in the first two papers to spectral coordinates, namely frequency, wavelength, velocity, and the radio
and optical conventional velocities. These are defined and methods of computing one type of
coordinate from a spectral axis gridded in another are given. A projection representative of optical
spectrometers is also defined. Coordinate reference frames may be specified. 

“Representation of distortions in FITS world coordinate systems”, Clabretta, M.R. et al. 2004
(hereafter Paper IV) is in preparation. It will define Distortion Correction Functions (DCFs) which
may be used to correct for instrumental defects including celestial coordinate warps (plate defects),
variation of actual frequency with celestial coordinate, refraction, and the like.

The set of WCS keywords usable are those supported by wcslib library http://www.atnf.csiro.au/people/mcalabre/WCS/

Checksum convention
-------------------
The CHECKSUM and DATASUM keywords that are embedded in the FITS header are used to verify the integrity of the HDU.

See http://fits.gsfc.nasa.gov/registry/checksum.html

==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 S   CHECKSUM  'ADFASASDLIEXV'           HDU checksum 
 S   DATASUM   '1929302939392            data unit checksum
==== ========  ========================  ============================================
