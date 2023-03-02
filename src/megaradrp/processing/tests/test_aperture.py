# from megaradrp.processing.aperture import ApertureExtractor
# from numina.core import DataFrame
# import shutil
# from tempfile import mkdtemp
# import numpy as np
# import astropy.io.fits as fits
#
# def test_Aperture_Extractor():
#
#     eq = 0.8 * np.ones((4112, 4096))
#     temporary_path = mkdtemp()
#
#     fits.writeto('%s/eq.fits' % temporary_path, eq, overwrite=True)
#
#     image =file(temporary_path + '/eq.fits')
#     obj = ApertureExtractor(None)
#     # result = obj._run(image)
#     # print (result)
#
#
# if __name__ == "__main__":
#     test_Aperture_Extractor()
