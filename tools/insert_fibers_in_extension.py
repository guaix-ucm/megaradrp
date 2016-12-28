import sys
import pkgutil
from StringIO import StringIO

import astropy.io.fits as fits

def create_fibers(fname, slit='LCB'):

    if slit == 'LCB':
        slith = 'lcb_default_header.txt'
    elif slit == 'MOS':
        slith= 'mos_default_header.txt'
    else:
        raise ValueError('Wrong slit %s' % (slit,))

    with fits.open(fname) as hdul:
        while "fibers" in hdul:
            del hdul["fibers"]
        data = pkgutil.get_data('megaradrp', slith)
        default_hdr = StringIO(data)
        hdr_fiber = fits.header.Header.fromfile(default_hdr)
        fibers = fits.ImageHDU(header=hdr_fiber, name='FIBERS')
        nhdul_l = [l for l in hdul]
        nhdul_l.append(fibers)
        thdulist = fits.HDUList(nhdul_l)
        thdulist.writeto(fname, clobber=True)
#        hdul.append(fibers)
#        hdul.flush()


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Insert FIBERS extension')
    for fname in sys.argv[1:]:
        try:
            print fname
            create_fibers(fname, 'MOS')
        except IOError as error:
            print error
