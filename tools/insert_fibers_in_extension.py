
from __future__ import print_function

import sys
import pkgutil
from six import StringIO

import astropy.io.fits as fits


def create_fibers(slit):

    if slit == 'LCB':
        slith = 'lcb_default_header.txt'
    elif slit == 'MOS':
        slith= 'mos_default_header.txt'
    else:
        raise ValueError('Wrong slit %s' % (slit,))

    data = pkgutil.get_data('megaradrp.instrument.configs', slith)
    default_hdr = StringIO(data)
    hdr_fiber = fits.header.Header.fromfile(default_hdr)
    fibers = fits.ImageHDU(header=hdr_fiber, name='FIBERS')
    return fibers


def insert_fibers(fname, ext):

    with fits.open(fname) as hdul:
        while "fibers" in hdul:
            del hdul["fibers"]
        nhdul_l = [l for l in hdul]
        nhdul_l.append(ext)
        thdulist = fits.HDUList(nhdul_l)
        thdulist.writeto(fname, clobber=True)


if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Insert FIBERS extension')
    parser.add_argument('filename', metavar='FILE', nargs='+',
                        help='files to insert the FIBERS extension')

    parser.add_argument('--slit', default='LCB', choices=['MOS', 'LCB'], help='insert FIBERS')

    args = parser.parse_args()

    fibers = create_fibers(args.slit)

    for fname in args.filename:
        try:
            print(fname)
            insert_fibers(fname, fibers)
        except IOError as error:
            print(error, file=sys.stderr)
