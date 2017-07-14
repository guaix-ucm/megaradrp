from __future__ import division
from __future__ import print_function

import argparse
from copy import deepcopy
import json
import numpy as np
from numpy.polynomial import Polynomial


def traces2ds9(json_file, ds9_file, rawimage, numpix=100, fibid_at=0,
               yoffset=0.0):
    """Transfor fiber traces from JSON to ds9-region format.

    Parameters
    ----------
    json_file : str
        Input JSON file name.
    ds9_file : file
        Handle to output file name in ds9-region format.
    rawimage : bool
        If True the traces must be generated to be overplotted on
        raw FITS images.
    numpix : int
        Number of abscissae per fiber trace.
    fibid_at : int
        Abscissae where the fibid is shown (default=0 -> not shown).
    yoffset : float
        Vertical offset (in pixels).

    """

    # offset between polynomial and image abscissae
    if rawimage:
        ix_offset = 51
    else:
        ix_offset = 1

    # open output file and insert header

    ds9_file.write('# Region file format: DS9 version 4.1\n')
    ds9_file.write('global color=green dashlist=2 4 width=2 '
                   'font="helvetica 10 normal roman" select=1 '
                   'highlite=1 dash=1 fixed=0 edit=1 '
                   'move=1 delete=1 include=1 source=1\n')
    ds9_file.write('physical\n')

    # read traces from JSON file and save region in ds9 file
    colorbox = ['#ff77ff','#4444ff']
    bigdict = json.loads(open(json_file).read())
    insmode = bigdict['tags']['insmode']
    ds9_file.write('#\n# insmode: {0}\n'.format(insmode))
    vph = bigdict['tags']['vph']
    ds9_file.write('# vph: {0}\n'.format(vph))
    uuid = bigdict['uuid']
    ds9_file.write('# uuid: {0}\n'.format(uuid))
    for fiberdict in bigdict['contents']:
        fibid = fiberdict['fibid']
        boxid = fiberdict['boxid']
        xmin = fiberdict['start']
        xmax = fiberdict['stop']
        coeff = np.array(fiberdict['fitparms'])
        ds9_file.write('#\n# fibid: {0}\n'.format(fibid))
        # skip fibers without trace
        if len(coeff) > 0:
            xp = np.linspace(start=xmin, stop=xmax, num=numpix)
            ypol = Polynomial(coeff)
            yp = ypol(xp)
            if rawimage:
                lcut = (yp > 2056.5)
                yp[lcut] += 100
            for i in range(len(xp)-1):
                x1 = xp[i] + ix_offset
                y1 = yp[i] + 1 + yoffset
                x2 = xp[i+1] + ix_offset
                y2 = yp[i+1] + 1 + yoffset
                ds9_file.write('line {0} {1} {2} {3}'.format(x1, y1, x2, y2))
                ds9_file.write(' # color={0}\n'.format(colorbox[boxid % 2]))
                if fibid_at != 0:
                    if x1 <= fibid_at <= x2:
                        ds9_file.write('text {0} {1} {{{2}}} # color=green '
                                       'font="helvetica 10 bold '
                                       'roman"\n'.format((x1+x2)/2,
                                                         (y1+y2)/2, fibid))
        else:
            print('Warning ---> Missing fiber:', fibid)


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='traces2ds9')
    # positional parameters
    parser.add_argument("json_file",
                        help="JSON file with fiber traces",
                        type=argparse.FileType('r'))
    parser.add_argument("ds9_file",
                        help="Output region file in ds9 format",
                        type=argparse.FileType('w'))
    # optional parameters
    parser.add_argument("--numpix",
                        help="Number of pixels/trace (default 100)",
                        default=100, type=int)
    parser.add_argument("--yoffset",
                        help="Vertical offset (+upwards, -downwards)",
                        default=0, type=float)
    parser.add_argument("--new_json",
                        help="New JSON file after applying specified yoffset",
                        type=argparse.FileType('w'))
    parser.add_argument("--rawimage",
                        help="FITS file is a RAW image (RSS assumed instead)",
                        action="store_true")
    parser.add_argument("--fibid_at",
                        help="Display fiber identification number at location",
                        default=0, type=int)

    args = parser.parse_args(args=args)

    traces2ds9(args.json_file.name, args.ds9_file, args.rawimage,
               args.numpix, args.fibid_at, args.yoffset)

    if args.new_json is not None:
        bigdict = json.loads(args.json_file.read())
        newdict = deepcopy(bigdict)
        for fiber in newdict['contents']:
            if len(fiber['fitparms']) > 0:
                fiber['fitparms'][0] += args.yoffset
        json.dump(newdict, args.new_json, indent=2)


if __name__ == "__main__":

    main()
