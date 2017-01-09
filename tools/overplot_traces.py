from __future__ import division
from __future__ import print_function

import argparse
import json
import numpy as np
from numpy.polynomial import Polynomial

from numina.array.display.ximshow import ximshow_file
from numina.array.display.pause_debugplot import pause_debugplot


def main(args=None):
    # parse command-line options
    parser = argparse.ArgumentParser(prog='overplot_traces')
    # positional parameters
    parser.add_argument("fits_file",
                        help="FITS image containing the spectra",
                        type=argparse.FileType('r'))
    parser.add_argument("traces_file",
                        help="JSON file with traces",
                        type=argparse.FileType('r'))
    # optional parameters
    parser.add_argument("--z1z2",
                        help="tuple z1,z2, minmax or None (use zscale)")
    parser.add_argument("--bbox",
                        help="bounding box tuple: nc1,nc2,ns1,ns2")
    parser.add_argument("--keystitle",
                        help="tuple of FITS keywords.format: " +
                             "key1,key2,...keyn.'format'")
    parser.add_argument("--geometry",
                        help="tuple x,y,dx,dy")

    args = parser.parse_args(args=args)

    ax = ximshow_file(args.fits_file.name,
                      args_z1z2=args.z1z2,
                      args_bbox=args.bbox,
                      args_keystitle=args.keystitle,
                      args_geometry=args.geometry,
                      show=False)

    # read traces from JSON file
    bigdict = json.loads(open(args.traces_file.name).read())
    for fiberdict in bigdict['contents']:
        xmin = fiberdict['start']
        xmax = fiberdict['stop']
        coeff = fiberdict['fitparms']
        # skip fibers without trace
        if len(coeff) > 0:
            num = int(float(xmax-xmin+1)+0.5)
            xp = np.linspace(start=xmin, stop=xmax, num=num)
            ypol = Polynomial(coeff)
            yp = ypol(xp)
            ax.plot(xp+1, yp+1, 'b-')

    import matplotlib.pyplot as plt
    plt.show(block=False)
    plt.pause(0.001)
    pause_debugplot(12)

if __name__ == "__main__":

    main()
