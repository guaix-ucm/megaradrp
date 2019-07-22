
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy.ndimage.filters as filt
from scipy.interpolate import LSQUnivariateSpline


conf = {"trim1": [[0,2056],[50,4146]],
     "trim2": [[2156,4212],[50,4146]],
     "overscan1": [[0,2056],[4146,4196]],
     "overscan1_corner": [[2056,2106],[4146,4196]],
     "overscan2": [[2156,4212],[0,50]],
     "overscan2_corner": [[2106,2156],[4146,4196]],
     "prescan1": [[0,2056],[0,50]],
     "prescan2": [[2156,4212],[4146,4196]],
     "middle1": [[2056,2106],[50,4146]],
     "middle2": [[2106,2156],[50,4146]],
     "gain1": 1.73,
     "gain2": 1.6
}


def to_slice(sec):
    sec1, sec2 = sec
    return slice(*sec1), slice(*sec2)


def to_str(sec, format='fits'):
    sec1, sec2 = sec
    return str([to_index(sec2), to_index(sec1)])


def to_index(ssec, format='fits'):
    a, b = ssec
    return [a + 1, b]


def to_pix(sec, axis=0):
    return np.arange(*sec[axis])


def plot2(data, name, knots, ax, s=0, r=-1):
    reg = conf[name]
    axis_u = 0
    axis_v = 1
    plot4(data, reg, axis_u, axis_v, knots, ax, s=s, r=r)


def plot3(data, name, knots, ax, s=0, r=-1):
    reg = conf[name]
    axis_u = 1
    axis_v = 0
    plot4(data, reg, axis_u, axis_v, knots, ax, s=s, r=r)


def plot4(data, reg, axis_u, axis_v, knots, ax, s=0, r=-1):

    u = to_pix(reg, axis=axis_u)
    region = to_slice(reg)
    v = data[region].mean(axis=axis_v)
    v = filt.median_filter(v, size=7)

    spl2 = LSQUnivariateSpline(u, v, knots, k=3)
    v_spl2 = spl2(u)

    ax.plot(u[s:r], v[s:r], label="D")
    ax.plot(u[s:r], v_spl2[s:r], label="S2")


if __name__ == '__main__':
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Check overscan')
    parser.add_argument('filename', metavar='FILE', nargs='+',
                        help='Check overscan')


    args = parser.parse_args()


    for fname in args.filename:
        print(fname)
        fname_base, ext = os.path.splitext(fname)

        hdulist = fits.open(fname)
        hdu = hdulist[0]
        data = hdu.data
        s = 0
        r = -1
        knots1 = [10, 100, 200, 300, 400, 500, 750, 1000, 1200, 1500, 1700, 2000]
        knots2 = [2200, 3000, 3500, 3900, 4000, 4100, 4200]

        knots1 = [125, 250, 375, 500, 1000, 1500]
        knots2 = [2500, 3000, 3500, 3600, 3700, 3800]

        knots1 = [1200]
        knots2 = [3100]
        knotsm = [2100]

        for regname, knots in zip(['overscan1', 'overscan2'], [knots1, knots2]):
            fig, axes = plt.subplots(1, 1)
            regs = to_str(conf[regname])
            axes.set_title("{}\n{} {}".format(fname, regname, regs))
            plot2(data, regname, knots, axes, s=s, r=r)
            plt.savefig("{}_{}.png".format(fname_base, regname))
            plt.close()

        for regname, knots in zip(['middle1', 'middle2'], [knotsm, knotsm]):
            fig, axes = plt.subplots(1, 1)
            regs = to_str(conf[regname])
            axes.set_title("{}\n{} {}".format(fname, regname, regs))
            plot3(data, regname, knots, axes, s=s, r=r)
            plt.savefig("{}_{}.png".format(fname_base, regname))
            plt.close()
