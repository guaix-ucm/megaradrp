
# Script to compute extraction weights


from __future__ import division


import copy
import numpy as np
from astropy.io import fits

import matplotlib
matplotlib.use('agg', warn=True)
import matplotlib.pyplot as plt

from astropy.modeling import fitting
from scipy.stats import norm

# Model

#from model import GaussBox, gauss_box_model

#
# Gaussian + Box
from scipy.stats import norm
from astropy.modeling.models import custom_model_1d
import math

M_SQRT_2_PI = math.sqrt(2 * math.pi)


def gauss_box_model(x, amplitude=1.0, mean=0.0, stddev=1.0, hpix=0.5):
    '''Integrate a gaussian profile.'''
    z = (x - mean) / stddev
    m2 = z + hpix / stddev
    m1 = z - hpix / stddev
    return amplitude * (norm.cdf(m2) - norm.cdf(m1))


def norm_pdf_t(x):
    return  np.exp(-0.5 *x*x) / M_SQRT_2_PI


def gauss_box_model_deriv(x, amplitude=1.0, mean=0.0, stddev=1.0, hpix=0.5):
    '''Integrate a gaussian profile.'''

    z = (x - mean) / stddev
    z2 = z + hpix / stddev
    z1 = z - hpix / stddev

    da = norm.cdf(z2) - norm.cdf(z1)

    fp2 = norm_pdf_t(z2)
    fp1 = norm_pdf_t(z1)

    dl = -amplitude / stddev * (fp2 - fp1)
    ds = -amplitude / stddev * (fp2 * z2 - fp1 * z1)
    dd = amplitude / stddev * (fp2 + fp1)

    return (da, dl, ds, dd)


GaussBox = custom_model_1d(gauss_box_model, func_fit_deriv=gauss_box_model_deriv)

## Model

def pixcont(i, x0, sig, hpix=0.5):
    '''Integrate a gaussian profile.'''
    z = (i - x0) / sig
    hpixs = hpix / sig
    z2 = z + hpixs
    z1 = z - hpixs
    return norm.cdf(z2) - norm.cdf(z1)


def g_profile(xl, l, s):
    '''A gaussian profile.'''
    z = (xl - l) / s
    return np.exp(-0.5 * z**2)


def fit1d_profile(xl, yl, init0, N, nloop=10, S=3):
    """Iterative fitting"""

    init = copy.deepcopy(init0)

    changes_a = np.zeros((N, nloop))
    changes_m = np.zeros((N, nloop))
    changes_s = np.zeros((N, nloop))

    #figure11 = plt.figure()
    #ax11 = figure11.add_subplot(111)
    #ax11.set_title('Fit')

    for il in range(nloop):

        # print 'loop', il, datetime.datetime.now()

        values = np.random.permutation(N)

        for val in values:

            m1 = max(0, int(init[val]['mean']) - 6 * S)
            m2 = int(init[val]['mean']) + 6 * S

            y = yl[m1:m2].copy()
            xt = xl[m1:m2]

            for peakid in range(max(0, val-S), min(N, val+S+1)):
                if peakid == val:
                    continue

                y -=  gauss_box_model(xt, **init[peakid])

            model = GaussBox(**init[val])
            model.mean.min = model.mean.value - 0.5
            model.mean.max = model.mean.value + 0.5
            #model.mean.fixed = True
            model.stddev.min = 1.0
            model.stddev.max = 2.0
            model.hpix.fixed = True


            #ax11.plot(xt, model(xt), 'r-')

            fitter = fitting.LevMarLSQFitter()
            model_fitted = fitter(model, xt, y)

            #ax11.plot(xt, model_fitted(xt), 'g-')
            #plt.show(block=False)

            na = model_fitted.amplitude.value
            nm = model_fitted.mean.value
            ns = model_fitted.stddev.value

            #print 'loop', il, 'peak', val
            #print 'amp', na, init2[val]['amplitude']
            #print 'mean', nm, init2[val]['mean']
            #print 'stddev', ns, init2[val]['stddev']

            changes_a[val, il] = na - init[val]['amplitude']
            changes_m[val, il] = nm - init[val]['mean']
            changes_s[val, il] = ns - init[val]['stddev']

            init[val]['amplitude'] = na
            init[val]['mean'] = nm
            init[val]['stddev'] = ns

    return init, (changes_a, changes_m, changes_s)


def calc_sparse_matrix(final, nrows, cut=1.0e-6, extra=10):
    from scipy.sparse import lil_matrix

    idxs = range(len(final))

#    g_ampl = np.array([final[i]['amplitude'] for i in idxs])
    g_mean = np.array([final[i]['mean'] for i in idxs])
    g_std = np.array([final[i]['stddev'] for i in idxs])

    # calc w
    begpix = np.ceil(g_mean - 0.5).astype('int')

    steps = np.arange(-extra, extra)
    ref = begpix + steps[:,np.newaxis]

    rr = gauss_box_model(ref, mean=g_mean, stddev=g_std)
    rrb = begpix - extra
    # Filter values below 'cut'
    rr[rr < cut] = 0.0

    # Calc Ws matrix
    block, nfib = rr.shape
    w_init = lil_matrix((nrows, nfib))

    for i in range(nfib):
        w_init[rrb[i]:rrb[i]+block, i] =  rr[:,i, np.newaxis]

    # Convert to CSR matrix
    wcol = w_init.tocsr()
    return wcol


def plot(xl, final, col, yl, changes_a, changes_m, changes_s, init_vals, centers, sigs):

    N = len(final)

    ym = np.zeros_like(xl, dtype='float')

    for i in range(N):
        yn = gauss_box_model(xl, **final[i])
        ym +=  yn

    if False:
        figure1 = plt.figure()
        ax1 = figure1.add_subplot(111)
        ax1.set_title('Values and model')
        ax1.plot(xl, yl, 'b*-')
        ax1.plot(xl, ym, 'r')

        ax12 = ax1.twinx()
        ax12.plot(xl, yl-ym, 'g')

    figure1 = plt.figure()
    ax1 = figure1.add_subplot(111)
    ax1.set_title('Values and model, col=%d' % (col,))
    ax1.plot(xl, yl, 'b*-')
    ax1.plot(xl, ym, 'r')

    figure11 = plt.figure()
    ax12 = figure11.add_subplot(111)
    ax12.plot(xl, yl-ym, 'g')

    #plt.show(block=False)

    figure2 = plt.figure()
    ax2 = figure2.add_subplot(111)
    ax2.set_title('Amplitude difference, col=%d' % (col,))
    for i in range(N):
        ax2.plot(changes_a[i], 'g*-')
    #plt.show(block=False)

    figure3 = plt.figure()
    ax3 = figure3.add_subplot(111)
    ax3.set_title('Center difference, col=%d' % (col,))
    for i in range(N):
        ax3.plot(changes_m[i], 'g*-')
    #plt.show(block=False)

    figure4 = plt.figure()
    ax4 = figure4.add_subplot(111)
    ax4.set_title('Stddev difference, col=%d' % (col,))
    for i in range(N):
        ax4.plot(changes_s[i] / 1.5, 'g*-')

    figure8 = plt.figure()
    plt.title('Amplitude output/input, col=%d' % (col,))
    plt.plot(range(N), [(final[a]['amplitude'] / init_vals[a]['amplitude']) for a in range(N)], 'b*-')

    figure9 = plt.figure()
    plt.title('Center output - input, col=%d' % (col,))
    plt.plot(range(N), [final[a]['mean']-centers[a] for a in range(N)], 'b*-')
    #plt.plot(range(N), centers, 'r*-')

    figure10 = plt.figure()
    plt.title('Stddev output/input, col=%d' % (col,))
    plt.plot(range(N), [final[a]['stddev']/sigs[a] for a in range(N)], 'b*-')
    #plt.plot(range(N), sigs, 'r*-')

    plt.show()


def plot_save(xl, final, col, yl, changes, init_vals, centers, sigs):

    changes_a, changes_m, changes_s = changes
    N = len(final)

    ym = np.zeros_like(xl, dtype='float')

    for i in range(N):
        yn = gauss_box_model(xl, **final[i])
        ym +=  yn


    figure1 = plt.figure()
    ax1 = figure1.add_subplot(111)
    ax1.set_title('Values and model, col=%d' % (col,))
    ax1.plot(xl, yl, 'b*-')
    ax1.plot(xl, ym, 'r')
    figure1.savefig('images2/fit_model-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure1)

    figure11 = plt.figure()
    ax12 = figure11.add_subplot(111)
    ax12.plot(xl, yl-ym, 'g')
    figure11.savefig('images2/fit_res-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure11)
    #plt.show(block=False)

    figure2 = plt.figure()
    ax2 = figure2.add_subplot(111)
    ax2.set_title('Amplitude difference, col=%d' % (col,))
    for i in range(N):
        ax2.plot(changes_a[i], 'g*-')
    #plt.show(block=False)
    figure2.savefig('images2/evol_ampl-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure2)

    figure3 = plt.figure()
    ax3 = figure3.add_subplot(111)
    ax3.set_title('Center difference, col=%d' % (col,))
    for i in range(N):
        ax3.plot(changes_m[i], 'g*-')
    figure3.savefig('images2/evol_center-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure3)
    #plt.show(block=False)

    figure4 = plt.figure()
    ax4 = figure4.add_subplot(111)
    ax4.set_title('Stddev difference, col=%d' % (col,))
    for i in range(N):
        ax4.plot(changes_s[i] / 1.5, 'g*-')
    figure4.savefig('images2/evol_stddev-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure4)

    figure8 = plt.figure()
    plt.title('Amplitude output/input, col=%d' % (col,))
    plt.plot(range(N), [(final[a]['amplitude'] / init_vals[a]['amplitude']) for a in range(N)], 'b*-')
    figure8.savefig('images2/io_ampl_ratio-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure8)

    figure9 = plt.figure()
    plt.title('Center output - input, col=%d' % (col,))
    plt.plot(range(N), [final[a]['mean']-centers[a] for a in range(N)], 'b*-')
    #plt.plot(range(N), centers, 'r*-')
    figure9.savefig('images2/io_center_diff-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure9)

    figure10 = plt.figure()
    plt.title('Stddev output/input, col=%d' % (col,))
    plt.plot(range(N), [final[a]['stddev']/sigs[a] for a in range(N)], 'b*-')
    #plt.plot(range(N), sigs, 'r*-')
    figure10.savefig('images2/io_stddev_ratio-%d.png' % (col,), bbox_inches='tight')
    plt.close(figure10)
    plt.show()


def calc_profile(data1, pols, col, sigma, start=0, doplots=False):
    print 'fitting column', col
    peaks = np.array([pol(col) for pol in pols])

    boxd = data1[:,col]

    centers = peaks[:]-start
    sigs = sigma * np.ones_like(centers)
    scale_sig = 0.25 # For sigma ~= 1.5, the peak is typically 0.25
    ecenters = np.ceil(centers -0.5).astype('int')

    N = len(centers)
    cmax = boxd.max()
    yl = boxd / cmax # Normalize to peak
    xl = np.arange(len(yl))

    init_vals = {}
    for i in range(N):
        init_vals[i] = {}
        init_vals[i]['amplitude'] = yl[ecenters[i]] / scale_sig
        #init_vals[i]['mean'] = ecenters[i]
        init_vals[i]['mean'] = centers[i]
        init_vals[i]['stddev'] = sigma

    final, changes = fit1d_profile(xl, yl, init_vals, N, nloop=10)

    if doplots:
        plot_save(xl, final, col, yl, changes, init_vals, centers, sigs)

    for i in range(N):
        final[i]['amplitude'] = final[i]['amplitude'] * cmax

    return final


if __name__ == '__main__':
    from multiprocessing import Pool
    import pickle
    import yaml
    # FIT the columns, store results in json final files

    data2 = fits.getdata('fiberflat_frame.fits')
    tracemap2 = yaml.load(open('master_traces.yaml'))
    pols2 = [np.poly1d(t['fitparms']) for t in tracemap2]

    start = 0
    sigma = 1.5 # Typical value for MEGARA
    nrows = data2.shape[0] # 4112
    cols = range(data2.shape[1]) # 4096

    def calc_all(col):

        import os.path
        prefix = 'fit1'

        fname = os.path.join(prefix, 'fit-%d.plk' % (col,))

        if os.path.isfile(fname):
            print 'fitting column', col, 'already done'
            return

        print 'fitting column', col
        final = calc_profile(data2, pols2, col, sigma, start=0)

        with open(fname, 'w') as fd:
            pickle.dump(final, fd)

        print 'create matrix for column', col
        wm = calc_sparse_matrix(final, nrows, cut=1.0e-6, extra=10)

        prefixw = 'plk1'
        oname = os.path.join(prefixw, 'tyu-%d.plk' % (col,))

        pickle.dump(wm, open(oname, 'w'))


    p = Pool(12)
    p.map(calc_profile, cols)



