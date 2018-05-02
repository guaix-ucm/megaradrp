import matplotlib.pyplot as plt

import astropy.io.fits as fits
import numpy as np
import yaml

import sys

fname = sys.argv[1]
fimage = sys.argv[2]


data = fits.getdata(fimage)
traces = yaml.load(open(fname))

plt.xlim([0, 4096])
plt.ylim([0, 4112])
plt.imshow(data)

for trace in traces:
    p = trace['fitparms']
    amx = np.poly1d(p)
    xx = np.arange(trace['start'], trace['stop'])
    y = amx(xx)


    plt.plot(xx, y, 'k')
plt.show()