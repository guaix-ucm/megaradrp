import matplotlib.pyplot as plt

import astropy.io.fits as fits
import numpy as np
import yaml

import sys


fimage = sys.argv[1]


data = fits.getdata(fimage)


#plt.plot(data[:, 1000], label='1000')
plt.plot(data[:, 2000], label='2000')
plt.plot(data[:, 2500], label='2500')
plt.plot(data[:, 3000], label='3000')

plt.legend()
plt.show()