#!/usr/bin/env python

from distutils.core import setup

setup(name='pymegara',
      version='0.1.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/~spr',
      license='GPLv3',
      description='Megara Data Reduction Pipeline',
      packages=[
                'megara', 'megara.recipes',
                ],
      install_requires=['numina'],
      classifiers=[
                   "Programming Language :: Python",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
)
