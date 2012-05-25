#!/usr/bin/env python

from setuptools import setup, find_packages

BASE_PKGS=find_packages('src', exclude=['drp', 'drp.*'])
NAMESPACE_PKGS = ['numina.pipelines', 'numina.pipelines.megara']
ALL_PKGS = BASE_PKGS + NAMESPACE_PKGS

setup(name='pymegara',
      version='0.1.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/~spr',
      license='GPLv3',
      description='Megara Data Reduction Pipeline',
      packages=ALL_PKGS,
      package_dir={'megara': 'src/megara', 'numina.pipelines': 'src/drp'},
      package_data={'megara': ['obsmodes.yaml', 'primary.txt']},
      install_requires=['numina>=0.7.0'],
      classifiers=[
                   "Programming Language :: Python :: 2.7",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
    long_description=open('README.txt').read()
)
