
from setuptools import find_packages
from setuptools import setup
import numpy


setup(name='megaradrp',
      version='0.5.dev0',
      include_dirs = [numpy.get_include()],
      url='http://guaix.fis.ucm.es/hg/megaradrp',
      license='GPLv3',
      description='MEGARA Data Reduction Pipeline',
      packages=find_packages(),
      package_data={'megaradrp': ['drp.yaml', 'primary.txt']},
      install_requires=[
         'numpy',
         'astropy >= 1.0',
         'scipy',
         'numina >= 0.13.0',
         ],
      zip_safe=False,
      entry_points = {
        'numina.pipeline.1': [
            'MEGARA = megaradrp.loader:megara_drp_load',
            ],
         'numina.storage.1': [
            'megara_default = megaradrp.loader:load_cli_storage',
            ]
        },
      classifiers=[
                   "Programming Language :: Python :: 2.7",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
    long_description=open('README.rst').read()
)
