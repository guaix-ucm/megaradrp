
from setuptools import find_packages
from setuptools import setup
import numpy


setup(name='megaradrp',
      version='0.5',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      include_dirs = [numpy.get_include()],
      url='https://github.com/guaix-ucm/megaradrp',
      license='GPLv3',
      description='MEGARA Data Reduction Pipeline',
      packages=find_packages(),
      package_data={'megaradrp': ['drp.yaml', 'primary.txt']},
      install_requires=[
         'numpy',
         'astropy >= 1.0',
         'scipy',
         'numina >= 0.14',
         ],
      zip_safe=False,
      entry_points = {
        'numina.pipeline.1': [
            'MEGARA = megaradrp.loader:load_drp',
            ],
        },
      classifiers=[
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3.3",
                   "Programming Language :: Python :: 3.4",
                   'Development Status :: 3 - Alpha',
                   "Environment :: Other Environment",
                   "Intended Audience :: Science/Research",
                   "License :: OSI Approved :: GNU General Public License (GPL)",
                   "Operating System :: OS Independent",
                   "Topic :: Scientific/Engineering :: Astronomy",
                   ],
    long_description=open('README.rst').read()
)
