
from setuptools import find_packages
from setuptools import setup, Extension
import numpy
# try to handle gracefully Cython
try:
    from Cython.Distutils import build_ext
    ext1 = Extension('megaradrp.trace._traces',
                     ['megaradrp/trace/traces.pyx',
                      'megaradrp/trace/Trace.cpp'],
                     language='c++')
    cmdclass = {'build_ext': build_ext}
except ImportError:
    print('We do not have Cython, just using the generated files')
    ext1 = Extension('megaradrp.trace._traces',
                     ['megaradrp/trace/traces.cpp',
                      'megaradrp/trace/Trace.cpp'],
                     language='c++')
    cmdclass = {}

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
      ext_modules=[ext1],
      cmdclass=cmdclass,
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
