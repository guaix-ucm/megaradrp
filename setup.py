
from setuptools import find_packages, setup


setup(name='megaradrp',
      version='0.4.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/hg/megaradrp',
      license='GPLv3',
      description='MEGARA Data Reduction Pipeline',
      packages=find_packages(where='src'),
      package_dir={'': 'src'},
      package_data={'megaradrp': ['drp.yaml', 'primary.txt']},
      install_requires=['numpy', 'astropy', 'scipy', 'numina>=0.13.0'],
      zip_safe=False,
      entry_points = {
        'numina.pipeline.1': [
            'megara = megaradrp.loader:megara_drp_load',
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
    long_description=open('README.txt').read()
)
