
from setuptools import find_packages, setup

BASE_PKGS=['megara', 'megara.drp', 'megara.drp.recipes']
NAMESPACE_PKGS = ['numina.pipelines', 'numina.pipelines.megara']
ALL_PKGS = BASE_PKGS + NAMESPACE_PKGS

setup(name='megaradrp',
      version='0.3.0',
      author='Sergio Pascual',
      author_email='sergiopr@fis.ucm.es',
      url='http://guaix.fis.ucm.es/projects/megara',
      license='GPLv3',
      description='MEGARA Data Reduction Pipeline',
      packages=ALL_PKGS,
      package_dir={'megara': 'src/megara', 'numina.pipelines': 'src/drp'},
      package_data={'megara.drp': ['drp.yaml', 'primary.txt']},
      install_requires=['numpy', 'astropy', 'scipy', 'numina>=0.12.0'],
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
