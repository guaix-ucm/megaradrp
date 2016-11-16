
from setuptools import find_packages
from setuptools import setup
import numpy


setup(
    name='megaradrp',
    version='0.6.dev1',
    author='Sergio Pascual',
    author_email='sergiopr@fis.ucm.es',
    include_dirs = [numpy.get_include()],
    url='https://github.com/guaix-ucm/megaradrp',
    license='GPLv3',
    description='MEGARA Data Reduction Pipeline',
    packages=find_packages(),
    package_data={
        'megaradrp': [
            'drp.yaml',
            'primary.txt',
            'lcb_default_header.txt'
        ],
        'megaradrp.instrument.configs': [
            'component-86a6e968-8d3d-456f-89f8-09ff0c7f7c57.json',
            'component-2e02e135-2325-47c9-9975-466b445b0b8b.json',
            'component-97d48545-3258-473e-9311-00e390999d52.json',
            'instrument-4fd05b24-2ed9-457b-b563-a3c618bb1d4c.json',
            'instrument-9a86b2b2-3f7d-48ec-8f4f-3780ec967c90.json',
            'instrument-66f2283e-3049-4d4b-8ef1-14d62fcb611d.json'
        ],
    },
    install_requires=[
        'numpy',
        'astropy >= 1.0',
        'scipy',
        'numina >= 0.14',
        'scikit-image',
    ],
    zip_safe=False,
    entry_points = {
        'numina.pipeline.1': [
        '   MEGARA = megaradrp.loader:load_drp',
            ],
        },
        classifiers=[
            "Programming Language :: Python :: 2.7",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            'Development Status :: 3 - Alpha',
            "Environment :: Other Environment",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: GNU General Public License (GPL)",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering :: Astronomy",
        ],
    long_description=open('README.rst').read()
)
