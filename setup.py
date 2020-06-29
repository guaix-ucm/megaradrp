
from setuptools import find_packages
from setuptools import setup


setup(
    name='megaradrp',
    version='0.10.1',
    author='Sergio Pascual',
    author_email='sergiopr@fis.ucm.es',
    url='https://github.com/guaix-ucm/megaradrp',
    license='GPLv3',
    description='MEGARA Data Reduction Pipeline',
    packages=find_packages(),
    package_data={
        'megaradrp': [
            'drp.yaml',
        ],
        'megaradrp.instrument.configs': [
            'primary.txt',
            'lcb_default_header.txt',
            'mos_default_header.txt',
            'component-*.json',
            'instrument-*.json',
            'properties-*.json',
        ],
        'megaradrp.schemas': [
            'baseimage.json',
            'basestruct.json'
        ]
    },
    install_requires=[
        'setuptools>=36.2.1',
        'numpy',
        'astropy >= 2',
        'scipy',
        'numina >= 0.22',
        'scikit-image',
        'enum34;python_version<"3.4"',
        'contextlib2;python_version<"3.5"',
        'jsonschema'
    ],
    extras_require={
        'DB': ['sqlalchemy', 'numinadb'],
        'test': ['pytest', 'pytest-remotedata']
    },
    zip_safe=False,
    entry_points={
        'numina.pipeline.1': [
            'MEGARA = megaradrp.loader:load_drp'
        ],
        'numinadb.extra.1': [
            'MEGARA = megaradrp.db [DB]'
        ],
        'console_scripts': [
            'megaradrp-overplot_traces = megaradrp.tools.overplot_traces:main',
            'megaradrp-heal_traces = megaradrp.tools.heal_traces:main',
            'megaradrp-cube = megaradrp.processing.cube:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        'Development Status :: 3 - Alpha',
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
    long_description=open('README.rst').read()
)
