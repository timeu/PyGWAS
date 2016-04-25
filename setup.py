from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path
import pygwas

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='PyGWAS',
    version=pygwas.__version__,
    description='A GWAS library',
    long_description=long_description,
    url='https://github.com/pypa/pygwas',
    author='Uemit Seren',
    author_email='uemit.seren@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
	'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],
    keywords='GWAS',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=[
        "numpy",
        "scipy",
        "h5py",
        "matplotlib >= 1.4.3"
    ],
    entry_points={
        'console_scripts': [
            'pygwas=pygwas.__main__:main'
        ],
    },
)

