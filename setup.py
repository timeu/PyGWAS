from distutils.core import setup

setup(
    name='PyGWAS',
    version='0.1.0',
    author='Uemit Seren',
    author_email='uemit.seren@gmi.oeaw.ac.at',
    packages=['pygwas','pygwas.core'],
    url='http://pypi.python.org/pypi/pygwas/',
    license='LICENSE.txt',
    description='GWAS library',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >=1.6.1",
        "scipy >=0.13.0",
        "h5py >=2.1.3",
        ],
    )
