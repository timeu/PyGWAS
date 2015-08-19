import pytest
from os import path
from pygwas.core import genotype

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",help="run slow tests")

def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getoption("--runslow"):
        pytest.skip("need --runslow option to run")

slow = pytest.mark.slow
resource_path = path.join(path.dirname(__file__), 'res')

@pytest.fixture(scope='module',params=["hdf5", slow("csv")])
def geno(request):

    if request.param == 'hdf5':
        return genotype.load_hdf5_genotype_data('%s/snps.hdf5' %resource_path)
    elif request.param == 'csv':
        return genotype.load_csv_genotype_data('%s/snps.csv' % resource_path)
    else:
        raise Exception('not supported')

