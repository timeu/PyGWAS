import pytest
from os import path
from pygwas.core import genotype
from pygwas.core import phenotype

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",help="run slow tests")

def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getoption("--runslow"):
        pytest.skip("need --runslow option to run")

slow = pytest.mark.slow
resource_path = path.join(path.dirname(__file__), 'res')

@pytest.fixture
def geno():
    return genotype.load_hdf5_genotype_data('%s/all_chromosomes_binary.hdf5' %resource_path)


@pytest.fixture
def kinship():
    return None


@pytest.fixture
def phenotype():
    return phenotype.parse_phenotype_file('%s/phenotype.csv' % resource_path)

@pytest.fixture
def statisticsArgs():
    return {'genotype_folder':resource_path,'genotype':'','type':'all'
        ,'file':'%s/phenotype.csv' %resource_path,'kinship':None}
