import pytest
from os import path
from pygwas.core import genotype
from pygwas.core import phenotype
from pygwas.core import result
import json

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true",help="run slow tests")

def pytest_runtest_setup(item):
    if 'slow' in item.keywords and not item.config.getoption("--runslow"):
        pytest.skip("need --runslow option to run")

slow = pytest.mark.slow
resource_path = path.join(path.dirname(__file__), 'res')
genes_for_enrichment = []
with open(path.join(resource_path,'genes.json'),'r') as fp:
    genes_for_enrichment = json.load(fp)

gwas_result = result.load_from_hdf5(path.join(resource_path,'gwas.hdf5'))

snps_for_enrichment = gwas_result.get_top_snps(1000)
snps_for_enrichment.sort(order=['scores'])
snps_for_enrichment = snps_for_enrichment[:1000]

@pytest.fixture
def geno():
    return genotype.load_hdf5_genotype_data('%s/all_chromosomes_binary.hdf5' %resource_path)


@pytest.fixture
def kinship():phenotype
    return None

@pytest.fixture
def phenotype():
    return phenotype.parse_phenotype_file('%s/phenotype.csv' % resource_path)

@pytest.fixture
def statisticsArgs():
    return {'genotype_folder':resource_path,'genotype':'','type':'all'
        ,'file':'%s/phenotype.csv' %resource_path,'kinship':None}

@pytest.fixture
def ld_filename():
    return '%s/ld.hdf5' % resource_path


@pytest.fixture
def genes():
    return {'gene_set':genes_for_enrichment[:]}


@pytest.fixture
def top_snps():
    return snps_for_enrichment[:]
