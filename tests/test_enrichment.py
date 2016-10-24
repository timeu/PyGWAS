import pytest
import numpy as np
from pygwas.core.enrichment import enrichment

class TestEnrichment:

    def test_enrichment(self,geno, genes, top_snps):
        pval = enrichment(genes,geno,top_snps,20000,10000)
        np.testing.assert_allclose(pval['gene_set'],0.055,rtol=0.05, atol=0)