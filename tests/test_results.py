import pytest
from pygwas.core import result
import numpy as np

class TestResult:

    def test_load_csv_pvalues(self, csv_pvalues):
        res = result.load_from_csv(csv_pvalues)
        self._assert_result(res)

    def test_load_csv_scores(self, csv_scores):
        res = result.load_from_csv(csv_scores)
        self._assert_result(res)

    def _assert_result(self, res):
        pvals_to_check = [2.2106989795385197e-08, 0.791800358928746, 0.731829069116948, 0.557893156575617, 0.490002833429499, 0.590525250208836, 1.0, 0.590525250208836, 0.233147318752054, 0.504440789460337]
        assert res.chromosomes == ['1', '1', '2', '2', '3', '3', '4', '4', '5', '5']
        assert res.chrs == ['1', '2', '3', '4', '5']
        np.testing.assert_allclose(res.min_pval, 2.2106989795385197e-08,rtol=1e-5, atol=0)
        assert res.bonferroni_threshold == 2.3010299956639813
        assert res.stats['med_pval'] == 0.59052525020883595
        assert res.stats['bh_thres_d']['thes_pval'] == 0.010000000000000002
        assert res.stats['ks_stats']['p_val'] == 0.30742245590503603
        assert res.stats['ks_stats']['D'] == 0.29000283342949901
        assert res.maf_dict['mafs'] == [0.0028328611898017, 0.0084985835694051, 0.0028328611898017, 0.0028328611898017, 0.0056657223796034, 0.141643059490085, 0.0028328611898017, 0.141643059490085, 0.113314447592068, 0.0028328611898017]
        assert res.maf_dict['macs'] == [1, 3, 1, 1, 2, 50, 1, 50, 40, 1]
        assert res.positions == [100003, 1000033, 1000084, 1000091, 1000114, 100013, 1000267, 100027, 1000383, 1000386]
        for i, pval in enumerate(res.pvals):
            np.testing.assert_allclose(pval, pvals_to_check[i],rtol=1e-5, atol=0)

