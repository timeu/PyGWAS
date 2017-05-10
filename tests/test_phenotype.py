import pytest
import numpy as np
from pygwas.core import phenotype

class TestPhenotype:

    def test_avg(self,small_phenotype):
        assert len(small_phenotype.ecotypes) == 16
        assert len(small_phenotype.values) == 16
        assert small_phenotype.has_replicates() == True
        small_phenotype.convert_to_averages()
        assert len(small_phenotype.ecotypes) == 4
        assert len(small_phenotype.values) == 4
        assert small_phenotype.ecotypes == [1,2,3,4]
        assert small_phenotype.values == [11.5, 21.75, 61.25, 85.25]
        assert small_phenotype.has_replicates() == False


    def test_replicates_are_converted_to_avg_before_transform(self,small_phenotype):
        small_phenotype.transform('log')
        assert len(small_phenotype.ecotypes) == 4
        assert len(small_phenotype.values) == 4
        assert small_phenotype.ecotypes == [1,2,3,4]
        assert small_phenotype.values == [4.484955974761633, 4.5943422435840375, 4.930314552983338, 5.090204331913628]

    def test_log_transformation(self,phenotype):
        assert phenotype.transformation is None
        phenotype.transform('log')
        assert phenotype.transformation == 'log'
        np.testing.assert_allclose(np.mean(phenotype.values),4.415795438044978,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),0.2659027593295565,rtol=0.05, atol=0)

    def test_sqrt_transformation(self,phenotype):
        assert phenotype.transformation is None
        phenotype.transform('sqrt')
        assert phenotype.transformation == 'sqrt'
        np.testing.assert_allclose(np.mean(phenotype.values),9.1766423930810301,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),1.2054970137664522,rtol=0.05, atol=0)

    def test_sqr_transformation(self,phenotype):
        assert phenotype.transformation is None
        phenotype.transform('sqr')
        assert phenotype.transformation == 'sqr'
        np.testing.assert_allclose(np.mean(phenotype.values),7829.101397792874,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),3889.6765923818066,rtol=0.05, atol=0)


    def test_ascombe_transformation(self,phenotype):
        assert phenotype.transformation is None
        phenotype.transform('ascombe')
        assert phenotype.transformation == 'ascombe'
        np.testing.assert_allclose(np.mean(phenotype.values),17.476328580491042,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),2.533952678322978,rtol=0.05, atol=0)


    def test_exp_transformation(self,phenotype):
        assert phenotype.transformation is None
        phenotype.transform('exp')
        assert phenotype.transformation == 'exp'
        np.testing.assert_allclose(np.mean(phenotype.values),1.0851750688952121e+56,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),1.9200593483252192e+57,rtol=0.05, atol=0)

    def test_boxcox_transformation(self,phenotype):
        assert phenotype.transformation is None
        phenotype.transform('box_cox')
        assert phenotype.transformation == 'box_cox'
        np.testing.assert_allclose(np.mean(phenotype.values),22.209859554452517,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),3.753060024281832,rtol=0.05, atol=0)

    def test_most_normal_transformation(self,phenotype, small_phenotype):
        assert phenotype.transformation is None
        trans_type, shapiro_pvalue = phenotype.most_normal_transformation()
        assert trans_type == 'box_cox'
        np.testing.assert_allclose(shapiro_pvalue,3.551881553633118e-14,rtol=0.05, atol=0)
        assert phenotype.transformation == 'box_cox'
        np.testing.assert_allclose(np.mean(phenotype.values),22.209859554452517,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(phenotype.values),3.753060024281832,rtol=0.05, atol=0)

        assert small_phenotype.transformation is None
        trans_type, shapiro_pvalue = small_phenotype.most_normal_transformation()
        assert trans_type == 'box_cox'
        np.testing.assert_allclose(shapiro_pvalue,0.7982271313667297,rtol=0.05, atol=0)
        assert small_phenotype.transformation == 'box_cox'
        np.testing.assert_allclose(np.mean(small_phenotype.values),4.78125,rtol=0.05, atol=0)
        np.testing.assert_allclose(np.std(small_phenotype.values),0.240035804620894,rtol=0.05, atol=0)

