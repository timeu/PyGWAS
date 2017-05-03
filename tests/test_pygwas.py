import pytest
from pygwas import pygwas
import numpy as np

class TestPyGWAS:

	def test_statistics(self,statisticsArgs):
		stats = pygwas.calculate_stats(statisticsArgs)
		assert len(stats) == 3
		assert stats['transformation'] == 'none'
		np.testing.assert_allclose(stats['shapiro'],4.432536287882552e-14,rtol=1e-5, atol=0)
		np.testing.assert_allclose(stats['pseudo_heritability'],0.89686416775572486,rtol=1e-5, atol=0)


	def test_statistics_only_shapiro(self, statisticsArgs):
		statisticsArgs['genotype_folder'] = None
		statisticsArgs['type'] = 'shapiro'
		stats = pygwas.calculate_stats(statisticsArgs)
		assert stats['transformation'] == 'none'
		assert len(stats) == 2
		np.testing.assert_allclose(stats['shapiro'],2.3893008814393317e-14,rtol=1e-5, atol=0)

	def test_statistics_only_shapiro_with_genotype(self, statisticsArgs):
		statisticsArgs['type'] = 'shapiro'
		stats = pygwas.calculate_stats(statisticsArgs)
		assert stats['transformation'] == 'none'
		assert len(stats) == 2
		np.testing.assert_allclose(stats['shapiro'],4.432536287882552e-14,rtol=1e-5, atol=0)

	def test_statistics_most_normal(self, statisticsArgs):
		statisticsArgs['transformation'] = 'most_normal'
		stats = pygwas.calculate_stats(statisticsArgs)
		assert stats['transformation'] == 'box_cox'
		assert len(stats) == 3
		np.testing.assert_allclose(stats['shapiro'],7.311389856176284e-14,rtol=1e-5, atol=0)

	def test_statistics_most_normal_only_shapiro(self, statisticsArgs):
		statisticsArgs['genotype_folder'] = None
		statisticsArgs['type'] = 'shapiro'
		statisticsArgs['transformation'] = 'most_normal'
		stats = pygwas.calculate_stats(statisticsArgs)
		assert stats['transformation'] == 'box_cox'
		assert len(stats) == 2
		np.testing.assert_allclose(stats['shapiro'],3.551881553633118e-14,rtol=1e-5, atol=0)


