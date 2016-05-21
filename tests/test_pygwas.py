import pytest
from pygwas import pygwas
import numpy as np

class TestPyGWAS:

	def test_statistics(self,statisticsArgs):
		stats = pygwas.calculate_stats(statisticsArgs)
		assert len(stats) == 2
		np.testing.assert_allclose(stats['shapiro'],4.432536287882552e-14,rtol=1e-5, atol=0)
		np.testing.assert_allclose(stats['pseudo_heritability'],0.89686416775572486,rtol=1e-5, atol=0)



