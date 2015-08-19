import pytest
from pygwas import pygwas

class TestPyGWAS:

	def test_statistics(self,statisticsArgs):
		stats = pygwas.calculate_stats(statisticsArgs)
		assert len(stats) == 2
		assert stats['shapiro'] == 4.432536287882552e-14
		assert stats['pseudo_heritability'] == 0.89686416775572486



