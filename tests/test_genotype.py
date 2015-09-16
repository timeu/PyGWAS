import pytest
import numpy


class TestGenotype:

	def test_get_snp_at(self,geno):
		snps = geno.get_snp_at(1,16036831)
		assert sum(snps) == 598
		assert snps[912] == 1

		snps = geno.get_snp_at(2,9539640)
		assert sum(snps) == 107
		assert snps[963] == 1

		snps = geno.get_snp_at(3,11283260)
		assert sum(snps) == 382
		assert snps[886] == 1

		snps = geno.get_snp_at(4,8734567)
		assert sum(snps) == 1035
		assert snps[890] == 0

		snps = geno.get_snp_at(5,14252117)
		assert sum(snps) == 306
		assert snps[727] == 1





	def test_get_snp_at_unknown_pos(self,geno):
		snps = geno.get_snp_at(1,16036830)
		assert snps is None

		snps = geno.get_snp_at(2,9539641)
		assert snps is None

		snps = geno.get_snp_at(3,11283261)
		assert snps is None

		snps = geno.get_snp_at(4,8734565)
		assert snps is None

		snps = geno.get_snp_at(5,14252114)
		assert snps is None


	def test_get_snps_from_pos_filtered(self,geno):
		accession_ix = range(0,1200)
		geno.filter_accessions_ix(accession_ix)
		chr_pos = [(1,572643),(3,11283260),(3,11283261),(1,16036830),(1,16036831),(5,14252117),(5,14252114),(4,8734565),(4,8734567),(2,9539640),(2,9539641)]

		retval = geno.get_snps_from_pos(chr_pos)
		assert len(retval) == 2
		assert len(retval[0]) == len(retval[1])
		assert len(retval[0]) == 6
		indices = retval[0]
		snps = retval[1]
		assert indices == [1000, 98313, 25105, 180388, 136915, 63848]

		assert sum(snps[0]) == 627
		assert snps[0][753] == 0

		assert sum(snps[1]) == 314
		assert snps[1][886] == 1

		assert sum(snps[2]) == 525
		assert snps[2][912] == 1

		assert sum(snps[3]) == 269
		assert snps[3][727] == 1

		assert sum(snps[4]) == 870
		assert snps[4][890] == 0

		assert sum(snps[5]) == 83
		assert snps[5][963] == 1


	def test_get_snps_from_pos(self,geno):
		chr_pos = [(1,572643),(3,11283260),(3,11283261),(1,16036830),(1,16036831),(5,14252117),(5,14252114),(4,8734565),(4,8734567),(2,9539640),(2,9539641)]
		retval = geno.get_snps_from_pos(chr_pos)
		assert len(retval) == 2
		assert len(retval[0]) == len(retval[1])
		assert len(retval[0]) == 6
		indices = retval[0]
		snps = retval[1]
		assert indices == [1000, 98313, 25105, 180388, 136915, 63848]

		assert sum(snps[0]) == 710
		assert snps[0][753] == 0

		assert sum(snps[1]) == 382
		assert snps[1][886] == 1

		assert sum(snps[2]) == 598
		assert snps[2][912] == 1

		assert sum(snps[3]) == 306
		assert snps[3][727] == 1

		assert sum(snps[4]) == 1035
		assert snps[4][890] == 0

		assert sum(snps[5]) == 107
		assert snps[5][963] == 1