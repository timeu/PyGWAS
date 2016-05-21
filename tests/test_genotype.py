import pytest
import numpy
from pygwas.core import phenotype


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


	def test_snp_iterator_chunked(self,geno):
		snps = next(geno.get_snps_iterator(2,True,chunk_size=1000))
		assert len(snps) == 1000
		snp = snps[0]
		assert sum(snp) == 1030
		assert snp[200] == 1
		assert snp[811] == 0
		
		#numpy.where(self.chrs == str(chr))[0][0]

	def test_snp_iterator_non_chunked(self,geno):
		snp = next(geno.get_snps_iterator(2,False))
		assert len(snp) == 1386
		assert sum(snp) == 1030
		assert snp[200] == 1
		assert snp[811] == 0
		

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


	def test_get_snps_filtered_by_ix(self,geno):
		index_filter = [10,500,1000]
		geno.filter_accessions_ix(index_filter)
		assert(len(geno.accessions) == 3)
		assert geno.accessions.tolist() == ['9352', '1872', '6188']
		
		
	def test_get_snps_filtered_by_id(self,geno):
		accession_filter = ['1872','9352','6188']
		geno.filter_accessions(accession_filter)
		assert(len(geno.accessions) == 3)
		assert geno.accessions.tolist() == ['9352', '1872', '6188']
		


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
    
	def test_coorindate_with_phenotype(self,geno):
		ets = ['9390','9356','6909','asdas']
		values = [1,2,3,4]
		pheno = phenotype.Phenotype(ets,values,'test')
		res = geno.coordinate_w_phenotype_data(pheno)
		assert len(geno.accessions) == 3
		assert geno.accessions.tolist() == ['9356', '6909', '9390']
		assert res['n_filtered_snps'] == 86497
		assert res['pd_indices_to_keep'] == [1,2,0]
