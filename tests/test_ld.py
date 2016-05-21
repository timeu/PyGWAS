import pytest
import numpy as np
from pygwas.core import genotype
from pygwas.core import ld

ld_check = np.asarray([[1.,0.07102647,0.22865889,0.15327988,0.04551902],[0.07102647,1.,0.05625208,0.04115092,0.01109007],[ 0.22865889,  0.05625208,1.,0.57805568,0.07132011],[0.15327988,0.04115092,0.57805568,1.,0.03045595],[0.04551902,0.01109007,0.07132011,0.03045595,1.]])
ld_check_filtered = np.asarray([[ 1.        ,  0.07673873,  0.21612092,  0.15915561,  0.0410255 ],
       [ 0.07673873,  1.        ,  0.04983246,  0.04086861,  0.00951472],
       [ 0.21612092,  0.04983246,  1.        ,  0.62678773,  0.05123246],
       [ 0.15915561,  0.04086861,  0.62678773,  1.        ,  0.02648462],
       [ 0.0410255 ,  0.00951472,  0.05123246,  0.02648462,  1.        ]])
chr_pos = [(1,6768),(1,6514),(1,6063),(1,6449),(1,6603)]

accession_filter = [8239,8369,6970,6064,6943,6961,8374,7514,8264,8274,8249,7520,9058,6973,7000,8387,7519,6965,6956,6958,8337,6969,100000,8354,7525,7515,8285,8420,6963]

class TestLD:
	
    def test_calculate_ld(self,geno):
        ld_data = geno.calculate_ld(chr_pos)
        np.testing.assert_allclose(ld_data,ld_check,rtol=1e-5, atol=0)
        

    def test_calculate_ld_from_filtered(self,geno):
        accession_ix = range(0,1200)
        geno.filter_accessions_ix(accession_ix)
        ld_data = geno.calculate_ld(chr_pos)
        np.testing.assert_allclose(ld_data,ld_check_filtered,rtol=1e-5, atol=0)



    def test_get_ld_for_snp(self,gwas_filename):
        ld_for_snp = ld.get_ld_for_snp(gwas_filename,'3',128748)
        assert len(ld_for_snp) == 2
        assert len(ld_for_snp['chr2']['snps']) == 47
        assert len(ld_for_snp['chr2']['snps']) == len(ld_for_snp['chr2']['r2'])
        
        assert ld_for_snp['chr2']['snps'][20] == 19662552
        np.testing.assert_allclose(ld_for_snp['chr2']['r2'][20],0.0046142023)
        
        assert len(ld_for_snp['chr3']['snps']) == 353
        assert len(ld_for_snp['chr3']['snps']) == len(ld_for_snp['chr3']['r2'])
        
        assert ld_for_snp['chr3']['snps'][150] == 121389
        np.testing.assert_allclose(ld_for_snp['chr3']['r2'][150],0.010669384)
        
    def test_get_ld_for_region(self,gwas_filename):
        ld_for_region = ld.get_ld_for_region(gwas_filename,'3',104414,240894)
        assert 'snps' in ld_for_region
        assert 'r2' in ld_for_region
        assert len(ld_for_region['snps']) == 184
        assert len(ld_for_region['snps']) == len(ld_for_region['r2'])
        
        assert ld_for_region['snps'][0] == 104414
        assert ld_for_region['snps'][90] == 189023
        assert ld_for_region['snps'][183] == 240894
        np.testing.assert_allclose(ld_for_region['r2'][90][0],0.05341605097055435,)
        np.testing.assert_allclose(ld_for_region['r2'][90][90],1.0)
        
    def test_calculate_ld_for_region_unfiltered(self,geno):
        accessions = []
        ld_data = ld.calculate_ld_for_region(geno,accessions,'3',128748,10)
        assert len(ld_data['snps']) == 20
        assert len(ld_data['r2']) == len(ld_data['snps'])
        assert ld_data['start'] == 123274
        assert ld_data['start'] == ld_data['snps'][0] 
        assert ld_data['end'] == 134344
        assert ld_data['end'] == ld_data['snps'][-1]
        assert len(ld_data['r2'][0]) == 1
        assert len(ld_data['r2'][-1]) == 20
        np.testing.assert_allclose(ld_data['r2'][10][10],1.0)
        np.testing.assert_allclose(ld_data['r2'][10][0],0.04109880944844489)
        
    # skip because nan should be there and later handled    
    #def test_calculate_ld_for_region_filtered_noNaNs(self,geno):
     #   ld_data = ld.calculate_ld_for_region(geno,accession_filter,'1',18922049,10)
      #  assert np.isnan(ld_data['r2'][5]).any() == False
                
        
        