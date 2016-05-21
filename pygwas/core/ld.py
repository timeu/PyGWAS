import numpy as np
import os
import h5py
from operator import itemgetter

def get_ld_for_snp(filename,chromosome,position):
    '''
    Retrieves r2 values for a specific SNP
    '''
    if not os.path.exists(filename):
        raise ValueError('file %s does not exist' % filename)
    f = h5py.File(filename,'r')
    try:
        if not 'ld' in f:
            raise ValueError('GWAS result does not contain ld dataset')
        g = f['ld']
        chr_region = g['ld_snps'].attrs['chr%s'%chromosome]
        snps = g['ld_snps'][:]
        snp_ix = chr_region[0] + np.searchsorted(snps[chr_region[0]:chr_region[1]],position)
        if snps[snp_ix] != position:
            raise ValueError('SNP %s can not be found' % position)
        ld_data = g['ld_data']
        r2_values_before_snp = ld_data[snp_ix]
        r2_values_after_snp =  []
        for row in ld_data[snp_ix+1:]:
            r2_values_after_snp.append(row[snp_ix]) 
        r2_values_after_snp = np.asarray(r2_values_after_snp)
        r2_values = np.concatenate([r2_values_before_snp,r2_values_after_snp])
        data = {}
        for chr in g['ld_snps'].attrs:
            chr_region = g['ld_snps'].attrs[chr]
            data[chr] = {'r2':r2_values[chr_region[0]:chr_region[1]].tolist(),'snps':snps[chr_region[0]:chr_region[1]].tolist()}
        return data
    finally:
        f.close()
        
def get_ld_for_region(filename,chromosome,start_pos,end_pos):
    '''
    Retrieves paire-wise r2 values for a region
    '''
    if not os.path.exists(filename):
        raise ValueError('file %s does not exist' % filename)
    f = h5py.File(filename,'r')
    try:
        if not 'ld' in f:
            raise ValueError('GWAS result does not contain ld dataset')
        g = f['ld']
        ld_data = g['ld_data']
        chr_region = g['ld_snps'].attrs['chr%s'%chromosome]
        snps = g['ld_snps'][chr_region[0]:chr_region[1]]
        start_ix = np.searchsorted(snps,start_pos)
        end_ix = np.searchsorted(snps,end_pos)+1 
        abs_start_ix = chr_region[0] + start_ix
        r2_values = [r2val[abs_start_ix:].tolist() for r2val in ld_data[abs_start_ix:chr_region[0] + end_ix]]
        return {'snps':snps[start_ix:end_ix].tolist(),'r2':r2_values}
    finally:
        f.close()
        
        
def calculate_ld_for_region(genotypeData,accessions,chromosome,position,num_snps=250):
    '''
    Calculates LD for a specific region
    '''
    genotypeData.filter_accessions(accessions)
    chr_region = genotypeData.chr_regions[genotypeData.get_chr_region_ix(chromosome)]
    abs_ix,ix,found = genotypeData.get_pos_ix(chromosome,position)
    min_ix = max(chr_region[0],abs_ix - num_snps)
    max_ix = min(abs_ix + num_snps, chr_region[1])
    positions = genotypeData.positions[min_ix:max_ix]
    chr_pos_list = zip([chromosome]*(max_ix-min_ix),positions)
    chr_pos_list = sorted(chr_pos_list,key=itemgetter(0,1))
    ld_data = genotypeData.calculate_ld(chr_pos_list)
    r2_values = []
    for i in range(0, len(ld_data)):
        r2_values.append(ld_data[i][:i + 1].tolist())
    return {'snps':positions.tolist(),'r2':r2_values,'start':int(positions[0]),'end':int(positions[-1])}