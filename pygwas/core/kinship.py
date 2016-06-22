"""
A collection of useful functions for manipulating and calculating kinship matrices.
"""
import scipy as sp
import sys
import os
import h5py
import logging
import itertools as it

log = logging.getLogger(__name__)

def calc_ibs_kinship(genotype, snp_dtype='int8', dtype='single',chunk_size=None):
    """
    Calculates IBS kinship
    
    data_format: two are currently supported, 'binary', and 'diploid_int'
    """
    num_snps = genotype.num_snps
    num_lines = len(genotype.accessions)
    if chunk_size == None:
        chunk_size = num_lines
    k_mat = sp.zeros((num_lines, num_lines), dtype=dtype)
    log.info('Starting calculation of IBS kinship')
    chunk_i = 0
    snps = genotype.get_snps_iterator(is_chunked=True,chunk_size=chunk_size)
    snps_data_format = genotype.data_format
    for snps_chunk in snps:
        chunk_i += 1
        snps_array = sp.array(snps_chunk, dtype=snp_dtype)
        snps_array = snps_array.T
        if snps_data_format == 'diploid_int':
            for i in range(num_lines):
                for j in range(i):
                    bin_counts = sp.bincount(sp.absolute(snps_array[j] - snps_array[i]))
                    if len(bin_counts) > 1:
                        k_mat[i, j] += (bin_counts[0] + 0.5 * bin_counts[1])
                    else:
                        k_mat[i, j] += bin_counts[0]
                    k_mat[j, i] = k_mat[i, j]
        elif snps_data_format == 'binary':
            sm = sp.mat(snps_array * 2.0 - 1.0)
            k_mat = k_mat + sm * sm.T
        else:
            raise NotImplementedError
        log.debug('%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * chunk_size) / num_snps))))
    if snps_data_format == 'diploid_int':
        k_mat = k_mat / float(num_snps) + sp.eye(num_lines)
    elif snps_data_format == 'binary':
        k_mat = k_mat / (2 * float(num_snps)) + 0.5
    log.info('Finished calculation')
    return k_mat


def calc_ibd_kinship(genotype, dtype='single',chunk_size=None):
    num_snps = genotype.num_snps
    num_lines = len(genotype.accessions)
    if chunk_size == None:
        chunk_size = num_lines
    k_mat = sp.zeros((num_lines, num_lines), dtype=dtype)
    log.info('Starting calculation of IBD Kinship')
    chunk_i = 0
    snps = genotype.get_snps_iterator(is_chunked=True,chunk_size=chunk_size)
    for snps_chunk in snps:
        chunk_i += 1
        snps_array = sp.array(snps_chunk)
        snps_array = snps_array.T
        norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
        x = sp.mat(norm_snps_array.T)
        k_mat += x.T * x
        log.debug('%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * chunk_size) / num_snps))))
    k_mat = k_mat / float(num_snps)
    log.info('Finished calculation')
    return k_mat



def prepare_k(k, k_accessions, accessions):
    if k_accessions == accessions:
        return sp.mat(k)
    indices_to_keep = []
    for acc in accessions:
        try:
            i = k_accessions.index(acc)
            indices_to_keep.append(i)
        except:
            continue
    k = k[indices_to_keep, :][:, indices_to_keep]
    return sp.mat(k)



def scale_k(k):
    c = sp.sum((sp.eye(len(k)) - (1.0 / len(k)) * sp.ones(k.shape)) * sp.array(k))
    scalar = (len(k) - 1) / c
    log.debug('Kinship scaled by: %0.4f' % scalar)
    k = scalar * k
    return k

def update_kinship(self, removed_snps, full_kinship, full_indivs, full_num_snps, retained_indivs, kinship_type='ibs',
                   snps_data_format='binary', snp_dtype='int8', dtype='single'):
    assert kinship_type == 'ibs', 'Only IBS kinships can be updated at the moment'
    #Cut full kinship
    cut_kinship = prepare_k(full_kinship, full_indivs, retained_indivs)
    num_lines = cut_kinship.shape[0]
    k_mat = sp.zeros((num_lines, num_lines), dtype=dtype)
    num_snps = len(removed_snps)
    snps_array = sp.array(removed_snps, dtype=snp_dtype)
    snps_array = snps_array.T
    if snps_data_format == 'diploid_int':
        for i in range(num_lines):
            for j in range(i):
                bin_counts = sp.bincount(sp.absolute(snps_array[j] - snps_array[i]))
                if len(bin_counts) > 1:
                    k_mat[i, j] += (bin_counts[0] + 0.5 * bin_counts[1])
                else:
                    k_mat[i, j] += bin_counts[0]
                k_mat[j, i] = k_mat[i, j]
    elif snps_data_format == 'binary':
        sm = sp.mat(snps_array * 2.0 - 1.0)
        k_mat = k_mat + sm * sm.T
    else:
        raise NotImplementedError
    if self.data_format == 'diploid_int':
        k_mat = k_mat / float(num_snps) + sp.eye(num_lines)
    elif self.data_format == 'binary':
        k_mat = k_mat / (2 * float(num_snps)) + 0.5

    updated_k = (cut_kinship * full_num_snps - k_mat * removed_snps) / (full_num_snps - removed_snps)
    return updated_k


def update_k_monomorphic(n_removed_snps, K, full_num_snps,dtype='single'):
    num_lines = K.shape[0]
    m = sp.ones((num_lines, num_lines), dtype=dtype) * n_removed_snps
    updated_k = (K * full_num_snps - m) / (full_num_snps - n_removed_snps)
    return updated_k





def load_kinship_from_file(kinship_file, accessions=None, scaled=True,n_removed_snps=None):
    assert os.path.isfile(kinship_file), 'File not found.'
    f = h5py.File(kinship_file,'r')
    k = f['kinship'][...]
    k_accessions = list(f['accessions'][...])
    n_snps = int(f['n_snps'][...])
    f.close()
    if accessions  != None:
        k = prepare_k(k, k_accessions, accessions)
    if n_removed_snps != None:
        k = update_k_monomorphic(n_removed_snps,k,n_snps)
    if scaled:
        k = scale_k(k)
    return {'k':k, 'accessions':k_accessions, 'n_snps':n_snps}


def save_kinship_to_file(kinship_file, kinship_mat, k_accessions, n_snps):
    f = h5py.File(kinship_file, 'w')
    f.create_dataset('kinship', data=kinship_mat)
    f.create_dataset('accessions', data=k_accessions)
    f.create_dataset('n_snps', data=n_snps)
    f.close()


def save_kinship_in_text_format(filename, k, accessions):
    with open(filename, 'w') as f:
	writer = csv.writer(f,delimiter=',')
        for acc, row in it.izip(accessions, k):
            f.writerow([acc,row.tolist()])





