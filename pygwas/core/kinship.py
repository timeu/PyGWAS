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

def calc_ibs_kinship(snps, snps_data_format='binary', snp_dtype='int8', dtype='single',
                     chunk_size=None):
    """
    Calculates IBS kinship
    
    data_format: two are currently supported, 'binary', and 'diploid_int'
    """
    num_snps = len(snps)
    log.debug('Allocating K matrix')
    num_lines = len(snps[0])
    if chunk_size == None:
        chunk_size = num_lines
    k_mat = sp.zeros((num_lines, num_lines), dtype=dtype)
    log.debug('Starting calculation')
    chunk_i = 0
    for snp_i in range(0, num_snps, chunk_size): #FINISH!!!
        chunk_i += 1
        snps_array = sp.array(snps[snp_i:snp_i + chunk_size], dtype=snp_dtype)
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
        sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * chunk_size) / num_snps))))
        sys.stdout.flush()
    if snps_data_format == 'diploid_int':
        k_mat = k_mat / float(num_snps) + sp.eye(num_lines)
    elif snps_data_format == 'binary':
        k_mat = k_mat / (2 * float(num_snps)) + 0.5
    return k_mat


def calc_ibd_kinship(snps, dtype='single'):
    num_snps = len(snps)
    n_indivs = len(snps[0])
    k_mat = sp.zeros((n_indivs, n_indivs), dtype=dtype)
    for chunk_i, i in enumerate(range(0, num_snps, n_indivs)):
        snps_array = sp.array(snps[i:i + n_indivs])
        snps_array = snps_array.T
        norm_snps_array = (snps_array - sp.mean(snps_array, 0)) / sp.std(snps_array, 0)
        x = sp.mat(norm_snps_array.T)
        k_mat += x.T * x
        sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, ((chunk_i + 1.0) * n_indivs) / num_snps))))
        sys.stdout.flush()
    k_mat = k_mat / float(num_snps)
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





def get_kinship(call_method_id=75, data_format='binary', method='ibs', n_removed_snps=None, remain_accessions=None,
                scaled=True, min_mac=5, sd=None, debug_filter=1, return_accessions=False):
    """
    Loads and processes the kinship matrix
    """
    import dataParsers as dp
    import env
    if method == 'ibd':
        if sd != None:
            k = sd.get_ibd_kinship_matrix()
            if scaled:
                k = scale_k(k)
            return k
        else:
            raise NotImplementedError('Currently only IBS kinship matrices are supported')
    elif method == 'ibs':
        if call_method_id:
            file_prefix = '%s%d/kinship_%s_%s' % (env.env['cm_dir'], call_method_id, method, data_format)
            kinship_file = file_prefix + '_mac%d.h5py' % min_mac
            if os.path.isfile(kinship_file):
                log.debug('Found kinship file: %s' % kinship_file)
                d = load_kinship_from_file(kinship_file, scaled=False)
                k = d['k']
                k_accessions = d['accessions']
                n_snps = d['n_snps']
            else:
                log.debug("Didn't find kinship file: %s, now generating one.." % kinship_file)
                try:
                    sd = dp.load_snps_call_method(call_method_id=call_method_id, data_format=data_format, min_mac=min_mac,
                                               debug_filter=debug_filter)
                except Exception:
                    if sd != None:
                        k = sd.get_ibs_kinship_matrix()
                    if scaled:
                        k = scale_k(k)
                    return k

                k = sd.get_ibs_kinship_matrix()
                k_accessions = sd.accessions
                n_snps = sd.num_snps()
                save_kinship_to_file(kinship_file, k, sd.accessions, n_snps)
            if n_removed_snps != None and remain_accessions != None:
                k = update_k_monomorphic(n_removed_snps, k, k_accessions, n_snps, remain_accessions,
                                                 kinship_type='ibs', dtype='single')
                if scaled:
                    k = scale_k(k)
                return k
            else:
                if scaled:
                    k = scale_k(k)
                if return_accessions:
                    return k, k_accessions
                else:
                    return k

    else:
        log.error('Method %s is not implemented' % method)
        raise NotImplementedError

