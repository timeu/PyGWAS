import genotype
import numpy as np
import random
import logging
import operator

logger = logging.getLogger(__name__)

class permute():
    def __init__(self):
        pass

    def load_vectors(self, target, match):
        ' load input as 0,1 vectors with length equal to total number of snps '
        self.target = target
        self.match = match
        if self.target.shape != self.match.shape:
            raise ValueError('Non matching number of snps')

    def rotate_and_count(self, step):
        ' Do one set of rotation and count number of matches'
        if step == 0:
            return np.sum(self.target*self.match)
        else:
            testtarget = np.hstack((self.target[step:],self.target[:step]))
            return np.sum(testtarget*self.match)

    def permute_p_val(self, num, key='gene_set'):
        ' Do multiple rotations to obtain a p-value '
        originalmatch = self.rotate_and_count(0)
        num_reach_original = 0
        randomset = random.sample(xrange(1,self.target.shape[0]), num)
        num = float(num)
        for ix,step in enumerate(randomset):
            if ix % (num / 10) == 0:
                logger.info('Permutation %s/%s' % (ix,num),extra={'progress':(5 + ix/num*95),'key':key})
            if self.rotate_and_count(step)>=originalmatch:
                num_reach_original += 1
        return float(num_reach_original)/num

def enrichment(genes_lists,genotype_data,top_snps,windowSize,numberOfPermutations):
    chr_regions = genotype_data.chr_regions
    pos_chr_list= np.asarray(zip(genotype_data.chromosomes,genotype_data.positions),dtype='int32')
    logger.info('Retrieve TOP SNP Matrix')
    top_snps_matrix = _get_snps_snps_matrix(pos_chr_list,top_snps,chr_regions)
    pvals = {}
    for key,genes in genes_lists.items():
        logger.info('Retrieve Gene Matrix (windowsize:%s) for %s' % (windowSize,key))
        gene_matrix = _get_gene_snps_matrix(pos_chr_list,genes,chr_regions,windowSize)
        per = permute()
        per.load_vectors(gene_matrix,top_snps_matrix)
        logger.info('Starting permutation for %s (#: %s)' % (key, numberOfPermutations),extra={'progress':5,'key':key})
        pvals[key] = per.permute_p_val(numberOfPermutations,key=key)
        logger.info('Finished permutation for %s ' % key,extra={'progress':95,'key':key})
    return pvals

def _get_snps_snps_matrix(snps,top_snps,chr_regions):
    # sort by chr and position
    #use lexsort
    top_snps.sort(order=['chr','positions'])
    chr_start = 0
    chr_start_ix = 0
    indices = []
    vector = np.zeros((len(snps),),dtype='int32')
    for snp in top_snps:
        chr = int(snp[0])
        if chr != chr_start:
            chr_start = chr
            chr_start_ix = chr_regions[chr-1][0]
            chr_end_ix = chr_regions[chr-1][1]
        ix = chr_start_ix + snps[chr_start_ix:chr_end_ix][:,1].searchsorted(snp[1])
        indices.append(ix)
    vector[indices] = 1
    return vector

def _get_gene_snps_matrix(snps,genes,chr_regions,windowsize):
    # sort by gene name and first position
    sorted_genes = sorted(genes, key=operator.itemgetter(3, 0))
    chr_start = 0
    chr_start_ix = 0
    indices = []
    vector = np.zeros((len(snps),),dtype='int32')
    for gene in sorted_genes:
        chr = int(gene[3][2])
        if chr != chr_start:
            chr_start = chr
            chr_start_ix = chr_regions[chr-1][0]
            chr_end_ix = chr_regions[chr-1][1]
        ix_start = chr_start_ix+ snps[chr_start_ix:chr_end_ix][:,1].searchsorted(gene[0]-windowsize,side='left')
        ix_end = chr_start_ix + snps[chr_start_ix:chr_end_ix][:,1].searchsorted(gene[1]+windowsize,side='right')
        vector[ix_start:ix_end] = 1
    return vector