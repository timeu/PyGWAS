import logging
import mtcorr
import statistics as stats
import math
import h5py
import numpy

log = logging.getLogger(__name__)

class GWASResult(object):


    def __init__(self,pvals,stats,method,transformation,genotype,additional_columns = None):
        self.pvals = pvals
        self.stats = stats
        self.method = method
        self.transformation = transformation
        self.genotype = genotype
        self.additional_columns = additional_columns
        self._calculate_stats_()


    def _calculate_stats_(self):
        log.info('Calculating Benjamini-Hochberg threshold',extra={'progress':90})
        #Calculate Benjamini-Hochberg threshold
        self.bh_thres_d = mtcorr.get_bhy_thres(self.pvals, fdr_thres=0.05)
        #Calculate Median p-value
        self.med_pval = stats.calc_median(self.pvals)
        #Calculate the Kolmogorov-Smirnov statistic
        self.ks_stats = stats.calc_ks_stats(self.pvals)
        self.quantiles_dict = stats.calculate_qqplot_data(self.pvals)

    
      
    
    
    def save_as_hdf5(self,hdf5_file):
        positions = self.genotype.positions
        chromosomes = self.genotype.chromosomes
        maf_dict = self.genotype.get_mafs()
        scores = map(lambda x:-math.log10(x), self.pvals)

        f = h5py.File(hdf5_file,'w') 

        # store quantiles
        quant_group = f.create_group('quantiles')
        quantiles_array = zip(self.quantiles_dict['exp_quantiles'],self.quantiles_dict['quantiles'])
        log_quantiles_array = zip(self.quantiles_dict['exp_log_quantiles'],self.quantiles_dict['log_quantiles'])
        quant_group.create_dataset('quantiles',(len(self.quantiles_dict['quantiles']), 2),'f8',data=quantiles_array)
        quant_group.create_dataset('log_quantiles',(len(self.quantiles_dict['log_quantiles']), 2),'f8',data=log_quantiles_array)
        
        #store pvalues
        pvals_group = f.create_group('pvalues')
        pvals_group.attrs['numberOfSNPs'] = len(scores)
        pvals_group.attrs['max_score'] = max(scores)
        pvals_group.attrs['analysis_method'] = self.method
        transformation = "raw"
        if self.transformation is not None:
            transformation = self.transformation
        pvals_group.attrs['transformation'] = transformation
        pvals_group.attrs['bonferroni_threshold'] = -math.log10(0.05 / len(scores))
        pvals_group.attrs['ks_stat'] = self.ks_stats['D']
        pvals_group.attrs['ks_pval'] = self.ks_stats['p_val']
        pvals_group.attrs['med_pval'] = self.med_pval
        pvals_group.attrs['bh_thres'] =-math.log10(self.bh_thres_d['thes_pval'])

        data = numpy.array(zip(chromosomes, positions, scores, maf_dict['mafs'], maf_dict['macs'],*self.additional_columns.values()))
        
        for chr in range(1,6):
            chr_group = pvals_group.create_group('chr%s' % chr)
            chr_data = data[numpy.where(data[:,0] == chr)]
            chr_data =chr_data[chr_data[:,2].argsort()[::-1]]
            positions = chr_data[:,1]
            chr_group.create_dataset('positions',(len(positions),),'i4',data=positions)
            scores = chr_data[:,2]
            chr_group.create_dataset('scores',(len(scores),),'f8',data=scores)
            mafs = chr_data[:,3]
            chr_group.create_dataset('mafs',(len(mafs),),'f8',data=mafs)
            macs = chr_data[:,4]
            chr_group.create_dataset('macs',(len(macs),),'i4',data=macs)

            if chr_data.shape[1] > 5: 
                for i,key in enumerate(self.additional_columns.keys()):
                    values = chr_data[:,5+i]
                    chr_group.create_dataset(key,values.shape,values.dtype,data=values) 
        f.close()
