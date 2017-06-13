import logging
import mtcorr
import statistics as stats
import math
import h5py
import numpy
import sys


log = logging.getLogger(__name__)


def load_from_hdf5(filename):
    f = h5py.File(filename,'r')
    quantiles_dict = {}
    stats =  {}
    if 'quantiles' in f:
        quantiles_dict['exp_quantiles'] = f['quantiles']['quantiles'][:,0].tolist()
        quantiles_dict['quantiles'] = f['quantiles']['quantiles'][:,1].tolist()
        quantiles_dict['exp_log_quantiles'] = f['quantiles']['log_quantiles'][:,0].tolist()
        quantiles_dict['log_quantiles'] = f['quantiles']['log_quantiles'][:,1].tolist()
        stats['quantiles_dict'] = quantiles_dict
    pvals_group = f['pvalues']
    method = pvals_group.attrs.get('analysis_method','')
    transformation = pvals_group.get('transformation','')
    if 'ks_stat' in pvals_group.attrs:
        stats['ks_stats'] = {'D':pvals_group.attrs['ks_stat']}
    if 'ks_pval' in pvals_group.attrs:
        stats['ks_stats']['p_val'] = pvals_group.attrs['ks_pval']
    if 'med_pval' in pvals_group.attrs:
        stats['med_pval'] =  pvals_group.attrs['med_pval']
    if 'bh_thres' in pvals_group.attrs:
        stats['bh_thres_d'] = {'thes_pval': math.pow(10,-pvals_group.attrs['bh_thres'])}
    chromosomes = []
    positions = []
    scores = []
    mafs = []
    macs = []
    additional_columns = {}
    chrs = map(lambda x:x[3:],f['pvalues'].keys())

    for ix,chr in enumerate(chrs):
        chr_group = pvals_group['chr%s'% chr]
        chromosomes.extend([chr]*len(chr_group['positions']))
        positions.extend(chr_group['positions'][:].tolist())
        scores.extend(chr_group['scores'][:].tolist())
        mafs.extend(chr_group['mafs'][:].tolist())
        macs.extend(chr_group['macs'][:].tolist())
        for i,key in enumerate(chr_group.keys()):
            if key not in ('positions','scores','mafs','macs'):
                values = chr_group[key][:].tolist()
                if key not in additional_columns:
                    additional_columns[key] = values
                else:
                    additional_columns[key].extend(values)
    f.close()
    scores = map(lambda x:math.pow(10,-1*x), scores)
    maf_dict = {'mafs':mafs,'macs':macs}
    return GWASResult(chrs,chromosomes,positions,scores,maf_dict,method,transformation,stats=stats,additional_columns=additional_columns)


def load_from_csv(filename):
    chromosomes =  []
    positions =  []
    pvals =  []
    mafs = []
    macs = []
    additional_columns = {}
    chrs = []
    chr = None
    is_pval = False
    with open(filename,'r') as f:
        header = f.readline().rstrip()
        add_header = header.split(",")[5:]
        for key in add_header:
            key = key.replace('"','')
            additional_columns[key] = []
        for row in f:
            fields = row.rstrip().split(",")
            if chr != fields[0]:
                chr = fields[0]
                chrs.append(chr)
            chromosomes.append(chr)
            positions.append(int(float(fields[1])))
            pvals.append(float(fields[2]))
            mafs.append(float(fields[3]))
            macs.append(int(float(fields[4])))
            if len(add_header) > 0:
                for i,key in enumerate(add_header):
                    key = key.replace('"','')
                    addit_value = None
                    if fields[(5+i)] != '':
                        addit_value = float(fields[(5+i)])
                    additional_columns[key].append(addit_value)
    is_pval = max(pvals) <= 1.0
    if is_pval is False:
        pvals = map(lambda x:math.pow(10,-1*x),pvals)
    return GWASResult(chrs,chromosomes,positions,pvals,{'mafs':mafs,'macs':macs},additional_columns = additional_columns)




class GWASResult(object):


    def __init__(self,chrs,chromosomes,positions,pvals,maf_dict,method = 'N/A',transformation = None,stats = None,additional_columns = None,step_stats = None):
        self.ix_with_bad_pvalues = ix_with_bad_pvalues = numpy.where(pvals == 0.0)[0]
        if len(ix_with_bad_pvalues) > 0:
            pvals[ix_with_bad_pvalues] = sys.float_info.min
        self.pvals = pvals
        self.method = method
        self.transformation = transformation
        self.chrs = chrs
        self.chromosomes = chromosomes
        self.positions = positions
        self.stats = stats
        self.maf_dict = maf_dict
        self.additional_columns = additional_columns
        self.step_stats = step_stats
        self.bonferroni_threshold = -math.log10(0.05 / len(pvals))
        self.min_pval = min(pvals)
        if not self.stats:
            self._calculate_stats_()



    def _calculate_stats_(self):
        log.info('Calculating Benjamini-Hochberg threshold',extra={'progress':90})
        #Calculate Benjamini-Hochberg threshold
        self.stats = {}
        self.stats['bh_thres_d'] = mtcorr.get_bhy_thres(self.pvals, fdr_thres=0.05)
        #Calculate Median p-value
        self.stats['med_pval'] = stats.calc_median(self.pvals)
        #Calculate the Kolmogorov-Smirnov statistic
        self.stats['ks_stats'] = stats.calc_ks_stats(self.pvals)
        self.stats['quantiles_dict'] = stats.calculate_qqplot_data(self.pvals)


    def get_top_snps(self,top_ratio=2500):
        data = numpy.core.records.fromrecords(zip(self.chromosomes, self.positions, self.pvals, self.maf_dict['mafs'], self.maf_dict['macs'],*self.additional_columns.values()),names='chr,positions,scores,mafs,macs')
        data_to_return=[]
        for ix,chr in enumerate(self.chrs):
            chr_data = data[numpy.where(data['chr'] == chr)]
            chr_data =chr_data[chr_data['scores'].argsort()[::]][:top_ratio]
            data_to_return.append(chr_data)
        return numpy.concatenate(data_to_return)


    def save_as_csv(self,csv_file):
        data = numpy.array(zip(self.chromosomes, self.positions, self.pvals, self.maf_dict['mafs'], self.maf_dict['macs'],*self.additional_columns.values()))
        data =data[numpy.lexsort((data[:,1],data[:,0]))]
        additional_column_headers = self.additional_columns.keys()
        header = ['chromosomes','positions','pvals','mafs','macs']
        header.extend(additional_column_headers)
        with open(csv_file,'w') as f:
            f.write(','.join(header)+"\n")
            for row in data:
               rows_to_write = row.tolist()
               rows_to_write[0] = int(rows_to_write[0])
               rows_to_write[1] = int(rows_to_write[1])
               rows_to_write[4] = int(rows_to_write[4])
               f.write(','.join(map(str,rows_to_write))+"\n")



    def save_as_hdf5(self,hdf5_file):
        positions = self.positions
        chromosomes = self.chromosomes
        maf_dict = self.maf_dict
        scores = map(lambda x:-math.log10(x), self.pvals)
        quantiles_dict = self.stats['quantiles_dict']
        f = h5py.File(hdf5_file,'w')

        # store quantiles
        quant_group = f.create_group('quantiles')
        quantiles_array = zip(quantiles_dict['exp_quantiles'],quantiles_dict['quantiles'])
        log_quantiles_array = zip(quantiles_dict['exp_log_quantiles'],quantiles_dict['log_quantiles'])
        quant_group.create_dataset('quantiles',(len(quantiles_dict['quantiles']), 2),'f8',data=quantiles_array)
        quant_group.create_dataset('log_quantiles',(len(quantiles_dict['log_quantiles']), 2),'f8',data=log_quantiles_array)

        #store pvalues
        pvals_group = f.create_group('pvalues')
        if len(self.ix_with_bad_pvalues) > 0:
            pvals_group.attrs['ix_with_bad_pvalues'] = self.ix_with_bad_pvalues
        pvals_group.attrs['numberOfSNPs'] = len(scores)
        pvals_group.attrs['max_score'] = max(scores)
        if self.method is not None:
            pvals_group.attrs['analysis_method'] = self.method
        transformation = "raw"
        if self.transformation is not None:
            transformation = self.transformation
        pvals_group.attrs['transformation'] = transformation
        pvals_group.attrs['bonferroni_threshold'] = self.bonferroni_threshold
        pvals_group.attrs['ks_stat'] = self.stats['ks_stats']['D']
        pvals_group.attrs['ks_pval'] = self.stats['ks_stats']['p_val']
        pvals_group.attrs['med_pval'] = self.stats['med_pval']
        pvals_group.attrs['bh_thres'] =-math.log10(self.stats['bh_thres_d']['thes_pval'])

        data = numpy.core.records.fromrecords(zip(chromosomes, positions, scores, maf_dict['mafs'], maf_dict['macs'],*self.additional_columns.values()),names='chr,positions,scores,mafs,macs')
        for ix,chr in enumerate(self.chrs):
            chr_group = pvals_group.create_group("chr%s" % chr)
            chr_data = data[numpy.where(data['chr'] == chr)]
            chr_data =chr_data[chr_data['scores'].argsort()[::-1]]
            positions = chr_data['positions']
            chr_group.create_dataset('positions',(len(positions),),'i4',data=positions)
            scores = chr_data['scores']
            chr_group.create_dataset('scores',(len(scores),),'f8',data=scores)
            mafs = chr_data['mafs']
            chr_group.create_dataset('mafs',(len(mafs),),'f8',data=mafs)
            macs = chr_data['macs']
            chr_group.create_dataset('macs',(len(macs),),'i4',data=macs)

            if len(chr_data.dtype) > 5:
                for i,key in enumerate(self.additional_columns.keys()):
                    values = chr_data['f%s'% (5+i)]
                    chr_group.create_dataset(key,values.shape,values.dtype,data=values)
        f.close()
