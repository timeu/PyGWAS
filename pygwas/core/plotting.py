import sys,os
import os.path
import numpy
import h5py
import math
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


SUPPORTED_FORMAT=('png','pdf')

def plot_gwas_result(hdf5_file,output_file,chrs=None,mac=15):
    format = os.path.splitext(output_file)[1][1:].strip().lower()
    if format not in SUPPORTED_FORMAT:
        raise Exception('%s not supported format'%format)
    if chrs is None:
        chrs = ['chr1','chr2','chr3','chr4','chr5']
    f = h5py.File(hdf5_file,'r')
    bh_thres,bonferroni_threshold,max_score,num_scores,min_score = _get_gwas_infos(f)
    chr_label = ''

    offset = 0 
    markersize=3
    color_map = {'chr1':'b', 'chr2':'g', 'chr3':'r', 'chr4':'c', 'chr5':'m'}
    plt.figure(figsize=(11, 3.8))
    plt.axes([0.045, 0.15, 0.99, 0.61])
    
    ticklist = []
    ticklabels = []
    for chr in chrs:
        chr_data = _get_chr_data(f,chr,mac)
        newPosList = [offset + pos for pos in chr_data['positions']]
        plt.plot(newPosList,chr_data['scores'],".", markersize=markersize, alpha=0.7, mew=0,color=color_map[chr])
        oldOffset = offset
        chr_end = chr_data['positions'][-1]
        offset =+newPosList[-1]
        for j in range(oldOffset, offset, 4000000):
            ticklist.append(j)
        for j in range(0, chr_end, 4000000):
            if j % 8000000 == 0 and j < chr_end - 4000000 :
                ticklabels.append(j / 1000000)
            else:
                ticklabels.append("")
    x_range = offset
    score_range = max_score - min_score
    padding = 0.05*(score_range)
    bonf_handle, = plt.plot([0, x_range], [bonferroni_threshold, bonferroni_threshold],color='r', linestyle="--",linewidth = 0.5)
    if bh_thres is not None:
        bh_handle, = plt.plot([0, x_range], [bh_thres, bh_thres], color='b', linestyle='--',linewidth = 0.5)
        plt.figlegend((bonf_handle,bh_handle),('Bonferroni','Benjamini Hochberg'),'upper right')
    else:
        plt.figlegend((bonf_handle,),('Bonferroni',),'upper right')
    plt.axis([-x_range * 0.01, x_range * 1.01, min_score - padding, max_score + padding])
    plt.xticks(ticklist, ticklabels, fontsize='x-small')
    plt.ylabel('$-log(p-$value$)$',size="large")
    plt.xlabel('Mb')
    if len(chrs) == 1:
         plt.title('Chromosome %s' % chr[3])

    f.close()
    if format == 'pdf':
        plt.savefig(output_file, format=format)
    else:
        plt.savefig(output_file, format=format,dpi=300,bbox_inches='tight')
    plt.clf()
    plt.close()
    return output_file




def _get_gwas_infos(f):
    attrs = f['pvalues'].attrs
    analysis_method = attrs.get('analysis_method',None)
    bh_thres = attrs.get('bh_thres',None)
    bonferroni_threshold = attrs.get('bonferroni_threshold',None)
    if bonferroni_threshold is None:
        bonferroni_threshold = attrs.get('bonferroniScore',None)
    if bonferroni_threshold is None:
        bonferroni_threshold = 1.0 / (num_snps * 20.0)
    max_score = attrs.get('max_score',None)
    if max_score is None:
    	max_score = attrs.get('maxScore',None)
    num_scores = 0
    min_score = 0
    for chr in ['chr1','chr2','chr3','chr4','chr5']:
        num_scores += len(f['pvalues'][chr]['positions'])
    return (bh_thres,bonferroni_threshold,max_score,num_scores,min_score)

def _get_chr_data(f,chr,mac=15):
    chr_data = f['pvalues'][chr]
    if 'macs' in chr_data:
        data = numpy.rec.fromarrays((chr_data['scores'][:],chr_data['positions'][:],chr_data['macs'][:]),dtype=[('scores', float), ('positions', int),('macs',int)])
    else:
        data = numpy.rec.fromarrays((chr_data['scores'][:],chr_data['positions'][:]),dtype=[('scores', float), ('positions', int)]) 
    data.sort(order='positions')
    if 'macs' in chr_data:
        ix = data['macs'] > mac
        return data[ix]
    return data

