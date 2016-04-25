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

def plot_gwas_result(gwas_result,output_file,chrs=None,mac=15):
    chr_map = None
    if chrs is not None:
        chr_map = set(chrs)
    chrs = gwas_result.chrs
    if chr_map is not None:
        for ix,chr in enumerate(chrs):
            if chr not in chr_map:
                chrs[ix] = None

    if len(filter(lambda x: x is not None,chrs)) == 0:
        raise ValueError('Chromosomes %s  not found' % chr_map)
    format = os.path.splitext(output_file)[1][1:].strip().lower()
    if format not in SUPPORTED_FORMAT:
        raise Exception('%s not supported format'%format)
    bh_thres,bonferroni_threshold,max_score,num_scores,min_score = _get_gwas_infos(gwas_result)
    chr_label = ''
    data = _get_data(gwas_result)
    offset = 0 
    markersize=3
    color_map = ['b', 'g', 'r', 'c', 'm']
    plt.figure(figsize=(11, 3.8))
    plt.axes([0.045, 0.15, 0.99, 0.61])
    
    ticklist = []
    ticklabels = []

    for ix,chr in enumerate(chrs):
        if chr is None:
            continue
        chr_data = _get_chr_data(data,chr,mac)
        newPosList = [offset + pos for pos in chr_data['positions']]
        color_ix = ix
        if color_ix >= len(color_map):
            color_ix = (ix+1) % color_ix
        color = color_map[color_ix]
        plt.plot(newPosList,chr_data['scores'],".", markersize=markersize, alpha=0.7, mew=0,color=color)
        oldOffset = offset
        chr_end = chr_data['positions'][-1] if len(chr_data['positions']) > 0 else 0
        offset =+ newPosList[-1] if len(newPosList) > 0 else 0
        for j in range(oldOffset, offset, 4000000):
            ticklist.append(j)
        for j in range(0, chr_end, 4000000):
            if j % 8000000 == 0 and j < chr_end - 4000000 :
                ticklabels.append(j / 1000000)
            else:
                ticklabels.append("")
    x_range = offset
    max_score = max([bonferroni_threshold,bh_thres,max_score])

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
    if format == 'pdf':
        plt.savefig(output_file, format=format)
    else:
        plt.savefig(output_file, format=format,dpi=300,bbox_inches='tight')
    plt.clf()
    plt.close()
    return output_file




def _get_gwas_infos(gwas_result):
    analysis_method = gwas_result.method
    bh_thres = gwas_result.stats['bh_thres_d'].get('thes_pval',None)
    if bh_thres is not None:
        bh_thres = -math.log10(bh_thres)
    bonferroni_threshold = gwas_result.bonferroni_threshold
    max_score = -math.log10(gwas_result.min_pval)
    num_scores = len(gwas_result.positions)
    min_score = 0
    return (bh_thres,bonferroni_threshold,max_score,num_scores,min_score)

def _get_data(gwas_result):
    scores = map(lambda x: -math.log10(x),gwas_result.pvals)
    if gwas_result.maf_dict is not None and 'macs' in gwas_result.maf_dict:
        data = numpy.rec.fromarrays((gwas_result.chromosomes,gwas_result.positions,scores,gwas_result.maf_dict['macs']),dtype=[('chromosomes','S2'),('positions', int),('scores', float),('macs',int)])
    else:
        data = numpy.rec.fromarrays((gwas_result.chromosomes,gwas_result.positions,scores),dtype=[('chromosomes','S2'),('positions', int),('scores', float)])
    return data

def _get_chr_data(data,chr,mac=15):
    chr_data = data[data['chromosomes'] == chr]
    chr_data.sort(order='positions')
    if len(chr_data.dtype) == 4:
        ix = chr_data['macs'] > mac
        return chr_data[ix]
    return chr_data

