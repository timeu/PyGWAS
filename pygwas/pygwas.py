#!/usr/bin/env python
# coding: utf-8
"""
    pygwas
    ~~~~~~~~~~~~~

    The main module for running Genome Wide Association studies

    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""

from __init__ import __version__,__updated__,__date__
import argparse
from core import kinship
from core import gwas
import logging, logging.config
from core import mtcorr
from core import statistics as stats
from core import phenotype
from core import genotype
from core.result import GWASResult
from core import enrichment
from core import result
from core import plotting
import numpy as np
import re
import os
import sys
import csv
import h5py
import json
from operator import itemgetter

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'default': {
            'format': '%(asctime)s %(levelname)s %(name)s %(message)s'
        },
    },
    'handlers': {
        'stdout':{
            'class' : 'logging.StreamHandler',
            'stream'  : 'ext://sys.stdout',
            'formatter': 'default',
        },
        'stderr':{
            'class' : 'logging.StreamHandler',
            'stream'  : 'ext://sys.stderr',
            'level':'ERROR',
            'formatter': 'default',
        },
    },
    'root': {
        'handlers': ['stdout','stderr'],
        'level': 'INFO',
    },
}

CHR_POS_PATTERN=re.compile(r"^(\w):(\d+)$")

HDF5_FILE_EXT = ('.hdf5','.h5','.hdf')

SUPPORTED_FILE_EXT =  HDF5_FILE_EXT  + ('.csv',)

logging.config.dictConfig(LOGGING)
log = logging.getLogger()

def get_parser(program_license,program_version_message):

    snp_help ="""Position/region for LD calculation.\n
    There are two ways to specify the region/positions:\n
      1.) GWAS result file (*.csv or *.hdf5). The top fraction (see -r/--range paramter) of SNPs will be taken and LD will be calculated for those positions\n
      2.) chr:position (i.e. 2:25212). Specify midpoint and take +/- range (see -r/--range paramter) to calculate LD\n
    """

    parser = argparse.ArgumentParser(description=program_license)
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument('-l', '--loglevel',dest='log_level', help='Set log level', default='INFO', choices=['DEBUG','INFO','WARNING','ERROR'])
    subparsers = parser.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')
    analysis_parser = subparsers.add_parser('run',help='Run a GWAS analysis')

    analysis_parser.add_argument("-t", "--transformation", dest="transformation", help="Apply a transformation to the data. Default[None]", default="none", choices=phenotype.SUPPORTED_TRANSFORMATIONS)
    analysis_parser.add_argument("-a", "--analysis_method", dest="analysis_method", help="analyis method to use",required=True,choices=["lm", "emma", "emmax", "kw", "ft", "emmax_anova", "lm_anova", "emmax_step", "lm_step","loc_glob_mm","amm"])
    analysis_parser.add_argument("-g", "--genotype", dest="genotype_folder", help="folder with the genotypes for the GWAS analysis", required=True,metavar="FOLDER")
    analysis_parser.add_argument("-k", "--kinship", dest="kinship", help="Specify the file containing the kinship matrix. (otherwise default file is used or it's generated.)", metavar="FILE" )
    analysis_parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file")
    analysis_parser.add_argument("-l","--calc-ld",default=False,dest="calc_ld",action="store_true",help="If set, will also calculte and store LD for the top SNPs (2500 per chromosome). If the output file is an hdf5 file it will add it to the hdf5 file. For a CSV file it will store it in a new file")
    analysis_parser.add_argument(dest="file", help="csv file containing phenotype values", metavar="FILE")
    analysis_parser.set_defaults(func=run)

    convert_parser = subparsers.add_parser('convert',help='Convert from HDF5 to CSV and the other way round')
    convert_parser.add_argument("-i", "--input_file", dest="inputfile", help="Name of the input file",required=True)
    convert_parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file",required=True)
    convert_parser.set_defaults(func=convert)

    plotter_parser = subparsers.add_parser('plot',help='Plot a GWAS result')
    plotter_parser.add_argument('-c','--chr',dest='chr',
	 help='Chromosome to plot. If not specified prints all chromosomes (Default:None)',default=None)
    plotter_parser.add_argument('-m','--macs',dest='macs',default=15,type=int,help='Minor Allele Count filter (Default: 15)')
    plotter_parser.add_argument('-s','--size',dest='marker_size',default=10,type=int,help='Size of the markers in the Manhattan plot (Default: 10)')
    plotter_parser.add_argument("-o",'--output',dest='output',required=True,help='The output image file')
    plotter_parser.add_argument(dest="file", help="GWAS result file (.hdf5 or .csv)", metavar="FILE")
    plotter_parser.set_defaults(func=plot)


    qq_plotter_parser = subparsers.add_parser('qqplot',help='Plot a QQ-plots for a GWAS result')
    qq_plotter_parser.add_argument("-o",'--output',dest='output',required=True,help='The output image file')
    qq_plotter_parser.add_argument(dest="file", help="GWAS result file (.hdf5 or .csv)", metavar="FILE")
    qq_plotter_parser.set_defaults(func=qq_plot)

    stats_parser = subparsers.add_parser('stats',help='Retrieve some stats')
    stats_parser.add_argument("-s", "--type", dest="type",required=True, help="type of the statistics to return",choices=["all","pseudo","shapiro"])
    stats_parser.add_argument("-t", "--transformation", dest="transformation",default='none', help="transformation to apply",choices=phenotype.SUPPORTED_TRANSFORMATIONS)
    stats_parser.add_argument("-g", "--genotype", dest="genotype_folder", help="folder with the genotypes for the GWAS analysis", metavar="FOLDER")
    stats_parser.add_argument("-k", "--kinship", dest="kinship", help="Specify the file containing the kinship matrix. (otherwise default file is used or it's generated.)", metavar="FILE" )
    stats_parser.add_argument(dest="file", help="csv file containing phenotype values",  metavar="FILE")
    stats_parser.set_defaults(func=calculate_stats)

    ld_parser = subparsers.add_parser('ld',help='Calculate and plot LD')
    ld_parser.add_argument(dest="genotype_file", help="genotype file", metavar="GENOTYPE FILE")
    ld_parser.add_argument(dest="output_file",help="output file (.hdf5 or .csv)",metavar="OUTPUT FILE")
    ld_parser.add_argument("-p","--positions",dest="positions",help=snp_help,required="True",metavar="position/FILE")
    ld_parser.add_argument("-r","--range",dest="range",help="Range of the LD block (Default: 2500 per chromosome)",default="2500",type=int)
    ld_parser.add_argument("-a", "--accessions", dest="acession_file", help="optional csv file with list of accession ids to filter", metavar="FILE")
    ld_parser.set_defaults(func=calc_ld)

    kinship_parser = subparsers.add_parser('kinship',help='Calculate Kinship')
    kinship_parser.add_argument(dest="genotype_file", help="genotype file", metavar="GENOTYPE FILE")
    kinship_parser.add_argument(dest="output_file",help="output file (.hdf5)",metavar="OUTPUT FILE")
    kinship_parser.add_argument("-t", "--type", dest="type",required=True, help="type of the kinship",choices=["ibs","ibd"])
    kinship_parser.set_defaults(func=calc_kinship)

    enrichment_parser = subparsers.add_parser('enrichment',help='Enrichment Analysis')
    enrichment_parser.add_argument(dest="genotype_folder", help="genotype folder", metavar="GENOTYPE-FOLDER")
    enrichment_parser.add_argument(dest="gwas_file",help="file with pvalues (.hdf5)",metavar="PVALUES-FILE")
    enrichment_parser.add_argument("-g", "--gene_files", dest="genes_file",help="files with genes + locations (multiple can be passed)",required=True,nargs='+')
    enrichment_parser.add_argument("-w", "--window_size", dest="window_size",default=20000, type=int, help="Window size around genes (default: 20kb)",)
    enrichment_parser.add_argument("-p", "--permutation_count", dest="permutation_count",default=10000, type=int, help="Number of permutations (default: 10.000)",)
    enrichment_parser.add_argument("-t", "--top_snps_count", dest="top_snps_count",default=1000, type=int, help="Number of top SNPs to check enrichment for (default: 1000)",)
    enrichment_parser.set_defaults(func=calc_enrichment)

    transformation_parser = subparsers.add_parser('transform',help='Transform phenotype')
    transformation_parser.add_argument(dest="file", help="csv file containing phenotype values",  metavar="FILE")
    transformation_parser.add_argument("-t", "--transformation", default='most_normal', dest="transformation",help="transformation to use (default: use the most normal one", choices=phenotype.SUPPORTED_TRANSFORMATIONS)
    transformation_parser.add_argument("-o",'--output',dest='output',required=True,help='The output of the transformed phenotype file')
    transformation_parser.set_defaults(func=transform_phenotype)



    return parser


def main():
    '''Command line options.'''
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = "The main module for running Genome Wide Association studies"
    program_license = '''%s

  Created by Ümit Seren on %s.
  Copyright 2012 Gregor Mendel Institute. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))
    # Process arguments
    parser = get_parser(program_license,program_version_message)
    args = vars(parser.parse_args())
    try:
        log_level = args['log_level']
        log.setLevel(log_level)
        args['func'](args)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        log.exception(e)
        return 2

def _get_indicies_(phen_acc,geno_acc):
    sd_indices_to_keep = set()
    pd_indices_to_keep = []
    for i, acc in enumerate(geno_acc):
        for j, et in enumerate(phen_acc):
            if str(et) == str(acc):
                sd_indices_to_keep.add(i)
                pd_indices_to_keep.append(j)
    sd_indices_to_keep = list(sd_indices_to_keep)
    sd_indices_to_keep.sort()
    return (sd_indices_to_keep,pd_indices_to_keep)


def calculate_stats(args):
    statistics = {}
    phenotype_file = args['file']
    genotype_folder = args['genotype_folder']
    stat_type =args['type']
    transformation = args['transformation']
    if 'phen_data' in args:
        phenData = args['phen_data']
    else:
        phenData = phenotype.parse_phenotype_file(phenotype_file)
    phenData.convert_to_averages()
    statistics['transformation'] = phenData.transform(transformation)
    if stat_type not in ('shapiro','all','pseudo'):
        raise Exception('%s not supported' % stat_type)
    elif stat_type in ('all','pseudo') and genotype_folder is None:
        raise Exception('%s for the stat_type shapiro and all, you need to provide specifiy the genotype_folder')

    if genotype_folder:
        genotypeData = _load_genotype_(genotype_folder)
        sd_indices_to_keep,pd_indices_to_keep = _get_indicies_(phenData.ecotypes,genotypeData.accessions)
        phenData.filter_ecotypes(pd_indices_to_keep)
        accessions = genotypeData.accessions[sd_indices_to_keep]
        if stat_type == 'all' or stat_type == 'pseudo':
            K = None
            kinship_file = args.get('kinship',None)
            if kinship_file is None:
                kinship_file = _get_kinship_file_(genotype_folder)
            K = kinship.load_kinship_from_file(kinship_file, accessions.tolist())['k']
            statistics['pseudo_heritability'] = gwas.calculate_pseudo_heritability(phenData.values,K)
    if stat_type == 'all' or stat_type == 'shapiro':
        statistics['shapiro'] = stats.calculate_sp_pval(phenData.values)

    print statistics
    return statistics


def convert(args):

    _,input_ext = os.path.splitext(args['inputfile'])
    _,output_ext = os.path.splitext(args['outputfile'])
    if input_ext == output_ext:
        raise Exception('use different file extension for input (%s) and output file (%s)' % (input_ext,output_ext))
    if input_ext not in SUPPORTED_FILE_EXT:
        raise Exception('The input file must have one of the supported extensions: %s' % str(SUPPORTED_FILE_EXT))
    if output_ext not in SUPPORTED_FILE_EXT:
        raise Exception('The output file must have one of the supported extensions: (%s)' % ', '.join(SUPPORTED_FILE_EXT))
    gwas_result = None
    if input_ext == '.csv':
        gwas_result = result.load_from_csv(args['inputfile'])
        gwas_result.save_as_hdf5(args['outputfile'])
    else:
        gwas_result = result.load_from_hdf5(args['inputfile'])
        gwas_result.save_as_csv(args['outputfile'])


def plot(args):
    _,ext = os.path.splitext(args['file'])
    if ext not in SUPPORTED_FILE_EXT:
        raise Exception('The input file must have one of the supported extensions: %s' % str(SUPPORTED_FILE_EXT))
    chrs = None
    marker_size = args['marker_size']
    if 'chr' in args and args['chr'] is not None and args['chr'] != '':
        chrs = [args['chr']]
    gwas_result = None
    if ext in HDF5_FILE_EXT:
        gwas_result = result.load_from_hdf5(args['file'])
    else:
        gwas_result = result.load_from_csv(args['file'])
    plotting.plot_gwas_result(gwas_result,args['output'],chrs,args['macs'],marker_size=marker_size)


def qq_plot(args):
    _,ext = os.path.splitext(args['file'])
    if ext not in SUPPORTED_FILE_EXT:
        raise Exception('The input file must have one of the supported extensions: %s' % str(SUPPORTED_FILE_EXT))
    gwas_result = None
    if ext in HDF5_FILE_EXT:
        gwas_result = result.load_from_hdf5(args['file'])
    else:
        gwas_result = result.load_from_csv(args['file'])
    plotting.plot_qq(gwas_result,args['output'])


def run(args):
    _,ext = os.path.splitext(args['outputfile'])
    if ext not in SUPPORTED_FILE_EXT:
        raise Exception('The output file must have one of the supported extensions: %s' % str(SUPPORTED_FILE_EXT))
    gwas_result = perform_gwas(args['file'],args['analysis_method'], args['genotype_folder'],args['transformation'],args['kinship'])
    if ext in HDF5_FILE_EXT:
        gwas_result.save_as_hdf5(args['outputfile'])
    else:
        gwas_result.save_as_csv(args['outputfile'])
    if args['calc_ld']:
        genotype_folder = args['genotype_folder']
        calc_ld_args = {'acession_file':args['file'],'positions':args['outputfile'],'genotype_file':os.path.join(genotype_folder,'all_chromosomes_binary.hdf5'),'range':2500}
        if ext in HDF5_FILE_EXT:
            calc_ld_args['output_file'] = args['outputfile']
        else:
            calc_ld_args['output_file'] = 'LD_%s' % args['outputfile']
        calc_ld(calc_ld_args)


def calc_ld(args):
    positions = args['positions']
    range = int(args['range'])
    is_chr_pos = CHR_POS_PATTERN.match(positions)
    genotypeData = genotype.load_hdf5_genotype_data(args['genotype_file'])
    accession_file = args.get('acession_file',None)
    if accession_file is not None:
        accessions = _load_accessions(accession_file)
        genotypeData.filter_accessions(accessions)
        genotypeData.filter_non_binary()
    if is_chr_pos:
        chr = is_chr_pos.group(1)
        position = int(is_chr_pos.group(2))
        log.info('Calculating LD for chr %s and position %s with +/- range of %s' % (chr,position, range))
        # get the positions
        abs_ix,ix,found = genotypeData.get_pos_ix(chr,position)
        min_ix = max(0,abs_ix - range)
        max_ix = min(abs_ix + range, genotypeData.genome_length)
        chr_pos_list = zip(genotypeData.chromosomes[min_ix:max_ix],genotypeData.positions[min_ix:max_ix])
    else:
        if os.path.exists(positions) == False:
            raise ValueError('%s path does not exist' % positions)
        # get the SNPs
        _,ext = os.path.splitext(positions)
        if ext not in SUPPORTED_FILE_EXT:
            raise Exception('The input file must have one of the supported extensions: %s' % str(SUPPORTED_FILE_EXT))
        if ext  in HDF5_FILE_EXT:
            gwas_result = result.load_from_hdf5(positions)
        else:
            gwas_result =  result.load_from_csv(positions)
        gwas_data = gwas_result.get_top_snps(range)
        chr_pos_list = zip(map(str,gwas_data['chr']),gwas_data['positions'])
    chr_pos_list = sorted(chr_pos_list,key=itemgetter(0,1))
    ld_data = genotypeData.calculate_ld(chr_pos_list)
    _save_ld_data(args['output_file'],ld_data,chr_pos_list)


def calc_kinship(args):
    genotypeData = genotype.load_hdf5_genotype_data(args['genotype_file'])
    type = args['type']
    if type == 'ibs':
        K = genotypeData.get_ibs_kinship_matrix(chunk_size = 10000)
    elif type == 'ibd':
        K = genotypeData.get_ibd_kinship_matrix(chunk_size = 10000)
    else:
        raise Exception('%s kinship type not supported' % type)
    kinship.save_kinship_to_file(args['output_file'],K,genotypeData.accessions,genotypeData.num_snps)


def calc_enrichment(args):
    genotype_folder = args['genotype_folder']
    genes_files = args['genes_file']
    gwas_file = args['gwas_file']
    log.info('Retrieving genes')
    genes = {}
    for genes_file in genes_files:
        with open(genes_file,'r') as f:
            genes[os.path.basename(genes_file)] = json.load(f)

    gwas_result = None
    _,input_ext = os.path.splitext(gwas_file)
    log.info('Opening GWAS file')
    if input_ext == '.csv':
        gwas_result = result.load_from_csv(gwas_file)
    else:
        gwas_result = result.load_from_hdf5(gwas_file)
    log.info('Loading genotype file')
    genotype_data = _load_genotype_(genotype_folder)
    top_snps = gwas_result.get_top_snps(args['top_snps_count'])
    top_snps.sort(order=['scores'])
    top_snps = top_snps[:args['top_snps_count']]
    pval = enrichment.enrichment(genes,genotype_data,top_snps,args['window_size'],args['permutation_count'])
    print json.dumps(pval)
    return pval



def perform_gwas(phenotype_file,analysis_method,genotype_folder,transformation=None,kinshipFile=None):
    phenData = phenotype.parse_phenotype_file(phenotype_file)  #load phenotype file
    additional_columns = {}
    genotypeData = _load_genotype_(genotype_folder)
    K = None
    n_filtered_snps = _prepare_data_(genotypeData,phenData,transformation)
    phen_vals = phenData.values
    if analysis_method in ['emma', 'emmax', 'emmax_anova', 'emmax_step', 'loc_glob_mm','amm']:
        #Load genotype file (in binary format)
        log.debug("Retrieving the Kinship matrix K.\n")
        if kinshipFile is None:   #Kinship file was supplied..
	        kinshipFile = _get_kinship_file_(genotype_folder)
        log.info('Loading kinship file: %s' % kinshipFile,extra={'progress':5})
        K = kinship.load_kinship_from_file(kinshipFile, genotypeData.accessions.tolist(),n_removed_snps=n_filtered_snps)['k']
        log.info('Done!')


    if analysis_method in ['kw']:
        res = gwas.kruskal_wallis(genotypeData, phen_vals)
    elif analysis_method in ['loc_glob_mm']:
        raise NotImplementedError
    elif analysis_method in ['emma']:
        res = gwas.emma(snps, phen_vals, K)
    elif analysis_method in ['emmax','amm']:
        d = gwas.emmax_step(phen_vals, genotypeData, K, [], emma_num=200)
        res = d['res']
    elif analysis_method in ['lm']:
        d = gwas.lin_reg_step(phen_vals, genotypeData, [])
        res = d['res']
    else:
        raise Exception('analysis method %s not supported' % analysis_method)

    pvals = res['ps']
    if analysis_method in ['lm', 'emma', 'emmax','amm']:
            additional_columns['genotype_var_perc'] = res['var_perc']
            if 'betas' in res:
                betas = map(list, zip(*res['betas']))
                additional_columns['beta0'] = betas[0]
                if len(betas) > 1:
                    additional_columns['beta1'] = betas[1]

    return GWASResult(genotypeData.chrs,genotypeData.chromosomes,genotypeData.positions,pvals,genotypeData.get_mafs(),method = analysis_method,transformation = phenData.transformation,additional_columns = additional_columns)


def transform_phenotype(args):
    phenotype_file = args['file']
    output_file = args['output']
    transformation = args['transformation']
    phenData = phenotype.parse_phenotype_file(phenotype_file)
    trans_type = phenData.transform(transformation)
    phenData.write_to_file(output_file)
    print trans_type
    return trans_type


def _save_ld_data(output_file,ld_data,chr_pos_list):
    _,ext = os.path.splitext(output_file)
    log.info('Saving LD data in %s ' % output_file)
    if ext not in SUPPORTED_FILE_EXT:
        raise Exception('The input file must have one of the supported extensions: %s' % str(SUPPORTED_FILE_EXT))
    if ext  in HDF5_FILE_EXT:
        _save_ld_data_as_hdf5(output_file,ld_data,chr_pos_list)
    else:
        _save_ld_data_as_csv(output_file,ld_data,chr_pos_list)

def _save_ld_data_as_hdf5(output_file,ld_data,chr_pos_list):
    f = h5py.File(output_file,'a',libvers='latest')
    chrs, positions = zip(*chr_pos_list)
    chr_regions,chrs = _get_chr_regions(chrs)
    min_chr_region = min(map(lambda x: x[1] - x[0], chr_regions))
    ld_group = f.create_group('ld')
    ld_snps = ld_group.create_dataset("ld_snps", (len(positions),),chunks=(min_chr_region,), dtype='i4',compression='lzf',data=positions)
    # required to store this way because JHDF5 does not support SIMPLE SPACE attribute arrays for strings
    chr_region_list = zip(map(str,chrs.tolist()),chr_regions)
    for chr, chr_region in chr_region_list:
        ld_snps.attrs["chr%s" % chr] = chr_region
    dt = h5py.special_dtype(vlen=np.dtype('float32'))
    ld_dset = ld_group.create_dataset('ld_data',(len(positions),),chunks=(min_chr_region,),compression='lzf', dtype=dt)
    for i in range(0, len(positions)):
        ld_dset[i] = ld_data[i][:i + 1]
    f.flush()
    f.close()

def _save_ld_data_as_csv(output_file,ld_data,chr_pos_list):
    with open(output_file,'w') as f:
        writer = csv.writer(f,delimiter=',')
        # write header
        header = ['']
        header.extend(map(lambda x: '%s-%s' % (x[0],x[1]),chr_pos_list))
        writer.writerow(header)
        for i,chr_pos in enumerate(chr_pos_list):
            row = [str(chr_pos)]
            row.extend(ld_data[i])
            writer.writerow(row)


def _get_chr_regions(chrs):
    grouped_chr = np.unique(chrs)
    chr_regions = []
    start = 0
    end = 0
    for chr in grouped_chr:
        end = np.searchsorted(chrs,chr)
        if end > 0:
            chr_regions.append((start,end))
        start = end
    chr_regions.append((start,len(chrs)))

    return chr_regions,grouped_chr

def _prepare_data_(genotypeData, phenData,transformation=None):
    """
    Coordinates phenotype and snps data for different mapping methods.
    """
    log.info('Converting replicates of phenotypes to averages')
    phenData.convert_to_averages()
    d = genotypeData.coordinate_w_phenotype_data(phenData)
    phenData.transform(transformation)
    return d['n_filtered_snps']


def _load_genotype_(folder):
    data_format = 'binary'
    hdf5_file = os.path.join(folder,'all_chromosomes_%s.hdf5' % data_format)
    if os.path.isfile(hdf5_file):
        return genotype.load_hdf5_genotype_data(hdf5_file)
    raise Exception('No Genotype files in %s folder were found.' % folder)

def _load_accessions(accession_file):
    if os.path.isfile(accession_file):
        with open(accession_file,'r') as f:
            accessions=[]
            reader = csv.reader(f)
            for row in reader:
                accessions.append(row[0])
            return accessions
    raise Exception('No accession file %s was found.' % accession_file)

def _get_kinship_file_(folder):
    return os.path.join(folder,'kinship_ibs_binary_mac5.h5py')


if __name__ == '__main__':
    sys.exit(main())

