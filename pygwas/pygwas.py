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
from core import result
from core import plotting
import numpy as np
import re
import os
import sys
import csv
import h5py
from operator import itemgetter
import ipdb

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

SUPPORTED_FILE_EXT =  ('.hdf5','.csv')

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
    subparsers = parser.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')
    analysis_parser = subparsers.add_parser('run',help='Run a GWAS analysis')

    analysis_parser.add_argument("-t", "--transformation", dest="transformation", help="Apply a transformation to the data. Default[None]", choices=["log", "sqrt", "exp", "sqr", "arcsin_sqrt", "box_cox"])
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
    plotter_parser.add_argument("-o",'--output',dest='output',required=True,help='The output image file')
    plotter_parser.add_argument(dest="file", help="GWAS result file (.hdf5 or .csv)", metavar="FILE")
    plotter_parser.set_defaults(func=plot)

    stats_parser = subparsers.add_parser('stats',help='Retrieve some stats')
    stats_parser.add_argument("-t", "--type", dest="type",required=True, help="type of the statistics to return",choices=["all","pseudo","shapiro"])
    stats_parser.add_argument("-g", "--genotype", dest="genotype_folder", help="folder with the genotypes for the GWAS analysis", required=True,metavar="FOLDER")
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
    phenotype_file = args['file']
    genotype_folder = args['genotype_folder']
    stat_type =args['type']
    phenData = phenotype.parse_phenotype_file(phenotype_file)  #load phenotype file
    genotypeData = _load_genotype_(genotype_folder)
    phenData.convert_to_averages()
    sd_indices_to_keep,pd_indices_to_keep = _get_indicies_(phenData.ecotypes,genotypeData.accessions)
    phenData.filter_ecotypes(pd_indices_to_keep)
    accessions = genotypeData.accessions[sd_indices_to_keep]
    K = None
    phen_vals = phenData.values
    kinship_file = args['kinship']
    if kinship_file is None:
        kinship_file = _get_kinship_file_(genotype_folder)
        K = kinship.load_kinship_from_file(kinship_file, accessions.tolist())['k']
    statistics = {}
    if stat_type == 'all':
        statistics['pseudo_heritability'] = gwas.calculate_pseudo_heritability(phen_vals,K)
        statistics['shapiro'] = stats.calculate_sp_pval(phen_vals)
    elif stat_type == 'pseudo':
        statistics['pseudo_heritability'] = gwas.calculate_pseudo_heritability(phen_vals,K)
    elif stat_type == 'kolmogorov':
        statistics['shapiro'] = stats.calculate_sp_pval(phen_vals)
    else: raise Exception('%s not supported' % stat_type)
    print statistics
    return statistics


def convert(args):

    _,input_ext = os.path.splitext(args['inputfile'])
    _,output_ext = os.path.splitext(args['outputfile'])
    if input_ext == output_ext:
        raise Exception('use different file extension for input (%s) and output file (%s)' % (input_ext,output_ext))
    if input_ext not in SUPPORTED_FILE_EXT:
        raise Exception('The input file must have one of the supported extensions: %s' % SUPPORTED_FILE_EXT)
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
        raise Exception('The input file must have one of the supported extensions: %s' % supported_extensions)
    chrs = None
    if 'chr' in args and args['chr'] is not None and args['chr'] != '':
        chrs = [args['chr']]
    gwas_result = None
    if ext == '.hdf5':
        gwas_result = result.load_from_hdf5(args['file'])
    else:
        gwas_result = result.load_from_csv(args['file'])
    plotting.plot_gwas_result(gwas_result,args['output'],chrs,args['macs'])


def run(args):
    _,ext = os.path.splitext(args['outputfile'])
    if ext not in SUPPORTED_FILE_EXT:
        raise Exception('The output file must have one of the supported extensions: %s' % supported_extensions)
    gwas_result = perform_gwas(args['file'],args['analysis_method'], args['genotype_folder'],args['transformation'],args['kinship'])
    if ext == '.hdf5':
        gwas_result.save_as_hdf5(args['outputfile'])
    else:
        gwas_result.save_as_csv(args['outputfile'])
    if args['calc_ld']:
        genotype_folder = args['genotype_folder']
        calc_ld_args = {'acession_file':args['file'],'positions':args['outputfile'],'genotype_file':os.path.join(genotype_folder,'all_chromosomes_binary.hdf5'),'range':2500}
        if ext == '.hdf5':
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
            raise Exception('The input file must have one of the supported extensions: %s' % supported_extensions)
        if ext  == '.hdf5':
            gwas_result = result.load_from_hdf5(positions)
        else:
            gwas_result =  result.load_from_csv(positions)          
        gwas_data = gwas_result.get_top_snps(range)
        chr_pos_list = zip(map(str,gwas_data['chr']),gwas_data['positions'])
    chr_pos_list = sorted(chr_pos_list,key=itemgetter(0,1)) 
    ld_data = genotypeData.calculate_ld(chr_pos_list)
    _save_ld_data(args['output_file'],ld_data,chr_pos_list)   
    

def perform_gwas(phenotype_file,analysis_method,genotype_folder,transformation=None,kinshipFile=None):
    phenData = phenotype.parse_phenotype_file(phenotype_file)  #load phenotype file
    additional_columns = {}
    genotypeData = _load_genotype_(genotype_folder)
    K = None
    n_filtered_snps = _prepare_data_(genotypeData,phenData)
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

    return GWASResult(genotypeData.chrs,genotypeData.chromosomes,genotypeData.positions,pvals,genotypeData.get_mafs(),method = analysis_method,transformation = transformation,additional_columns = additional_columns)


def _save_ld_data(output_file,ld_data,chr_pos_list):
    _,ext = os.path.splitext(output_file)
    log.info('Saving LD data in %s ' % output_file)
    if ext not in SUPPORTED_FILE_EXT:
        raise Exception('The input file must have one of the supported extensions: %s' % supported_extensions)
    if ext  == '.hdf5':
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

def _prepare_data_(genotypeData, phenData,with_replicates=False):
    """
    Coordinates phenotype and snps data for different mapping methods.
    """
    if not with_replicates:
        log.info('Converting replicates of phenotypes to averages')
        phenData.convert_to_averages()
    d = genotypeData.coordinate_w_phenotype_data(phenData)
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
