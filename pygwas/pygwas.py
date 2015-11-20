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
import os
import sys


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

SUPPORTED_FILE_EXT =  ('.hdf5','.csv')

logging.config.dictConfig(LOGGING)
log = logging.getLogger()

def get_parser(program_license,program_version_message):
    parser = argparse.ArgumentParser(description=program_license)
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    subparsers = parser.add_subparsers(title='subcommands',description='Choose a command to run',help='Following commands are supported')
    analysis_parser = subparsers.add_parser('run',help='Run a GWAS analysis')

    analysis_parser.add_argument("-t", "--transformation", dest="transformation", help="Apply a transformation to the data. Default[None]", choices=["log", "sqrt", "exp", "sqr", "arcsin_sqrt", "box_cox"])
    analysis_parser.add_argument("-a", "--analysis_method", dest="analysis_method", help="analyis method to use",required=True,choices=["lm", "emma", "emmax", "kw", "ft", "emmax_anova", "lm_anova", "emmax_step", "lm_step","loc_glob_mm","amm"])
    analysis_parser.add_argument("-g", "--genotype", dest="genotype_folder", help="folder with the genotypes for the GWAS analysis", required=True,metavar="FOLDER")
    analysis_parser.add_argument("-k", "--kinship", dest="kinship", help="Specify the file containing the kinship matrix. (otherwise default file is used or it's generated.)", metavar="FILE" )
    analysis_parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file")
    analysis_parser.add_argument(dest="file", help="csv file containing phenotype values", metavar="FILE")
    analysis_parser.set_defaults(func=run)

    convert_parser = subparsers.add_parser('convert',help='Convert from HDF5 to CSV and the other way round')
    convert_parser.add_argument("-i", "--input_file", dest="inputfile", help="Name of the input file",required=True)
    convert_parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file",required=True)
    convert_parser.set_defaults(func=convert)

    plotter_parser = subparsers.add_parser('plot',help='Plot a GWAS result')
    plotter_parser.add_argument('-c','--chr',dest='chr',
	 help='Chromosome to plot. If not specified prints all chromosomes (Default:None)',
	choices=['chr1','chr2','chr3','chr4','chr5'],default=None)
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
    raise Exception('No Genotype files in %s folder were found.' % file_prefix)

def _get_kinship_file_(folder):
    return os.path.join(folder,'kinship_ibs_binary_mac5.h5py')


if __name__ == '__main__':
    sys.exit(main())
