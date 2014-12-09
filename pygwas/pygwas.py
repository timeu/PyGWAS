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

logging.config.dictConfig(LOGGING)
log = logging.getLogger()

def get_parser(program_license,program_version_message):
    parser = argparse.ArgumentParser(description=program_license)
    parser.add_argument("-t", "--transformation", dest="transformation", help="Apply a transformation to the data. Default[None]", choices=["log", "sqrt", "exp", "sqr", "arcsin_sqrt", "box_cox"])
    parser.add_argument("-a", "--analysis_method", dest="analysis_method", help="analyis method to use",required=True,choices=["lm", "emma", "emmax", "kw", "ft", "emmax_anova", "lm_anova", "emmax_step", "lm_step","loc_glob_mm","amm"])
    parser.add_argument("-g", "--genotype", dest="genotype", help="genotype dataset to be used in the GWAS analysis (run with option -l to display list of available genotype datasets)", required=True, type=int,metavar="INTEGER" )
    parser.add_argument("-f", "--genotype_folder", dest="genotype_folder", help="Folder where the genotypes are located",metavar="FOLDER")

    parser.add_argument("-k", "--kinship", dest="kinship", help="Specify the file containing the kinship matrix. (otherwise default file is used or it's generated.)", metavar="FILE" )
    parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument(dest="file", help="csv file containing phenotype values", metavar="FILE")
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
        result = perform_gwas(args['file'],args['analysis_method'], args['genotype_folder'],args['genotype'],args['transformation'],args['kinship'])
        result.save_as_hdf5(args['outputfile'])
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        log.exception(e)
        return 2


def perform_gwas(phenotype_file,analysis_method,genotype_folder,genotype,transformation=None,kinshipFile=None):
    phenData = phenotype.parse_phenotype_file(phenotype_file)  #load phenotype file
    additional_columns = {}
    genotypeData = _load_genotype_(genotype_folder,genotype)
    K = None
    n_filtered_snps = _prepare_data_(genotypeData,phenData)
    phen_vals = phenData.values
    if analysis_method in ['emma', 'emmax', 'emmax_anova', 'emmax_step', 'loc_glob_mm','amm']:
        #Load genotype file (in binary format)
        log.debug("Retrieving the Kinship matrix K.\n")
        if kinshipFile is None:   #Kinship file was supplied..
	    kinshipFile = _get_kinship_file_(genotype_folder,genotype)
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
    
    return GWASResult(pvals,None,analysis_method,transformation,genotypeData,additional_columns)



def _prepare_data_(genotypeData, phenData,with_replicates=False):
    """
    Coordinates phenotype and snps data for different mapping methods.
    """
    if not with_replicates:
        log.info('Converting replicates of phenotypes to averages')
        phenData.convert_to_averages()
    d = genotypeData.coordinate_w_phenotype_data(phenData)
    return d['n_filtered_snps']

def _get_folder_(folder,genotype_id):
    return os.path.join(folder, str(genotype_id))

def _load_genotype_(folder,genotype_id):
    data_format = 'binary'
    file_prefix = _get_folder_(folder,genotype_id)
    hdf5_file = os.path.join(file_prefix,'all_chromosomes_%s.hdf5' % data_format)
    if os.path.isfile(hdf5_file):
        return genotype.load_hdf5_genotype_data(hdf5_file)
    raise Exception('No Genotype files in %s folder were found.' % file_prefix)

def _get_kinship_file_(folder,genotype_id):
    return os.path.join(_get_folder_(folder,genotype_id),'kinship_ibs_binary_mac5.h5py') 


if __name__ == '__main__':
    sys.exit(main())
