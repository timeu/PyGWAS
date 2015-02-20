#!/usr/bin/env python
# coding: utf-8
"""
    plotter
    ~~~~~~~~~~~~~

    plotter is a wrapper script for plotting GWAS results

    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""

from __init__ import __version__,__updated__,__date__
from core import plotting
import sys,os
import argparse
import logging
import logging.config

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
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument('-c','--chr',dest='chr',
	 help='Chromosome to plot. If not specified prints all chromosomes (Default:None)',
	choices=['chr1','chr2','chr3','chr4','chr5'],default=None)
    parser.add_argument('-m','--macs',dest='macs',default=15,type=int,help='Minor Allele Count filter (Default: 15)')
    parser.add_argument("-o",'--output',dest='output',required=True,help='The output image file')
    parser.add_argument(dest="file", help="HDF5 GWAS result file", metavar="FILE")
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
        chrs = None
        if 'chr' in args and args['chr'] is not None and args['chr'] != '':
            chrs = [args['chr']]
        plotting.plot_gwas_result(args['file'],args['output'],chrs,args['macs'])
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        log.exception(e)
        return 2

if __name__ == "__main__":
    sys.exit(main())
