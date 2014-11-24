#!/usr/bin/env python
# coding: utf-8
"""
    pygwas.pygwas
    ~~~~~~~~~~~~~

    The main module for running Genome Wide Association studies

    :copyright: year by my name, see AUTHORS for more details
    :license: license_name, see LICENSE for more details
"""

from __init__ import __version__,__updated__,__date__
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from core import kinship
from core import gwas
import logging
from core import mtcorr
from core import statistics as stats
from core import phenotype
from core import genotype
import os
import sys


FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
logging.basicConfig(format=FORMAT,level=logging.DEBUG)
log = logging.getLogger()

def get_parser(program_license,program_version_message):
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-t", "--transformation", dest="transformation", help="Apply a transformation to the data. Default[None]", choices=["log", "sqrt", "exp", "sqr", "arcsin_sqrt", "box_cox"])
    parser.add_argument("-a", "--analysis_method", dest="analysis_method", help="analyis method to use",required=True,choices=["lm", "emma", "emmax", "kw", "ft", "emmax_anova", "lm_anova", "emmax_step", "lm_step","loc_glob_mm","amm"])
    parser.add_argument("-k", "--kinship", dest="kinship", help="Specify the file containing the kinship matrix. (otherwise default file is used or it's generated.)", metavar="FILE" )
    parser.add_argument("-s", "--kinship_type", dest="kinship_type", help="Type of kinship calculated. Possible types are ibs (default) or ibd ", choices=["ibs", "ibd"],default="ibs")
    parser.add_argument("-o", "--output_file", dest="outputfile", help="Name of the output file")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument(dest="file", help="csv file containing phenotype values", metavar="FILE")
    return parser


def main(): 
    '''Command line options.'''

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by Ümit Seren on %s.
  Copyright 2012 Gregor Mendel Institute. All rights reserved.
  
  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0
  
  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Process arguments
        parser = get_parser(program_license,program_version_message)
        args = vars(parser.parse_args())
        messenger = StdoutMessenger()
        if args.queue: 
            messenger = ProgressMessenger(args.queue_host,5672,'admin','eastern')
        
        messenger.update_status(progress=0.0, task_status='Loading phenotype data')
        phenData = phenotypeData.parse_phenotype_file(args.file,False)  #load phenotype file
        phen_ids = phenData.phen_dict.keys()  # get phenotype ids
        #If not analysis plots... then GWAS
        for phen_id in phen_ids:
            phenotype_name = phenData.get_name(phen_id)
            messenger.update_status(progress=0.0, task_status='Loading phenotype data')
            print "Performing GWAS for phenotype: %s, phenotype_id: %s" % (phenotype_name, phen_id)
            _perform_gwas_(phen_id, phenData, args.analysis_method, args.transformation,args.genotype,args.kinship_type,args.kinship,messenger,args.outputfile)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


if __name__ == '__main__':
    sys.exit(main())
