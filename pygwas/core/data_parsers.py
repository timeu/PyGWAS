import logging
import csv
import numpy
import scipy
import numpy as np

log = logging.getLogger(__name__)


#Standard missing value is N (Used to be NA)
missing_val = 'N'
#Standard nt decoding, is using the IUPAC alphabet
nt_decoder = {'A':'A',
          'C':'C',
          'G':'G',
          'T':'T',
          'AG':'R',
          'AC':'M',
          'GT':'K',
          'CT':'Y',
          'AT':'W',
          'CG':'S',
          'Y':'Y',
          'R':'R',
          'W':'W',
          'S':'S',
          'K':'K',
          'M':'M',
          'D':'D',
          'H':'H',
          'V':'V',
          'B':'B',
          'X':'X', #Unknown base(s)
          'N':'X', #Unknown base(s)
          '-':'-', #Indel 
          '|':missing_val}


#An int decoder is useful for processing the data efficiently
nt_int_decoder = {'A':1,
          'C':2,
          'G':3,
          'T':4,
          'AG':5,
          'AC':6,
          'GT':7,
          'CT':8,
          'AT':9,
          'CG':10,
          'Y':11,
          'R':12,
          'W':13,
          'S':14,
          'K':15,
          'M':16,
          'D':17,
          'H':18,
          'V':19,
          'B':20,
          'X':21, #Unknown base(s)
          'N':21, #Unknown base(s)
          '-':22, #Indel 
          '|':0}


def parse_genotype_csv_file(csv_file,format='binary', missingVal='NA'):
    log.info('Parsing Genotype file (CSV) %s in format %s' % (csv_file, format))
    retval ={}

    with open(csv_file) as f:
        dialect = csv.Sniffer().sniff(f.read(40))
        f.seek(0)
        reader = csv.reader(f, dialect)
        header = reader.next()
        if header[0] != 'Chromosome' or (header[1] != 'Positions' and header[1] != 'Position'): 
            raise Exception('First two columns must be in form  Chromosome, Positions')
        dtype = _get_dtype_from_format(format)
        accessions = map(lambda x: x.strip(),header[2:])
        log.info('%s accessions found %s' % (len(accessions),accessions))
        first_line = reader.next()
        snps = []
        positions = []
        old_chr = first_line[0]
        positions.append(int(first_line[1]))
        snps.append(_get_snps_from_row(format,first_line[2:]))
        data = []
        retval = {'format':format,"accessions":accessions}
        for row in reader:
            chr = row[0]
            if chr != old_chr: #Then save snps_data
                log.info('Chromosome %s finished. %s SNPs found' % (old_chr,len(snps)))
                data.append({'snps':snps,'positions':positions,'chr':old_chr})
                positions = []
                snps = []
                old_chr = chr
            positions.append(int(row[1]))
            snps.append(_get_snps_from_row(format,row[2:]))
        data.append({'snps':snps,'positions':positions,'chr':chr})
    retval['data'] = data
    log.info('Finished parsing Genotype file')
    return retval     
   

def _get_snps_from_row(format,snps):
    dtype = _get_dtype_from_format(format)
    if format == 'S1':
        map(lambda x:nt_decoder[x],snps)
    return np.array(snps,dtype=dtype)

def _get_dtype_from_format(format):
   dtype = 'int8'
   if format == 'float':
       dtype = 'float32'
   elif format == 'nucleotides':
       dtype ='S1'
   return dtype  
        


