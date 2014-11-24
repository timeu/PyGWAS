import logging
import kinship
import bisect
import itertools as iter
import phenotype
import data_parsers
import h5py
import numpy
import scipy
import pdb
from abc import ABCMeta, abstractmethod, abstractproperty

log = logging.getLogger(__name__)


def load_hdf5_genotype_data(hdf5_file):
    return HDF5Genotype(hdf5_file)


def load_csv_genotype_data(csv_files,format='binary'):
    log.info("Loading Genotype file")
    snps = []
    positions = []
    chr_regions = []
    chrs = []
    accessions = None
    if type(csv_files) is list:
        log.info("Genotype split in %s files" % len(csv_file))
        for i,csv_file in enumerate(csv_files):
            log.info("Loading %s " % csv_file)
            data = data_parser.parse_genotype_csv_file(csv_file,format)
            if accessions == None:
                accessions = data['accessions']
            if data['accessions'] != accessions:
                raise Exception('Accessions must match')
            for chr_data in data['data']:
                num_snps = len(chr_data['positions'])
                chr_regions.append((len(positions),len(positions)+num_snps)) 
                chrs.append(chr_data['chr'])
                positions.extend(chr_data['positions'])
                snps.extend(chr_data['snps'])
                log.info("Finished loading %s SNPs for Chr %s" % (len(positions),chr_data['chr']))
    else:
        data = data_parsers.parse_genotype_csv_file(csv_files,format)
        accessions = data['accessions']
        for chr_data in data['data']:
            num_snps = len(chr_data['positions'])
            chr_regions.append((len(positions),len(positions)+num_snps)) 
            chrs.append(chr_data['chr'])
            positions.extend(chr_data['positions'])
            snps.extend(chr_data['snps'])
            log.info("Finished loading %s SNPs for Chr %s" % (len(positions),chr_data['chr'])) 
    log.info("Finished loading Genotype file: %s SNPs %s accessions %s chromosomes" % (len(positions),len(accessions),len(chrs)))
    return Genotype(snps,positions, accessions,chr_regions,chrs,format)




class AbstractGenotype(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def data_format(self):
        pass

    @abstractmethod
    def get_snps_iterator(self,chr=None,is_chunked=False):
        pass

    @abstractmethod
    def get_snp_at(self,chr,position):
        pass

    @abstractproperty
    def positions(self):
        pass

    @abstractproperty
    def accessions(self):
        pass

    @property
    def chromosomes(self):
        chromosomes = []
        for i,chr_region in enumerate(self.chr_regions):
			chromosomes.extend([(i+1)] * (chr_region[1] - chr_region[0]))			
        return chromosomes

    @abstractproperty
    def chrs(self):
        pass

    @abstractproperty
    def num_snps(self):
        pass

    @abstractproperty
    def chr_regions(self):
        pass

    def get_chr_pos_from_index(self,ix):
        return (self.get_chr_from_index(ix),self.positions[ix])

    def get_chr_from_index(self,ix):
        for i,chr_region in enumerate(self.chr_regions):
            if chr_region[0] <= ix and chr_region[1] >= ix:
                return i+1
        raise Exception('Index %s outside of chr_regions' %ix)

    def get_mafs(self):
        macs = []
        mafs = []
        num_nts = len(self.accessions)
        for snp in self.get_snps_iterator():
            l = scipy.bincount(snp)
            mac = min(l)
            macs.append(mac)
            mafs.append(mac / float(num_nts))
        return {"macs":macs, "mafs":mafs}

    @abstractproperty
    def genome_length(self):
        pass

    @abstractmethod
    def filter_accessions_ix(self,accessions_ix): 
        pass

    @abstractmethod
    def filter_snps_ix(self,accessions_ix): 
        pass
  

    def coordinate_w_phenotype_data(self, phenotype, coord_phen=True):
        """
        Deletes accessions which are not common, and sorts the accessions, removes monomorphic SNPs, etc.
        """
        log.debug("Coordinating SNP and Phenotype data.")
        ets = phenotype.ecotypes
        sd_indices_to_keep = set()
        pd_indices_to_keep = []
        for i, acc in enumerate(self.accessions):
            for j, et in enumerate(ets):
                if str(et) == str(acc):
                    sd_indices_to_keep.add(i)
                    pd_indices_to_keep.append(j)

        sd_indices_to_keep = list(sd_indices_to_keep)
        sd_indices_to_keep.sort()
        #Filter accessions which do not have phenotype values (from the genotype data).
        log.debug("Filtering accessions")
        self.filter_accessions_ix(sd_indices_to_keep)

        if coord_phen:
            num_values = phenotype.num_vals
            log.debug("Filtering phenotype data.")
            #Removing accessions that don't have genotypes or phenotype values
            phenotype.filter_ecotypes(pd_indices_to_keep) 
            ets = phenotype.ecotypes
            log.debug("Out of %d, leaving %d values." % (num_values, len(ets)))
        if self.data_format == 'binary':
            log.debug('Filtering non-binary SNPs')
            (total_num,removed_num) = self.filter_non_binary()
            log.debug('Removed %d non-binary SNPs out of %d SNPs' % (removed_num, total_num))
        elif self.data_format in ['int', 'diploid_int']:
            log.debug('Filtering monomorhpic SNPs')
            total_num = 0
            removed_num = 0
            (total_num,removed_num) = self.filter_monomorphic_snps()
            log.debug('Removed %d monomorphic SNPs out of %d SNPs' % (removed_num, total_num))
        return {'pd_indices_to_keep':pd_indices_to_keep, 'n_filtered_snps':removed_num}

    def get_ibs_kinship_matrix(self, debug_filter=1, snp_dtype='int8', dtype='single'):
        """
        Calculate the IBS kinship matrix. 
        (un-scaled)
        
        Currently it works only for binary kinship matrices.
        """
        log.debug('Starting kinship calculation')
        snps = self.getSnps()
        return kinship.calc_ibs_kinship(snps, snps_data_format=self.data_format, snp_dtype=snp_dtype,
                                        dtype=dtype)

    def get_ibd_kinship_matrix(self, debug_filter=1, dtype='single'):
        log.debug('Starting IBD calculation')
        snps = self.getSnps(debug_filter)
        cov_mat = kinship.calc_ibd_kinship(snps, len(self.accessions), dtype=dtype)
        log.debug('Finished calculating IBD kinship matrix')
        return cov_mat

   

    def save_as_csv(self):
        pass


    def save_as_hdf5(self, hdf5_file):
        log.info('Writing genotype to HDF5 file %s' %hdf5_file)
        if self.data_format in ['binary', 'diploid_int']:
            h5file = h5py.File(hdf5_file, 'w')
            num_snps = self.num_snps
            num_accessions = len(self.accessions)
            h5file.create_dataset('accessions', data=self.accessions,shape=(num_accessions,))
            h5file.create_dataset('positions', data=self.positions,shape=(num_snps,),dtype='i4')
            h5file['positions'].attrs['chrs'] = self.chrs
            h5file['positions'].attrs['chr_regions']  = self.chr_regions
            h5file.create_dataset('snps', shape=(num_snps, len(self.accessions)),
                              dtype='int8', compression='lzf', chunks=((1000, num_accessions)),data=list(self.get_snps_iterator()))
            h5file['snps'].attrs['data_format'] = self.data_format
            h5file['snps'].attrs['num_snps'] = num_snps
            h5file['snps'].attrs['num_accessions'] = num_accessions
            h5file.close()
            log.info("Finished writing genotype to HDF5 file")
        else:
            raise NotImplementedError
 

    def filter_monomorphic_snps(self):
        """
        Removes SNPs from the data which are monomorphic.
        """
        snps_ix = []
        num_snps = self.num_snps
        for i,snps in enumerate(self.get_snps_iterator()):
            bc = scipy.bincount(snps)
            if len(bc) > 1:
                snps_ix.append(i)
        numRemoved = len(self.positions) - len(snps_ix)
        self.filter_snps_ix(snps_ix)
        log.info("Removed %d monomoprhic SNPs, leaving %d SNPs in total." % (numRemoved, num_snps))
        return (num_snps,numRemoved)

    def filter_non_binary(self):
        """
        Removes all but binary SNPs.  (I.e. monomorphic, tertiary and quaternary alleles SNPs are removed.)
        """
        num_snps = self.num_snps
        snps_ix = []
        for i,snps in enumerate(self.get_snps_iterator()):
            bc = scipy.unique(snps)
            if len(bc) != 2:
                snps_ix.append(i)
        numRemoved = len(snps_ix)
        self.filter_snps_ix(snps_ix)
        log.info("Removed %d non-binary SNPs, leaving %d SNPs in total." % (numRemoved, self.num_snps))
        return (num_snps,numRemoved)
    

class Genotype(AbstractGenotype):

    """
    A class that encompasses multiple _SnpsData_ chromosomes objects (chromosomes), and can deal with them as a whole.

    This object should eventually replace the snpsdata lists.. 
    """

    def __init__(self, snps,positions, accessions,chr_regions, chrs,data_format=None):
        self._snps = snps
        self._positions = positions
        self._chr_regions = chr_regions
        self._chrs = chrs
        self._accessions = accessions
        self._data_format = data_format # binary, diploid_ints, floats, int


    @property
    def data_format(self):
        return self._data_format

    def get_snps_iterator(self,chr=None,is_chunked=False):
        for snp in self._snps:
            yield snp

    def get_snp_at(self,chr,position):
        chr_ix = self._chrs.index(chr)
        chr_start_ix = self.chr_regions[chr_ix][0]
        snp_ix = bisect.bisect(self._positions, position) - 1
        return self.snps[snp_ix]

    @property
    def accessions(self):
        return self._accessions

    
    @property
    def positions(self):
        return self._positions

    @property
    def chrs(self):
        return self._chrs

    @property
    def num_snps(self):
        return len(self._snps)

    @property
    def chr_regions(self):
        return self._chr_regions

    @property
    def genome_length(self):
        return self.chr_regions[-1][0]

    def filter_snps_ix(self,snps_ix): 
        self._snps = snps[snps_ix]
        self._positions = snps[snps_ix]

 

    def filter_accessions_ix(self,indicesToKeep): 
        """
        Removes accessions from the data.
        """
        num_accessions  = len(self.accessions)
        newAccessions = []
        newArrayIds = []
        for i in indicesToKeep:
            newAccessions.append(self.accessions[i])
        for i in range(len(self.snps)):
            snp = self.snps[i]
            newSnp = []
            for j in indicesToKeep:
                newSnp.append(snp[j])
            self.snps[i] = newSnp
        self.accessions = newAccessions
        log.debug("Removed %d accessions, leaving %d in total." % (num_accessions - len(indicesToKeep), len(indicesToKeep)))


class HDF5Genotype(AbstractGenotype):

    def __init__(self,hdf5_file): 
        self.h5file = h5py.File(hdf5_file, 'r')
        self.filter_snps = None
        self.filter_accessions = None

    def __del__(self):
        if self.h5file != None:
           self.h5file.close()

    def get_snps(self):
        return self.h5file['snps']

    @property
    def data_format(self):
        return self.h5file['snps'].attrs['data_format']

    
    def get_snp_at(self,chr,position):
        chr_ix = self._chrs.index(chr)
        chr_start_ix = self.chr_regions[chr_ix][0]
        snp_ix = bisect.bisect(self._positions, position) - 1
        return self.snps[snp_ix]

    @property
    def accessions(self):
        if self.filter_accessions == None or len(self.filter_accessions) == 0:
            return self.h5file['accessions'][:]
        else:
            return self.h5file['accessions'][self.filter_accessions]


    def _get_snps_(self, start=0,end=None,chunk_size=1000): 
        """
        An generator/iterator for SNP chunks.
        """
        if end is None:
            end = self.original_num_snps
        for i in xrange(start,end,chunk_size):
            stop_i = min(i + chunk_size, end)
            if self.filter_accessions == None or len(self.filter_accessions) == 0:
                if self.filter_snps == None:
                    yield self.h5file['snps'][i:stop_i]
                else:
                    yield self.h5file['snps'][self.filter_snps[i:stop_i],:]
            else:
                if self.filter_snps == None:
                    yield self.h5file['snps'][i:stop_i,self.filter_accessions]
                else:
                    filter_chunk = self.filter_snps[i:stop_i]
                    snps_chunk = self.h5file['snps'][i:stop_i, self.filter_accessions]
                    yield snps_chunk[filter_chunk]

    def get_snps_iterator(self,chr=None,is_chunked=False):
        """
        Returns an generator containing a chunked generator. 
        If chr is passed the generator will only iterate over the specified chr
        """
        start = 0
        end = None
        if chr != None:
			# use unfiltered chr_regions because filtering happens in _get_snps_
            chr_region = self.h5file['positions'].attrs['chr_regions'][chr-1]
            start = chr_region[0]
            end = chr_region[1]
        for snp_chunk in self._get_snps_(start=start,end=end):
            if is_chunked:
                yield snp_chunk
            else:
                for snp in snp_chunk:
                    yield snp

    @property
    def original_num_snps(self):
        return self.h5file['snps'].attrs['num_snps'] 
   
    @property
    def positions(self):
        if self.filter_snps != None:
            return self.h5file['positions'][self.filter_snps]
        return self.h5file['positions'][:]

    def get_snps_from_pos(self,chr_pos):
        """
        Returns a list of snps based on a list of (chr,pos) tuples
        """
        snps = []
        indices = []
        if chr_pos == None or len(chr_pos) == 0:
            return (indices,snps)
        chr_pos_ix = map(list,zip(*sorted(zip(chr_pos,range(len(chr_pos))))))
        # group by chr for efficient sequentielly iteration over snps generator
        for chr,positions in iter.groupby(chr_pos_ix[0],lambda x:x[0]):
            pos_indices = []
            it = self.get_snps_iterator(chr)
            for position in positions:
                pos_ix = self._get_pos_ix_(chr,position[1])
                pos_indices.append(pos_ix[1])
                indices.append(pos_ix[0])
            for i,ix in enumerate(pos_indices):
                previous_ix = 0 if i == 0 else pos_indices[i-1] +1 
                snps.append(next(iter.islice(it,ix-previous_ix,None),None))
        return (map(list,zip(*sorted(zip(chr_pos_ix[1],indices))))[1],map(list,zip(*sorted(zip(chr_pos_ix[1],snps))))[1])


    def _get_pos_ix_(self,chr,position):
        """
        Returns the index of chr,position using bisect
        """
        chr_region = self.chr_regions[(chr-1)]
        positions = self.positions[chr_region[0]:chr_region[1]]
        i = bisect.bisect(positions, position) - 1
        return (chr_region[0] + i,i)
    

    @property
    def chrs(self):
        return self.h5file['positions'].attrs['chrs']

    @property
    def num_snps(self):
        if self.filter_snps != None:
            return self.filter_snps.sum()
        return self.original_num_snps
        
    @property
    def chr_regions(self):
        if self.filter_snps != None:
            return self.filtered_chr_regions
        return self.h5file['positions'].attrs['chr_regions']

    @property
    def genome_length(self):
        return self.chr_regions[-1][0]


    def filter_accessions_ix(self,indicesToKeep): 
        """
        Removes accessions from the data.
        """
        num_accessions = len(self.accessions)
        self.filter_accessions = indicesToKeep
        log.debug("Removed %d accessions, leaving %d in total." % (num_accessions - len(indicesToKeep), len(indicesToKeep)))


    def filter_snps_ix(self,snps_ix): 
        if snps_ix == None or len(snps_ix) == 0:
            self.filter_snps = None
            self.filtered_chr_regions = None
        else:
            self.filter_snps = numpy.ones((self.num_snps,),dtype=numpy.bool)
            self.filter_snps[snps_ix] = 0
            self.filtered_chr_regions = self._get_filtered_regons()


    def _get_filtered_regons(self):
        filtered_chr_regions = []
        if self.filter_snps == None: 
            return None
        start_ix = 0
        end_ix = 0
        for chr_region in self.h5file['positions'].attrs['chr_regions']:
            end_ix = start_ix + self.filter_snps[chr_region[0]:chr_region[1]].sum()
            filtered_chr_regions.append((start_ix,end_ix))
            start_ix = end_ix
        return filtered_chr_regions
            
            




   






















