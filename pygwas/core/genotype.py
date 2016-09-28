import logging
import kinship
import bisect
import itertools as iter
import phenotype
import data_parsers
import h5py
import numpy
import scipy
import csv
import pdb
import itertools
from operator import itemgetter
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
        log.info("Genotype split in %s files" % len(csv_files))
        for i,csv_file in enumerate(csv_files):
            log.info("Loading %s " % csv_file)
            data = data_parsers.parse_genotype_csv_file(csv_file,format)
            if accessions is None:
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

def calculate_ld(snps):
        #filter non binary snps
        snps_t = scipy.transpose(snps)
        snps_stand = scipy.transpose((snps_t - scipy.mean(snps, 1)) / scipy.std(snps, 1))
        r2_values =scipy.dot(snps_stand, scipy.transpose(snps_stand))
        r2_values *= (1.0 / snps.shape[1])
        r2_values **= 2
        return r2_values


class AbstractGenotype(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def data_format(self):
        pass

    @abstractmethod
    def get_snps_iterator(self,chr=None,is_chunked=False,chunk_size=1000):
        pass

    def get_snp_at(self,chr,position):
        chr_ix = chr -1
        chr_region = self.chr_regions[chr_ix]
        i = bisect.bisect_left(self.positions[chr_region[0]:chr_region[1]], position)
        if i != chr_region[1]:
            snp_ix = i + chr_region[0]
            if self.positions[snp_ix] == position:
                return self.snps[snp_ix]
        return None

    def get_chr_region_ix(self,chr):
        # some genotype datasets have chromosomes as integers instead of strings
        if self.chrs.dtype.kind == 'i':
            return numpy.where(self.chrs == int(chr))[0][0]
        return numpy.where(self.chrs == str(chr))[0][0]

    @abstractproperty
    def positions(self):
        pass

    @abstractproperty
    def accessions(self):
        pass

    @abstractproperty
    def snps(self):
        pass

    @property
    def chromosomes(self):
        chromosomes = []
        for i,chr_region in enumerate(self.chr_regions):
			chromosomes.extend([self.chrs[i]] * (chr_region[1] - chr_region[0]))
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

    def get_snps_from_pos(self,chr_pos):
        """
        Returns a list of snps based on a list of (chr,pos) tuples
        """
        snps = []
        indices = []
        if chr_pos is None or len(chr_pos) == 0:
            return (indices,snps)
        chr_pos_ix = map(list,zip(*sorted(zip(chr_pos,range(len(chr_pos))))))
        # group by chr for efficient sequentielly iteration over snps generator
        it_ix = 0
        filtered_chr_pos_ix = [[],[]]
        for chr,positions in iter.groupby(chr_pos_ix[0],lambda x:x[0]):
            pos_indices = []
            it = self.get_snps_iterator(chr)
            for position in positions:
                pos_ix = self.get_pos_ix(chr,position[1])
                if pos_ix[2] == False:
                    it_ix+=1
                    continue
                filtered_chr_pos_ix[0].append(chr_pos_ix[0][it_ix])
                filtered_chr_pos_ix[1].append(chr_pos_ix[1][it_ix])
                pos_indices.append(pos_ix[1])
                indices.append(pos_ix[0])
                it_ix+=1
            for i,ix in enumerate(pos_indices):
                previous_ix = 0 if i == 0 else pos_indices[i-1] +1
                snps.append(next(iter.islice(it,ix-previous_ix,None),None))
        return (map(list,zip(*sorted(zip(filtered_chr_pos_ix[1],indices))))[1],map(list,zip(*sorted(zip(filtered_chr_pos_ix[1],snps))))[1])


    def get_pos_ix(self,chr,position):
        """
        Returns the index of chr,position using bisect
        """
        chr_region = self.chr_regions[self.get_chr_region_ix(chr)]
        positions = self.positions[chr_region[0]:chr_region[1]]
        i = bisect.bisect_left(positions, position)
        abs_ix = i + chr_region[0]
        found = False
        if abs_ix != len(self.positions) and self.positions[abs_ix] == position:
            found = True
        return (abs_ix,i,found)


    def get_chr_region_from_index(self,ix):
        for i,chr_region in enumerate(self.chr_regions):
            if chr_region[0] <= ix and chr_region[1] >= ix:
                return i, chr_region
        return None

    def get_chr_pos_from_index(self,ix):
        return (self.get_chr_from_index(ix),self.positions[ix])

    def get_chr_from_index(self,ix):
        for i,chr_region in enumerate(self.chr_regions):
            if chr_region[0] <= ix and chr_region[1] >= ix:
                return int(self.chrs[i])
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


    def filter_accessions(self,accession_ids):
        sd_indices_to_keep = set()
        pd_indices_to_keep = []
        for i, acc in enumerate(self.accessions):
            for j, et in enumerate(accession_ids):
                if str(et) == str(acc):
                    sd_indices_to_keep.add(i)
                    pd_indices_to_keep.append(j)

        sd_indices_to_keep = list(sd_indices_to_keep)
        sd_indices_to_keep.sort()
        self.filter_accessions_ix(sd_indices_to_keep)
        return sd_indices_to_keep,pd_indices_to_keep


    @abstractmethod
    def filter_accessions_ix(self,accessions_ix):
        pass

    @abstractmethod
    def filter_snps_ix(self,accessions_ix):
        pass

    @abstractmethod
    def convert_data_format(self,target_format='binary'):
        pass

    def coordinate_w_phenotype_data(self, phenotype, coord_phen=True):
        """
        Deletes accessions which are not common, and sorts the accessions, removes monomorphic SNPs, etc.
        """
        log.debug("Coordinating SNP and Phenotype data.")
        ets = phenotype.ecotypes

        #Filter accessions which do not have phenotype values (from the genotype data).
        log.debug("Filtering accessions")
        sd_indices_to_keep,pd_indices_to_keep =  self.filter_accessions(ets)

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

    def get_ibs_kinship_matrix(self, debug_filter=1, snp_dtype='int8', dtype='single',chunk_size=None):
        """
        Calculate the IBS kinship matrix.
        (un-scaled)

        Currently it works only for binary kinship matrices.
        """
        log.debug('Starting kinship calculation')
        return kinship.calc_ibs_kinship(self,chunk_size=chunk_size)

    def get_ibd_kinship_matrix(self, debug_filter=1, dtype='single',chunk_size=None):
        log.debug('Starting IBD calculation')
        return kinship.calc_ibd_kinship(self,chunk_size=chunk_size)
        log.debug('Finished calculating IBD kinship matrix')
        return cov_mat



    def save_as_csv(self,csv_file,chunk_size=1000):
        log.info('Writing genotype to CSV file %s' % csv_file)
        with open(csv_file,'w') as csvfile:
            csv_writer = csv.writer(csvfile,delimiter=',')
            header = ['Chromosome','Positions']
            header.extend(self.accessions)
            csv_writer.writerow(header)
            snp_iterator = self.get_snps_iterator(is_chunked=True,chunk_size=chunk_size)
            chromosomes = self.chromosomes
            positions = self.positions
            num = len(positions)
            log_iteration_count = num/(chunk_size*10)
            for i,snps in enumerate(snp_iterator):
                current_ix = i* chunk_size
                if i % log_iteration_count == 0:
                    log.info('Line %s/%s of Genotype written' % (current_ix,num))
                rows = numpy.hstack((numpy.asarray(zip(chromosomes[current_ix:current_ix+chunk_size],positions[current_ix:current_ix+chunk_size])),snps))
                csv_writer.writerows(rows.tolist())
        log.info('Finished writing genotype ')

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
            bc = scipy.unique(snps)
            if len(bc) == 1:
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
        num_accessions = len(self.accessions)
        # Faster (2.2 ms) than doing numpy.bincount (91.1ms per loop)
        # TODO chunk size is hardcoded
        for i,snps in enumerate(self.get_snps_iterator(is_chunked=True)):
            sm = numpy.sum(snps,axis=1)
            snps_ix.extend(numpy.where( (sm == 0) | (sm == num_accessions))[0]+i*1000)
        numRemoved = len(snps_ix)
        self.filter_snps_ix(snps_ix)
        log.info("Removed %d non-binary SNPs, leaving %d SNPs in total." % (numRemoved, self.num_snps))
        return (num_snps,numRemoved)

    def calculate_ld(self,chr_pos):
        # sort chr_pos first
        chr_pos = sorted(chr_pos,key=itemgetter(0,1))
        indices,snps = self.get_snps_from_pos(chr_pos)
        return calculate_ld(numpy.vstack(snps))


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




    @property
    def accessions(self):
        return self._accessions

    @property
    def snps(self):
        return self._snps

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

    def get_snps_iterator(self,chr=None,is_chunked=False,chunk_size=1000):
        start = 0
        end = None
        if chr is not None:
            chr_region = self.chr_regions[self.get_chr_region_ix(chr)]
            start = chr_region[0]
            end = chr_region[1]
        if is_chunked:
            for i in xrange(start,end,chunk_size):
                stop_i = min(i + chunk_size, end)
                yield self._snps[i:stop_i]
        else:
            for snp in self._snps[start:end]:
                yield snp


    def convert_data_format(self,target_format='binary',reference_ecotype='6909'):
        """
        Converts the underlying raw data format to a binary one, i.e. A,C,G,T,NA,etc. are converted to 0,1,-1
        """
        log.info('Converting from %s to %s' % (self.data_format,target_format))
        if self.data_format == target_format:
            log.warning("Data appears to be already in %s format!" % target_format)
        else:
            if self.data_format == 'nucleotides' and target_format == 'binary':
                self._convert_snps_to_binary(reference_ecotype = reference_ecotype)
            else:
                raise NotImplementedError

    def _convert_snps_to_binary(self,reference_ecotype='6909'):
        missingVal = 'NA'
        decoder = {missingVal:-1} #Might cause errors somewhere???!!!
        coding_fun = scipy.vectorize(lambda x: decoder[x], otypes=['int8'])

        if reference_ecotype in self.accessions:
            ref_i = self.accessions.index(reference_ecotype)
        else:
            ref_i = 0
            log.warning("Given reference ecotype %s wasn't found, using %s as 0-reference." % \
                    (reference_ecotype, self.accessions[ref_i]))
        snps = []
        num = len(self.positions)
        positions = []
        chr_region_removed_snps = [0] * len(self.chr_regions)
        for snp_i, (snp, pos) in enumerate(itertools.izip(self._snps, self.positions)):
            chr_region_ix,chr_region = self.get_chr_region_from_index(snp_i)
            if snp_i % (num / 10) == 0:
                log.info('Converted %s/%s of SNPs.' % (snp_i,num))
            unique_nts = scipy.unique(snp).tolist()
            if missingVal in unique_nts:
                if len(unique_nts) != 3:
                    chr_region_removed_snps[chr_region_ix] +=1
                    continue #Skipping non-binary SNP
                else:
                    unique_nts.remove(self.missingVal)
            else:
                if len(unique_nts) != 2:
                    chr_region_removed_snps[chr_region_ix] +=1
                    continue #Skipping non-binary SNP
            if snp[ref_i] != missingVal:
                ref_nt = snp[ref_i]
                decoder[ref_nt] = 0
                unique_nts.remove(ref_nt)
                decoder[unique_nts[0]] = 1
            else:
                decoder[unique_nts[0]] = 0
                decoder[unique_nts[1]] = 1
            snps.append(coding_fun(snp))
            positions.append(pos)
        log.info('Removed %d non-binary SNPs out of %d, when converting to binary SNPs.'\
            % (len(self.positions) - len(positions), len(self.positions)))
        assert len(snps) == len(positions), 'Somthing odd with the lengths.'
        chr_regions = self.chr_regions[:]
        sum_removed_snps = 0
        for i,num_snps in enumerate(chr_region_removed_snps):
            sum_removed_snps +=num_snps
            if i == 0:
                chr_regions[0] = (0,chr_regions[0][1] - num_snps)
            else:
                chr_regions[i] = (chr_regions[i-1][1],chr_regions[i][1] - sum_removed_snps)
        self._chr_regions = chr_regions

        self._snps = snps
        self._positions = positions
        self._data_format = 'binary'


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
        # TODO update chr_regions
        log.debug("Removed %d accessions, leaving %d in total." % (num_accessions - len(indicesToKeep), len(indicesToKeep)))


class HDF5Genotype(AbstractGenotype):

    def __init__(self,hdf5_file):
        self.h5file = h5py.File(hdf5_file, 'r')
        self.filter_snps = None
        self.accession_filter = None

    def __del__(self):
        if self.h5file is not None:
           self.h5file.close()

    def get_snps(self):
        return self.h5file['snps']

    @property
    def snps(self):
        return self.h5file['snps']

    @property
    def data_format(self):
        return self.h5file['snps'].attrs['data_format']



    @property
    def accessions(self):
        if self.accession_filter is None or len(self.accession_filter) == 0:
            return self.h5file['accessions'][:]
        else:
            return self.h5file['accessions'][self.accession_filter]

    def convert_data_format(self,target_format='binary'):
        raise NotImplementedError

    def _get_snps_(self, start=0,end=None,chunk_size=1000):
        """
        An generator/iterator for SNP chunks.
        """
        if end is None:
            end = self.original_num_snps
        for i in xrange(start,end,chunk_size):
            stop_i = min(i + chunk_size, end)
            if self.accession_filter is None or len(self.accession_filter) == 0:
                if self.filter_snps is None:
                    yield self.h5file['snps'][i:stop_i]
                else:
                    yield self.h5file['snps'][self.filter_snps[i:stop_i],:]
            else:
                if self.filter_snps is None:
                    # first read the entire row and then filter columns (7.91ms vs 94.3ms)
                    yield self.h5file['snps'][i:stop_i][:,self.accession_filter]
                else:
                    filter_chunk = self.filter_snps[i:stop_i]
                    # first read the entire row and then filter columns (7.91ms vs 94.3ms)
                    snps_chunk = self.h5file['snps'][i:stop_i][:,self.accession_filter]
                    yield snps_chunk[filter_chunk]

    def get_snps_iterator(self,chr=None,is_chunked=False,chunk_size=1000):
        """
        Returns an generator containing a chunked generator.
        If chr is passed the generator will only iterate over the specified chr
        """
        start = 0
        end = None
        if chr is not None:
			# use unfiltered chr_regions because filtering happens in _get_snps_
            chr_region = self.h5file['positions'].attrs['chr_regions'][self.get_chr_region_ix(chr)]
            start = chr_region[0]
            end = chr_region[1]
        for snp_chunk in self._get_snps_(start=start,end=end,chunk_size=chunk_size):
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
        if self.filter_snps is not None:
            return self.h5file['positions'][self.filter_snps]
        return self.h5file['positions'][:]




    @property
    def chrs(self):
        return self.h5file['positions'].attrs['chrs']

    @property
    def num_snps(self):
        if self.filter_snps is not None:
            return self.filter_snps.sum()
        return self.original_num_snps

    @property
    def chr_regions(self):
        if self.filter_snps is not None:
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
        self.accession_filter = indicesToKeep
        log.debug("Removed %d accessions, leaving %d in total." % (num_accessions - len(indicesToKeep), len(indicesToKeep)))


    def filter_snps_ix(self,snps_ix):
        if snps_ix is None or len(snps_ix) == 0:
            self.filter_snps = None
            self.filtered_chr_regions = None
        else:
            self.filter_snps = numpy.ones((self.num_snps,),dtype=numpy.bool)
            self.filter_snps[snps_ix] = 0
            self.filtered_chr_regions = self._get_filtered_regons()


    def _get_filtered_regons(self):
        filtered_chr_regions = []
        if self.filter_snps is None:
            return None
        start_ix = 0
        end_ix = 0
        for chr_region in self.h5file['positions'].attrs['chr_regions']:
            end_ix = start_ix + self.filter_snps[chr_region[0]:chr_region[1]].sum()
            filtered_chr_regions.append((start_ix,end_ix))
            start_ix = end_ix
        return filtered_chr_regions





























