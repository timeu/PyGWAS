import itertools as it
import sys
import numpy
import scipy as sp
import logging

log = logging.getLogger(__name__)

SUPPORTED_TRANSFORMATIONS = ("none","log", "sqrt", "exp", "sqr", "arcsin_sqrt", "box_cox","ascombe", "most_normal")

class Phenotype(object):
    """
    A class that encapsulates phenotype values and provides basic functionality for these.

    This is an update of an older class.
    """

    def __init__(self, ecotypes,values,name = None, phenotype_id = None):
        if len(values) != len(ecotypes):
            raise Exception('number of values must match with number of ecotypes')
        self._values  = values
        self._ecotypes = ecotypes
        self._transformation = None
        self._raw_values = None
        self._name = name
        self._phenotype_id = phenotype_id

    @property
    def num_vals(self):
        return len(self.values)

    @property
    def name(self):
        return self._name

    @property
    def transformation(self):
        return self._transformation

    @transformation.setter
    def transformation(self,transformation):
        self._transformation = transformation

    @property
    def phenotype_id(self):
        return self._phenotype_id

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self,values):
        self._values = values

    @property
    def ecotypes(self):
        return self._ecotypes

    @ecotypes.setter
    def ecotypes(self,ecotypes):
        self._ecotypes = ecotypes

    @property
    def raw_values(self):
        return self._raw_values

    @raw_values.setter
    def raw_values(self,raw_values):
        self._raw_values = raw_values


    # move to kinship
    def get_pseudo_heritability(self,K):
        """
        Returns the REML estimate of the heritability.

        methods: 'avg' (averages), 'repl' (replicates)
        """
        from scipy import stats
        import linear_models as lm
        lmm = lm.LinearMixedModel(self.phen_vals)
        if len(self.values) == len(set(self.values)):
            lmm.add_random_effect(K)
        else:
            Z = self.get_incidence_matrix()
            lmm.add_random_effect(Z * K * Z.T)
        r1 = lmm.get_REML()
        ll1 = r1['max_ll']
        rlm = lm.LinearModel(self.values)
        ll0 = rlm.get_ll()
        lrt_stat = 2 * (ll1 - ll0)
        pval = stats.chi2.sf(lrt_stat, 1)
        return {'pseudo_heritability':r1['pseudo_heritability'], 'pval':pval}


    # TODO move to kinship
    def get_blup(self,K):
        """
        Returns the REML estimate for the BLUP and the pseudo-heritability.
        """
        from scipy import stats
        import linear_models as lm
        phen_vals = self.values
        lmm = lm.LinearMixedModel(phen_vals)
        if len(phen_vals) == len(set(phen_vals)):
            lmm.add_random_effect(K)
        else:
            Z = self.get_incidence_matrix()
            lmm.add_random_effect(Z * K * Z.T)
        r1 = lmm.get_REML()
        ll1 = r1['max_ll']
        rlm = lm.LinearModel(phen_vals)
        ll0 = rlm.get_ll()
        lrt_stat = 2 * (ll1 - ll0)
        pval = stats.chi2.sf(lrt_stat, 1)

        #Now the BLUP.
        y_mean = sp.mean(lmm.Y)
        Y = lmm.Y - y_mean
        p_herit = r1['pseudo_heritability']
        delta = (1 - p_herit) / p_herit
#        if K_inverse == None:
#            K_inverse = K.I
#        M = (sp.eye(K.shape[0]) + delta * K_inverse)
#        u_blup = M.I * Y
        M = (K + delta * sp.eye(K.shape[0]))
        u_mean_pred = K * (M.I * Y)
        blup_residuals = Y - u_mean_pred
        return {'pseudo_heritability':r1['pseudo_heritability'], 'pval':pval, 'u_blup':u_mean_pred, 'blup_residuals':blup_residuals}

    def has_replicates(self):
        ets = map(int, self.ecotypes)
        num_vals = len(ets)
        num_ets = len(set(ets))
        return num_vals != num_ets

    # TODO move to linear_models
    def get_broad_sense_heritability(self):
        """
        Estimates the broad sense heritability from replicates.
        """
        import linear_models as lm
        ets = map(int, self.ecotypes)
        num_vals = len(ets)
        num_ets = len(set(ets))
        if num_vals == num_ets:
            log.warning("Can't estimate the broad sense heritability when replicates are missing.")
        else:
            avg_repl_num = float(num_vals) / num_ets
            log.info('Average number of replicates:%s'% avg_repl_num)
        values = self.values
        mod = lm.LinearModel(values)
        res = mod.anova_f_test([sp.array(ets)])
        bs_herit = res['var_perc'][0]
        bs_herit_pval = res['ps'][0]
        log.info('Heritability:%s' % bs_herit)
        log.info('Heritability (different from 0) p-value :%s'% bs_herit_pval)
        bs_avg_herit = 1.0 - (1.0 - bs_herit) / avg_repl_num
        log.info('Estimated mean value heritability:%s' % bs_avg_herit)
        return {'bs_herit':bs_herit, 'bs_pid':bs_pid, 'bs_avg_herit':bs_avg_herit, 'bs_herit_pval':bs_herit_pval}


    def _perform_transform(self,values,transformation):
        if not self.transformation:
            self.raw_values = self.values
            self.transformation = transformation
        else:
            self.transformation = '%s(' + self.transformation + ')' % transformation
        self.values = values.tolist()

    def _log_transform(self, method='standard'):
        a = sp.array(self.values)
        if method == 'standard':
            vals = sp.log((a - min(a)) + 0.1 * sp.var(a))
        else:
            vals = sp.log(a)
        self._perform_transform(vals,"log")
        return True

    def _sqrt_transform(self, method='standard'):
        a = sp.array(self.values)
        if method == 'standard':
            vals = sp.sqrt((a - min(a)) + 0.1 * sp.var(a))
        else:
            vals = sp.sqrt(a)
        self._perform_transform(vals,"sqrt")
        return True


    def _ascombe_transform(self):
        a = sp.array(self.values)
        vals = 2.0 * sp.sqrt(a + 3.0 / 8.0)
        self._perform_transform(vals,"ascombe")
        return True


    def _sqr_transform(self,  method='standard'):
        a = sp.array(self.values)
        if method == 'standard':
            vals = ((a - min(a)) + 0.1 * sp.var(a)) * ((a - min(a)) + 0.1 * sp.var(a))
        else:
            vals = a * a
        self._perform_transform(vals,"sqr")
        return True

    def _exp_transform(self, method='standard'):
        a = sp.array(self.values)
        if method == 'standard':
            vals = sp.exp((a - min(a)) + 0.1 * sp.var(a))
        else:
            vals = sp.exp(a)
        self._perform_transform(vals,"exp")
        return True

    def _arcsin_sqrt_transform(self, verbose=False):
        a = sp.array(self.values)
        if min(a) < 0 or max(a) > 1:
            log.debug('Some values are outside of range [0,1], hence skipping transformation!')
            return False
        else:
            vals = sp.arcsin(sp.sqrt(a))
        self._perform_transform(vals,"arcsin")
        return True

    def _box_cox_transform(self, verbose=False, method='standard'):
        """
        Performs the Box-Cox transformation, over different ranges, picking the optimal one w. respect to normality.
        """
        from scipy import stats
        a = sp.array(self.values)
        if method == 'standard':
            vals = (a - min(a)) + 0.1 * sp.var(a)
        else:
            vals = a
        sw_pvals = []
        lambdas = sp.arange(-2.0, 2.1, 0.1)
        for l in lambdas:
            if l == 0:
                vs = sp.log(vals)
            else:
                vs = ((vals ** l) - 1) / l
            r = stats.shapiro(vs)
            if sp.isfinite(r[0]):
                pval = r[1]
            else:
                pval = 0.0
            sw_pvals.append(pval)
        i = sp.argmax(sw_pvals)
        l = lambdas[i]
        if l == 0:
            vs = sp.log(vals)
        else:
            vs = ((vals ** l) - 1) / l
        self._perform_transform(vs,"box_cox")
        log.debug('optimal lambda was %0.1f' % l)
        return True



    def transform(self, trans_type="most_normal", method='standard', verbose=False):
        log.info('Transforming phenotypes: %s' % trans_type)
        if self.has_replicates():
            log.info('Replicates found. Conveting to averages before applying transformation')
            self.convert_to_averages()
        if trans_type == 'sqrt':
            self._sqrt_transform(method=method)
        elif trans_type == 'ascombe':
            self._ascombe_transform()
        elif trans_type == 'log':
            self._log_transform(method=method)
        elif trans_type == 'sqr':
            self._sqr_transform(method=method)
        elif trans_type == 'exp':
            self._exp_transform(method=method)
        elif trans_type == 'arcsin_sqrt':
            self._arcsin_sqrt_transform()
        elif trans_type == 'box_cox':
            self._box_cox_transform(verbose=verbose)
        elif trans_type == 'most_normal':
            trans_type, shapiro_pval = self.most_normal_transformation(verbose=verbose)
        elif trans_type is None or trans_type == 'none':
            pass
        else:
            raise Exception('Transformation %s unknown' % trans_type)
        return trans_type


    def revert_to_raw_values(self):
        if not self.transformation:
            log.warning('Phenotype values are already raw..')
        else:
            self.transformation = None
            self.values = self.raw_values


    def most_normal_transformation(self,trans_types=SUPPORTED_TRANSFORMATIONS,
                perform_trans=True, verbose=False):
        """
        Performs the transformation which results in most normal looking data, according to Shapiro-Wilk's test
        """
        from scipy import stats
        shapiro_pvals = []
        for trans_type in trans_types:
            if trans_type == 'most_normal':
                continue
            if trans_type != 'none':
                if not self.transform(trans_type=trans_type):
                    continue
            phen_vals = self.values
            #print 'sp.inf in phen_vals:', sp.inf in phen_vals
            if sp.inf in phen_vals:
                pval = 0.0
            else:
                r = stats.shapiro(phen_vals)
                if sp.isfinite(r[0]):
                    pval = r[1]
                else:
                    pval = 0.0
            shapiro_pvals.append(pval)
            if trans_type != 'none':
                self.revert_to_raw_values()
        argmin_i = sp.argmax(shapiro_pvals)
        trans_type = trans_types[argmin_i]
        shapiro_pval = shapiro_pvals[argmin_i]
        if perform_trans:
            self.transform(trans_type=trans_type)
        log.info("The most normal-looking transformation was %s, with a Shapiro-Wilk's p-value of %.2E" % \
                (trans_type, shapiro_pval))
        return trans_type, shapiro_pval


    def normalize_values(self):
        a = sp.array(self.values)
        v = sp.var(self.get_avg_values(), ddof=1)
        vals = a / sp.sqrt(v)
        self.values = vals.tolist()


    def na_outliers(self, iqr_threshold):
        raise NotImplementedError


    def filter_ecotypes(self, indices_to_keep, random_fraction=1):
        """
        Removes the ecotypes from all data.
        """
        import random
        el = []
        vl = []
        if self.transformation:
            rvl = []
        if random_fraction < 1:
            indices = range(len(self.ecotypes))
            indices_to_keep = sorted(random.sample(indices, int(len(self.ecotypes) * random_fraction)))
        for i in indices_to_keep:
            el.append(self.ecotypes[i])
            vl.append(self.values[i])
            if self.transformation:
                rvl.append(self.raw_values[i])
        self.ecotypes = el
        self.values = vl
        if self.transformation:
            self.raw_values = rvl

    def filter_ecotypes_2(self, ecotypes_to_keep, pids=None):
        unique_ets = set()
        el = []
        vl = []
        if self.transformation:
            rvl = []
        for et in ecotypes_to_keep:
            if et in self.ecotypes:
                i = self.ecotypes.index(et)
                el.append(self.ecotypes[i])
                vl.append(self.values[i])
                if self.transformation:
                    rvl.append(self.raw_values[i])
                unique_ets.add(et)
        self.ecotypes = el
        self.values = vl
        if self.transformation:
            self.raw_values = rvl
        return list(unique_ets)


    def order_ecotypes(self, ets_map):
        ets = []
        vals = []
        if self.transformation:
            rvl = []
        for i in ets_map:
            ets.append(self.ecotypes[i])
            vals.append(self.values[i])
            if self.transformation:
                rvl.append(self.raw_values[i])
        self.ecotypes = ets
        self.values = vals
        if self.transformation:
            self.raw_values = rvl





    def get_incidence_matrix(self):
        ets = sp.array(self.ecotypes)
        unique_ets = []
        i = 0
        while i < len(ets):
            et = ets[i]
            unique_ets.append(et)
            while i < len(ets) and ets[i] == et: #The ecotypes are assumed to be sorted
                i += 1
        Z = sp.int8(sp.mat(ets).T == sp.mat(unique_ets))
        return Z


    def _get_ecotype_value_dict(self):
        d = {}
        for et in set(self.ecotypes):
            d[et] = {'values':[], 'rep_num':0}

        for et, val in it.izip(self.ecotypes, self.values):
            d[et]['values'].append(val)
            d[et]['rep_num'] += 1
        return d



    def get_avg_value(self):
        """
        Returns the average values, along with the ecotypes, and rep_number
        """
        d = self._get_ecotype_value_dict()
        ecotypes = d.keys()
        avg_vals = []
        rep_nums = []
        for et in d:
            avg_vals.append(sp.mean(d[et]['values']))
            rep_nums.append(d[et]['rep_num'])
        return {'ecotypes':ecotypes, 'values':avg_vals, 'rep_nums': rep_nums}



    def convert_to_averages(self):
        """
        Replaces phenotypes with replicates with their averages.
        """
        avg_values = self.get_avg_value()
        self.ecotypes = avg_values['ecotypes']
        self.values = avg_values['values']



    def write_to_file(self, file_name, delim=','):
        """
        Writes the object to a file.. (in the new format)
        """
        import csv
        with open(file_name,'w') as f:
            writer = csv.writer(f,delimiter=delim)
            writer.writerow(['ecotype_id', 'value', 'replicate_id'])
            ets_vals = zip(self.ecotypes, self.values)
            ets_vals.sort()
            last_et = -1
            for (et, val) in ets_vals:
                if et != last_et:
                    repl_id = 1
                else:
                    repl_id += 1
                writer.writerow([et,val,repl_id])


    def is_binary(self):
        return len(sp.unique(self.values)) == 2

    def is_constant(self):
        return len(sp.unique(self.values)) == 1

    def is_near_constant(self, min_num_diff=10):
        vals = sp.array(self.values)
        if sp.std(vals) > 0:
            vals = 50 * (vals - sp.mean(vals)) / sp.std(vals)
            vals = vals - vals.min() + 0.1
            b_counts = sp.bincount(sp.array(sp.around(vals), dtype='int'))
            b = b_counts.max() > len(vals) - min_num_diff
            return b
        else:
            return True

def parse_phenotype_file(file_name, delim=','):
    import csv
    """
    Parses a phenotype file, and returns a new phenotype_data object.

    """
    phen_data = {}
    with open(file_name, 'rU') as f:
        reader = csv.reader(f,delimiter=delim)
        ets = []
        values = []
        header = reader.next()
        name = header[1]
        for row in reader:
            try:
                values.append(float(row[1]))
                ets.append(row[0].strip())
            except ValueError:
                log.warning('could not parse row %s' %row)
    return Phenotype(ets,values,name)


















