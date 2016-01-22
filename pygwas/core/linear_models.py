"""
Contains functions to perform various linear regression schemes, such as simple linear regression, and mixed models.

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com

"""

import scipy as sp
from scipy import linalg
from scipy import stats
from scipy import optimize
import genotype
import pdb
import cPickle
import os
import sys
import time
import itertools as it
import math
import os
import kinship
import logging
import warnings

log = logging.getLogger(__name__)




def cholesky(V):
    """
    A ad hoc cholesky decomposition that attempts to use R through rpy2, 
    if the matrix V is singular. 
    """
    H_sqrt = sp.mat(linalg.cholesky(V))
    return H_sqrt


def qr_decomp(X):
    """
    QR decomposition.  Adaptations to changes between scipy versions.
    """
    ver_list = tuple(map(int, (sp.__version__).split('.')[:2]))
    if ver_list >= (0, 9):
        return linalg.qr(X, mode='economic')  # Do the QR-decomposition for the Gram-Schmidt process.            
    else:
        return linalg.qr(X, econ=True)  # Do the QR-decomposition for the Gram-Schmidt process.
    



class LinearModel(object):
    """
    A simple linear model
    """
    def __init__(self, Y=None):
        """
        The fixed effects should be a list of fixed effect lists (SNPs)
        """
        self.n = len(Y)
        self.Y = sp.matrix(Y).T
        self.X = sp.ones((self.n, 1))  # The intercept
        self.p = 1
        self.beta_est = None
        self.cofactors = []



    def add_factor(self, x, lin_depend_thres=1e-8):
        """
        Adds an explanatory variable to the X matrix.
        """
        # Checking whether this new cofactor in linearly independent.
        new_x = sp.array(x)
        new_x.shape = len(x)
        (beta, rss, rank, sigma) = linalg.lstsq(self.X, new_x)
        if float(rss) < lin_depend_thres:
            log.warning('A factor was found to be linearly dependent on the factors already in the X matrix.  Hence skipping it!')
            return False
        new_x.shape = (self.n, 1)
        self.X = sp.hstack([self.X, new_x])
        self.cofactors.append(x)
        self.p += 1
        return True


    def set_factors(self, factors, include_intercept=True):
        """
        Set the cofactors.
        """
        self.p = 0
        if include_intercept:
            self.X = sp.matrix(sp.ones((self.n, 1), dtype='single'))  # The intercept
            self.p = 1
            if len(factors) > 0:
                self.X = sp.hstack([self.X, sp.matrix(factors, dtype='single').T])
            self.p += len(factors)

        else:
            self.X = sp.matrix(factors, dtype='single').T
            self.p = len(factors)



    def get_hat_matrix(self):
        self.X_squared_inverse = (self.X.T * self.X).I
        self.hat_matrix = self.X * self.X_squared_inverse * self.X.T
        return self.hat_matrix



    def least_square_estimate(self):
        """
        Via Normal equations, get LSEs
        """
        self.X_squared_inverse = (self.X.T * self.X).I
        self.beta_est = self.X_squared_inverse * self.X.T * self.Y
        return self.beta_est

    def get_residuals(self):
        """
        Calculates and returns the residuals as a column vector.
        """
        # if self.beta_est==None:
        self.least_square_estimate()
        self.residuals = self.Y - self.X * self.beta_est
        return self.residuals


    def general_f_test(self, A, c):
        """
        A general F-test implementation.
        Where the hypothesis is A*beta=c, a constraint.
        
        Here A is a matrix, and c a column vector        
        """
        # if not self.residuals:
        self.get_residuals()
        q, p = A.shape
        assert p == self.p, 'Shape of A is wrong!'
        B = (A * self.beta_est - c)
        M = A * (self.X_squared_inverse) * A.T
        f_stat = (B.T * M.I * B) / ((self.residuals.T * self.residuals) / (self.n - self.p))
        p_value = 1 - stats.f.cdf(f_stat, q, self.n - self.p)
        return p_value, f_stat



    def get_rss(self, dtype='single'):
        """
        Returns the RSS
        """
        h0_X = sp.mat(self.X, dtype=dtype)
        (betas, rss, r, s) = linalg.lstsq(h0_X, self.Y)
        return rss


    def get_ll(self, rss=None, dtype='single'):
        """
        Returns the log-likelihood
        """
        if not rss:
            rss = self.get_rss(dtype)
        return (-self.n / 2) * (1 + sp.log(2 * sp.pi) + rss / self.n)


    def fast_f_test(self, genotype, verbose=True, Z=None, with_betas=False, progress_file_writer=None):
        """
        LM implementation 
        Single SNPs
                    
        """
        snps = genotype.get_snps_iterator(is_chunked=True)
        dtype = 'single'
        q = 1  # Single SNP is being tested
        p = len(self.X.T) + q
        n = self.n
        n_p = n - p
        num_snps = genotype.num_snps


        log.info('Starting regression scan',extra={'progress':25})
        if not progress_file_writer == None:
            progress_file_writer.update_progress_bar(progress=0.25, task_status='Starting regression scan')
        h0_X = sp.mat(self.X, dtype=dtype)
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, self.Y)
        Y = sp.mat(self.Y - h0_X * h0_betas, dtype=dtype)
        h0_betas = map(float, list(h0_betas))

        if not with_betas:
            (Q, R) = qr_decomp(h0_X)  # Do the QR-decomposition for the Gram-Schmidt process.
            Q = sp.mat(Q, dtype=dtype)
            Q2 = Q * Q.T
            M = sp.mat((sp.eye(n) - Q2), dtype=dtype)
        else:
            betas_list = [h0_betas] * num_snps

        rss_list = sp.repeat(h0_rss, num_snps)
        chunk_size = len(Y)
        i = 0
        for snps_chunk in snps:
            if len(snps_chunk) == 0:
                continue
            snps_chunk = sp.matrix(snps_chunk)
            if with_betas:
                Xs = snps_chunk
            else:
                Xs = sp.mat(snps_chunk, dtype=dtype) * M
            for j in range(len(Xs)):
                if with_betas:
                    (betas, rss, p, sigma) = linalg.lstsq(sp.hstack([h0_X, Xs[j].T]), Y, \
                                    overwrite_a=True)
                    if not rss:
                        log.debug('No predictability in the marker, moving on...')
                        continue
                    betas_list[i] = map(float, list(betas))
                else:
                    (betas, rss, p, sigma) = linalg.lstsq(Xs[j].T, Y, overwrite_a=True)
                rss_list[i] = rss[0]
                
                if (i+1) % (num_snps / 10) == 0:
                    perc = 100.0 * i / num_snps
                    log.info('Performing regression (completed:%d %%)' % perc,extra={'progress':25 + 80*perc/100})
                if not progress_file_writer == None:
                   progress_file_writer.update_progress_bar(task_status='Performing regression (completed: %d %%)' % (100.0 * i / num_snps))
                i += 1

        rss_ratio = h0_rss / rss_list
        var_perc = 1 - 1 / rss_ratio
        f_stats = (rss_ratio - 1) * n_p / float(q)
        p_vals = stats.f.sf(f_stats, q, n_p)

        res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
            'h0_rss':h0_rss, 'h0_betas':h0_betas}
        if with_betas:
            res_d['betas'] = betas_list
        return res_d


    def two_snps_ftest(self, snps, verbose=True):
        """
        Test a pair of SNPs.
        """
        num_snps = len(snps)

        ftest_res = self.fast_f_test(snps)
        full_rss = ftest_res['h0_rss']
        h0_X = self.X
        Y = self.Y  # The transformed outputs.

        # Contructing result containers        
        p3_f_stat_array = sp.zeros((num_snps, num_snps))
        p3_p_val_array = sp.ones((num_snps, num_snps))
        p4_f_stat_array = sp.zeros((num_snps, num_snps))
        p4_p_val_array = sp.ones((num_snps, num_snps))
        f_stat_array = sp.zeros((num_snps, num_snps))
        p_val_array = sp.ones((num_snps, num_snps))
        rss_array = sp.repeat(full_rss, num_snps * num_snps)
        rss_array.shape = (num_snps, num_snps)
        var_perc_array = sp.zeros((num_snps, num_snps))
        haplotype_counts = [[{} for j in range(i + 1)] for i in range(num_snps)]

        # Fill the diagonals with the single SNP emmax
        for i, snp in enumerate(snps):
            hap_set, hap_counts, haplotypes = genotype.get_haplotypes([snp], self.n,
                                        count_haplotypes=True)
            d = {'num_haps':hap_counts}
            for hap, hap_c in zip(hap_set, hap_counts):
                d[hap] = hap_c
            haplotype_counts[i][i] = d
            p_val_array[i, i] = ftest_res['ps'][i]
            p3_p_val_array[i, i] = p_val_array[i, i]
            p4_p_val_array[i, i] = p_val_array[i, i]
            f_stat_array[i, i] = ftest_res['f_stats'][i]
            p3_f_stat_array[i, i] = f_stat_array[i, i]
            p4_f_stat_array[i, i] = f_stat_array[i, i]
            rss_array[i, i] = ftest_res['rss'][i]
            var_perc_array[i, i] = ftest_res['var_perc'][i]


        identical_snp_count = 0
        no_interaction_count = 0
        for i, snp1 in enumerate(snps):
            snp1 = snps[i]
            for j in range(i):
                snp2 = snps[j]
                if i == j: continue  # Skip diagonals..

                # Haplotype counts 
                hap_set, hap_counts, haplotypes = genotype.get_haplotypes([snp1, snp2], self.n,
                                            count_haplotypes=True)
                groups = set(haplotypes)
                d = {'num_haps':len(hap_counts)}
                for hap, hap_c in zip(hap_set, hap_counts):
                    d[hap] = hap_c
                haplotype_counts[i][j] = d

                # Fill the upper right part with more powerful of two tests.

                if ftest_res['ps'][i] < ftest_res['ps'][j]:
                    rss_array[j, i] = ftest_res['rss'][i]
                    max_i = i
                else:
                    rss_array[j, i] = ftest_res['rss'][j]
                    max_i = j

                if d['num_haps'] == 2:
                    identical_snp_count += 1
                    continue
                elif d['num_haps'] == 3:
                    n_p = self.n - 3
                    no_interaction_count += 1
                    # Do ANOVA
                    l = []
                    for g in groups:
                        l.append(sp.int8(haplotypes == g))
                    X = sp.mat(l)
                    (betas, rss, p, sigma) = linalg.lstsq(X.T, Y)
                    rss_array[i, j] = rss[0]
                    var_perc_array[i, j] = 1 - rss / full_rss
                    f_stat = (rss_array[j, i] / rss - 1) * n_p  # FINISH
                    p_val = stats.f.sf([f_stat], 1, n_p)[0]
                    p3_f_stat_array[j, i] = f_stat
                    p3_f_stat_array[i, j] = f_stat
                    p3_p_val_array[j, i] = p_val
                    p3_p_val_array[i, j] = p_val

                    f_stat = ((full_rss - rss) / 2) / (rss / n_p)
                    f_stat_array[j, i] = f_stat
                    p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]


                elif d['num_haps'] == 4:  # With interaction
                    n_p = self.n - 3
                    # Do additive model
                    snp_mat = sp.mat([snp1, snp2]).T  # Transformed inputs
                    X = sp.hstack([h0_X, snp_mat])
                    (betas, rss, rank, s) = linalg.lstsq(X, Y)
                    f_stat = (rss_array[j, i] / rss - 1) * n_p  # Compared to only one SNP
                    p_val = stats.f.sf([f_stat], 1, n_p)[0]
                    rss_array[i, j] = rss
                    p3_f_stat_array[j, i] = f_stat
                    p3_p_val_array[j, i] = p_val

#                    v_f_stat_array[j, i] = f_stat
#                    v_p_val_array[j, i] = stats.f.sf([f_stat], 1, n_p)[0]

                    f_stat = ((full_rss - rss) / 2) / (rss / n_p)  # Compared to only the intercept
                    f_stat_array[j, i] = f_stat
                    p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]

                    # Generate the interaction, and test it.
                    i_snp = snp1 & snp2
                    snp_mat = sp.mat([i_snp]).T
                    X = sp.hstack([h0_X, sp.mat([snps[max_i]]).T, snp_mat])
                    (betas, rss, rank, s) = linalg.lstsq(X, Y)
                    f_stat = (rss_array[j, i] / rss - 1) * n_p  # Compared to only one SNP
                    p_val = stats.f.sf([f_stat], 1, n_p)[0]
                    p3_f_stat_array[i, j] = f_stat
                    p3_p_val_array[i, j] = p_val


                    # full model p-value
                    n_p = self.n - 4
                    l = []
                    for g in groups:
                        l.append(sp.int8(haplotypes == g))
                    X = sp.mat(l)
                    (betas, rss, p, sigma) = linalg.lstsq(X.T, Y)

                    f_stat = ((rss_array[j, i] - rss) / 2) / (rss / n_p)  # Compared to only one SNP
                    p_val = stats.f.sf([f_stat], 2, n_p)[0]
                    p4_f_stat_array[j, i] = f_stat
                    p4_p_val_array[j, i] = p_val

                    f_stat = (rss_array[i, j] / rss - 1) * n_p  # Compared to two SNPs
                    p4_f_stat_array[i, j] = f_stat
                    p4_p_val_array[i, j] = stats.f.sf([f_stat], 1, n_p)[0]

                    f_stat = ((full_rss - rss) / 3) / (rss / n_p)  # Compared to only intercept
                    f_stat_array[i, j] = f_stat
                    p_val_array[i, j] = stats.f.sf([f_stat], 3, n_p)[0]
                    rss_array[j, i] = rss

            if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()

        print no_interaction_count, identical_snp_count

        #FINISH res dict!!!
        res_dict = {'p3_ps':p3_p_val_array, 'p3_f_stats':p3_f_stat_array, 'p4_ps':p4_p_val_array,
            'p4_f_stats':p4_f_stat_array, 'rss':rss_array, 'var_perc':var_perc_array,
            'haplotype_counts':haplotype_counts,
            'f_stats':f_stat_array, 'ps':p_val_array}
        return res_dict



    def anova_f_test(self, snps, dtype='int8'):
        """
        A standard ANOVA, using a F-test
        """
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X, self.Y)
        num_snps = len(snps)
        rss_list = sp.repeat(h0_rss, num_snps)
        h0_betas = map(float, list(h0_betas)) + [0.0]
        betas_list = [h0_betas] * num_snps
        var_perc = sp.zeros(num_snps)
        f_stats = sp.zeros(num_snps)
        dfs = sp.zeros(num_snps)
        p_vals = sp.ones(num_snps)
        n = self.n
        p_0 = len(self.X.T)

        for i, snp in enumerate(snps):
            groups = sp.unique(snp)
            q = len(groups) - 1  # Null model has 1 df.
            p = p_0 + q
            n_p = n - p
            x = sp.zeros((len(groups), n), dtype=dtype)
            for g_i, g in enumerate(groups):
                x[g_i] = sp.int32(snp == g)
            (betas, rss, p, sigma) = linalg.lstsq(sp.mat(x).T, self.Y)

            if not rss:
                print 'No predictability in the marker, moving on...'
                continue
            rss_list[i] = rss[0]
            betas_list[i] = map(float, list(betas))
            rss_ratio = h0_rss / rss
            var_perc[i] = 1 - 1 / rss_ratio
            f_stat = (rss_ratio - 1) * n_p / float(q)
            p_vals[i] = stats.f.sf([f_stat], q, n_p)[0]
            f_stats[i] = f_stat
            dfs[i] = n_p
            if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()

        if num_snps >= 10:
            sys.stdout.write('\n')

        return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
            'var_perc':var_perc, 'dfs':dfs}



    def anova_f_test_w_missing(self, snps):
        """
        A standard ANOVA, using a F-test
        
        Handles SNPs w. missing data...
        
        Warning: This is a very time consuming approach. A faster, and probably a better, 
                 approach would be to impute the missing values a priori 

        """
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(self.X, self.Y)
        num_snps = len(snps)
        rss_list = sp.repeat(h0_rss, num_snps)
        h0_betas = map(float, list(h0_betas)) + [0.0]
        betas_list = [h0_betas] * num_snps
        var_perc = sp.zeros(num_snps)
        f_stats = sp.zeros(num_snps)
        dfs = sp.zeros(num_snps)
        p_vals = sp.ones(num_snps)
        n = self.n
        p_0 = len(self.X.T)

        for i, snp in enumerate(snps):
            groups = sp.unique(snp)
            q = len(groups) - 1  # Null model has 1 df.
            p = p_0 + q
            n_p = n - p
            x = []
            for g in groups:
                x.append(sp.int8(snp == g))
            (betas, rss, p, sigma) = linalg.lstsq(sp.mat(x).T, self.Y)

            if not rss:
                log.debug('No predictability in the marker, moving on...')
                continue
            rss_list[i] = rss[0]
            betas_list[i] = map(float, list(betas))
            rss_ratio = h0_rss / rss
            var_perc[i] = 1 - 1 / rss_ratio
            f_stat = (rss_ratio - 1) * n_p / float(q)
            p_vals[i] = stats.f.sf([f_stat], q, n_p)[0]
            f_stats[i] = f_stat
            dfs[i] = n_p
            if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()

        if num_snps >= 10:
            sys.stdout.write('\n')

        return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
            'var_perc':var_perc, 'dfs':dfs}





    def test_explanatory_variable(self, x):
        """
        Returns a p-value for whether adding this variable to the model explains the model better.
        
        Hopefully a sped-up version of a specific F-test.
        """

        # THIS CAN BE SPED-UP MORE IF WE CHECK WHETHER self.X IS A VECTOR.  
        # AND USE t-TEST. 
        res_1 = self.get_residuals()

        X_mat = sp.hstack([self.X, sp.matrix([[v] for v in x])])
        X_sq = X_mat.T * X_mat
        try:
            X_sq_inv = X_sq.I
        except Exception, err_str:
            log.error(err_str)
            raise Exception('Andskotinn!!')

        res_2 = self.Y - X_mat * X_sq_inv * X_mat.T * self.Y
        rss_1 = res_1.T * res_1
        rss_2 = res_2.T * res_2
        f_stat = (rss_1 - rss_2) / (rss_2 / (self.n - self.p + 1))
        p_value = stats.f.sf(f_stat, 1, self.n - self.p + 1)
        return p_value, f_stat






class LinearMixedModel(LinearModel):
    """
    A class for linear mixed models
    """
    def __init__(self, Y=None, dtype='single'):
        """
        The fixed effects should be a list of fixed effect lists (SNPs)
        """
        self.n = len(Y)
        self.y_var = sp.var(Y, ddof=1, dtype=dtype)
        self.Y = sp.matrix(Y, dtype=dtype)
        self.Y.shape = (self.n, 1) 
        self.X = sp.matrix(sp.ones((self.n, 1), dtype=dtype))  # The intercept
        self.p = 1
        self.beta_est = None
        self.cofactors = []

        # A list of random effect type, and the cov matrix.
        self.random_effects = [('normal', sp.matrix(sp.identity(self.n)))]  # The first random effect is the IID error.


    def add_random_effect(self, cov_matrix=None, effect_type='normal'):
        if effect_type != 'normal':
            raise Exception('Currently, only Normal random effects are allowed.')
        self.random_effects.append((effect_type, kinship.scale_k(cov_matrix)))


    def set_random_effect(self, cov_matrix_list, effect_types=None):
        self.random_effects = [('normal', sp.matrix(sp.identity(self.n)))]
        for cov_matrix in cov_matrix_list:
            self.add_random_effect(cov_matrix=kinship.scale_k(cov_matrix))


    def _get_eigen_L_(self, K=None, dtype='single'):
        if K is None:
            K = self.random_effects[1][1]
        if sp.__version__ < '0.8':
            K = sp.mat(K, dtype=dtype)
        evals, evecs = linalg.eigh(K)
        evals = sp.array(evals, dtype=dtype)
        return {'values':evals, 'vectors':sp.mat(evecs, dtype=dtype).T}



    def _get_eigen_R_(self, X=None, K=None, hat_matrix=None, dtype='single'):
        if X is None:
            X = self.X
        q = X.shape[1]
        if not hat_matrix:
            X_squared_inverse = linalg.pinv(X.T * X)  # (X.T*X).I
            hat_matrix = X * X_squared_inverse * X.T
        if K is None:
            K = self.random_effects[1][1]
        S = sp.mat(sp.identity(self.n)) - hat_matrix  # S=I-X(X'X)^{-1}X'
        M = sp.mat(S * (K + self.random_effects[0][1]) * S, dtype='double')
        if sp.__version__ < '0.8':
            M = sp.mat(M, dtype=dtype)
        evals, evecs = linalg.eigh(M, overwrite_a=True)  # eigen of S(K+I)S
        eig_values = sp.array(evals[q:], dtype=dtype) - 1  # Because of S(K+I)S?
        return {'values':eig_values, 'vectors':(sp.mat(evecs, dtype=dtype).T[q:])}


    def _rell_(self, delta, eig_vals, sq_etas):
        num_eig_vals = len(eig_vals)
        c_1 = 0.5 * num_eig_vals * (sp.log(num_eig_vals / (2.0 * sp.pi)) - 1)
        v = eig_vals + delta
        res = c_1 - 0.5 * (num_eig_vals * sp.log(sp.sum(sq_etas.flatten() / v)) + sp.sum(sp.log(v)))
        return res  # log-likelihoods (eq. 7 from paper)


    def _redll_(self, delta, eig_vals, sq_etas):
        num_eig_vals = len(eig_vals)
        v1 = eig_vals + delta
        v2 = sq_etas.flatten() / v1
        res = (num_eig_vals * sp.sum(v2 / v1) / sp.sum(v2) - sp.sum(1.0 / v1))
        return res  # diffrentiated log-likelihoods (*2) (eq. 9 from paper)


    def _ll_(self, delta, eig_vals, eig_vals_L, sq_etas):
        n = self.n
        c_1 = 0.5 * n * (sp.log(n / (2.0 * sp.pi)) - 1)
        v1 = eig_vals + delta
        v2 = eig_vals_L + delta
        res = c_1 - 0.5 * (n * sp.log(sp.sum(sq_etas.flatten() / v1)) + sp.sum(sp.log(v2)))
        return res  # log-likelihoods (eq. 6 from paper)


    def _dll_(self, delta, eig_vals, eig_vals_L, sq_etas):
        num_eig_vals = len(eig_vals)
        v1 = eig_vals + delta
        v2 = sq_etas.flatten() / v1
        v3 = eig_vals_L + delta
        res = (self.n * sp.sum(v2 / v1) / sp.sum(v2) - sp.sum(1.0 / v3))
        return res  # diffrentiated log-likelihoods (*2) (eq. 8 from paper)



    def get_REML(self, ngrids=100, llim= -10, ulim=10, esp=1e-6, eig_L=None, eig_R=None):
        """
        Get REML estimates for the effect sizes, as well as the random effect contributions.
        
        This is EMMA
        """
        if not eig_L:
            K = self.random_effects[1][1]
            eig_L = self._get_eigen_L_(K)

        # Get the variance estimates..
        res = self.get_estimates(eig_L, ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='REML',
                    eig_R=eig_R)
        # res = self.get_estimates(ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='REML', K=K)
        res['eig_L'] = eig_L
        return res



    def get_ML(self, ngrids=100, llim= -10, ulim=10, esp=1e-6, eig_L=None, eig_R=None, H=None, H_inv=None, H_sqrt_inv=None, dtype='single'):
        """
        Get ML estimates for the effect sizes, as well as the random effect contributions.
        
        H is the (full) covariance matrix, which speeds up calculations if it's available.
        """
        if H is None:
            if not eig_L:
                K = self.random_effects[1][1]
                eig_L = self._get_eigen_L_(K)
            res = self.get_estimates(eig_L, ngrids=ngrids, llim=llim, ulim=ulim, esp=esp, method='ML',
                    eig_R=eig_R)
        else:
            # If the variance matrix is given, then we calculate the likelihood using that.
            evals, evecs = linalg.eigh(H)
            X_t = sp.mat(H_sqrt_inv * self.X, dtype=dtype)
            Y_t = H_sqrt_inv * self.Y  # The transformed outputs.
            (betas, mahalanobis_rss, rank, hs) = linalg.lstsq(X_t, Y_t)
            rss = sp.sum(sp.array(self.Y - sp.dot(self.X, betas)) ** 2)
            betas = map(float, list(betas))
            n = Y_t.shape[0]
            ll = -0.5 * (n * sp.log(2 * sp.pi) + sp.sum(sp.log(evals)) + mahalanobis_rss)
            assert len(mahalanobis_rss) > 0, 'WTF?'
            res = {'ll': ll, 'rss':rss, 'mahalanobis_rss':mahalanobis_rss} 
        return res


    def get_estimates_3(self, xs=None, ngrids=[5, 5, 5, 5, 5], llim= -3, ulim=3, method='REML',
                verbose=False, dtype='single', debug=False):
        """
        Handles two K matrices, and one I matrix.
        
        Methods available are 'REML', and 'ML'        
        """

        if verbose:
            print 'Retrieving %s variance estimates' % method
        if xs is not None:
            X = sp.hstack([self.X, xs])
        else:
            X = self.X
#        xs1 = [[] for i in range(len(ngrids))]
#        ys1 = [[] for i in range(len(ngrids))]
#        xs2 = []
#        ys2 = []
        for it_i in range(len(ngrids)):
            delta = float(ulim - llim) / ngrids[it_i]
            # print llim, ulim
            # print delta
            log_k_ratio = llim
            lls = []
            res_list = []
            for i in range(ngrids[it_i] + 1):
                # xs1[it_i].append(log_k_ratio)
                # xs2.append(log_k_ratio)
                k_ratio = sp.exp(log_k_ratio)
                a = k_ratio / (k_ratio + 1.0)
                K = a * self.random_effects[1][1] + (1 - a) * self.random_effects[2][1]
                eig_L = self._get_eigen_L_(K=K)
                # Now perform EMMA
                res_dict = self.get_estimates(eig_L, K=K, xs=xs, ngrids=10, method=method, llim= -8, ulim=8)
                res_list.append(res_dict)
                lls.append(res_dict['max_ll'])
                # ys1[it_i].append(res_dict['max_ll'])
                # ys2.append(res_dict['pseudo_heritability'])
                log_k_ratio += delta
            max_ll_i = sp.argmax(lls)
            # print 'max_ll_i:', max_ll_i
            # print 'delta:', delta
            # Update the ulim and llim
            ulim = llim + delta * (max_ll_i + 1)
            llim = llim + delta * (max_ll_i - 1)
#        opt_log_k_ratio = llim + delta * max_ll_i
        res_dict = res_list[max_ll_i]
        opt_k_ratio = sp.exp(log_k_ratio)
        a = opt_k_ratio / (opt_k_ratio + 1)
        opt_k = kinship.scale_k(a * self.random_effects[1][1] + (1 - a) * self.random_effects[2][1])
        res_dict = self.get_estimates(eig_L, K=opt_k, xs=xs, ngrids=20, method=method, llim= -8, ulim=8)
        res_dict['opt_k'] = opt_k
#        res_dict['opt_k_ratio'] = opt_k_ratio
        res_dict['perc_var1'] = a * res_dict['pseudo_heritability']
        res_dict['perc_var2'] = (1 - a) * res_dict['pseudo_heritability']
#        if debug:
#            import pylab
#            pylab.figure()
#            for i in range(3):
#                pylab.plot(xs1[i], ys1[i], marker='o', ls='None', alpha=0.4)
#            pylab.xlabel('log(K ratio)')
#            pylab.ylabel('log-likelihood')
#            pylab.savefig('/Users/bjarni.vilhjalmsson/tmp/%s.png' % str(debug))
#            pylab.figure()
#            pylab.plot(xs2, ys2, marker='o', ls='None', alpha=0.5)
#            pylab.xlabel('log(K ratio)')
#            pylab.ylabel('pseudo_heritability')
#            pylab.savefig('/Users/bjarni.vilhjalmsson/tmp/%s.png' % str(debug))
        return res_dict



    def get_estimates(self, eig_L, K=None, xs=None, ngrids=50, llim= -10, ulim=10, esp=1e-6,
                return_pvalue=False, return_f_stat=False, method='REML', verbose=False,
                dtype='single', eig_R=None, rss_0=None):
        """
        Get ML/REML estimates for the effect sizes, as well as the random effect contributions.
        Using the EMMA algorithm (Kang et al., Genetics, 2008).
        
        Methods available are 'REML', and 'ML'        
        """
        if verbose:
            print 'Retrieving %s variance estimates' % method
        if xs is not None:
            X = sp.hstack([self.X, xs])
        else:
            X = self.X

        if not (eig_R and xs is not None):
            eig_R = self._get_eigen_R_(X=X, K=K)
        q = X.shape[1]  # number of columns
        n = self.n
        p = n - q
        m = ngrids + 1

        etas = sp.array(eig_R['vectors'] * self.Y, dtype=dtype)
        sq_etas = etas * etas
        log_deltas = (sp.arange(m, dtype=dtype) / ngrids) * (ulim - llim) + llim  # a list of deltas to search
        deltas = sp.exp(log_deltas)
        assert len(deltas) == m, 'Number of deltas is incorrect.'
        eig_vals = sp.array(eig_R['values'], dtype=dtype)
        assert len(eig_vals) == p, 'Number of eigenvalues is incorrect.'

        lambdas = sp.reshape(sp.repeat(eig_vals, m), (p, m)) + sp.reshape(sp.repeat(deltas, p), (m, p)).T
        s1 = sp.sum(sq_etas / lambdas, axis=0)
        if method == 'REML':
            if verbose: print 'Calculating grid REML values'
            s2 = sp.sum(sp.log(lambdas), axis=0)
            lls = 0.5 * (p * (sp.log((p) / (2.0 * sp.pi)) - 1 - sp.log(s1)) - s2)
            s3 = sp.sum(sq_etas / (lambdas * lambdas), axis=0)
            s4 = sp.sum(1 / lambdas, axis=0)
            dlls = 0.5 * (p * s3 / s1 - s4)
        elif method == 'ML':
            if verbose: print 'Calculating grid ML values'
            # Xis < -matrix(eig.L$values, n, m) + matrix(delta, n, m, byrow=TRUE)
            eig_vals_L = sp.array(eig_L['values'], dtype=dtype)
            xis = sp.reshape(sp.repeat(eig_vals_L, m), (n, m)) + \
                sp.reshape(sp.repeat(deltas, n), (m, n)).T
            # LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))    
            # dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))    

            s2 = sp.sum(sp.log(xis), axis=0)
            lls = 0.5 * (n * (sp.log((n) / (2.0 * sp.pi)) - 1 - sp.log(s1)) - s2)
            s3 = sp.sum(sq_etas / (lambdas * lambdas), axis=0)
            s4 = sp.sum(1 / xis, axis=0)
            dlls = 0.5 * (n * s3 / s1 - s4)

        max_ll_i = sp.argmax(lls)
        max_ll = lls[max_ll_i]

        last_dll = dlls[0]
        last_ll = lls[0]
        zero_intervals = []
        for i in range(1, len(dlls)):
            if dlls[i] < 0 and last_dll > 0:
                zero_intervals.append(((lls[i] + last_ll) * 0.5, i))
            last_ll = lls[i]
            last_dll = dlls[i]

        if len(zero_intervals) > 0:
            if verbose: print 'Found a zero interval... now performing Newton-Rhapson alg.'
            opt_ll, opt_i = max(zero_intervals)
            opt_delta = 0.5 * (deltas[opt_i - 1] + deltas[opt_i])
            # Newton-Raphson
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    if method == 'REML':
                        new_opt_delta = optimize.newton(self._redll_, opt_delta, args=(eig_vals, sq_etas), tol=esp,
                                                        maxiter=100)
                    elif method == 'ML':
                        new_opt_delta = optimize.newton(self._dll_, opt_delta, args=(eig_vals, eig_vals_L, sq_etas),
                                                        tol=esp, maxiter=100)
            except Exception, err_str:
                if verbose:
                    print 'Problems with Newton-Raphson method:', err_str
                    print "Using the maximum grid value instead."
                    print 'opt_i:', opt_i
                    print 'Grid optimal delta:', opt_delta
                new_opt_delta = opt_delta
            # Validating the delta
            if opt_i > 1 and deltas[opt_i - 1] - esp < new_opt_delta < deltas[opt_i] + esp:
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            # Cheking lower boundary
            elif opt_i == 1 and 0.0 < new_opt_delta < deltas[opt_i] + esp:
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            # Cheking upper boundary
            elif opt_i == len(deltas) - 1 and new_opt_delta > deltas[opt_i - 1] - esp \
                        and not sp.isinf(new_opt_delta):
                opt_delta = new_opt_delta
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            else:
                if verbose:
                    print 'Local maximum outside of suggested possible areas?'
                    print 'opt_i:', opt_i
                    print 'Grid optimal delta:', opt_delta
                    print "Newton's optimal delta:", new_opt_delta
                    print 'Using the maximum grid value instead.'

            if verbose: print 'Done with Newton-Rahpson'
            if method == 'REML':
                opt_ll = self._rell_(opt_delta, eig_vals, sq_etas)
            elif method == 'ML':
                opt_ll = self._ll_(opt_delta, eig_vals, eig_vals_L, sq_etas)

            if opt_ll < max_ll:
                opt_delta = deltas[max_ll_i]
        else:
            if verbose: print 'No zero-interval was found, taking the maximum grid value.'
            opt_delta = deltas[max_ll_i]
            opt_ll = max_ll

        if verbose: print 'Finishing up.. calculating H_sqrt_inv.'
        l = sq_etas / (eig_vals + opt_delta)
        opt_vg = sp.sum(l) / p  # vg   
        opt_ve = opt_vg * opt_delta  # ve

        H_sqrt_inv = sp.mat(sp.diag(1.0 / sp.sqrt(eig_L['values'] + opt_delta)), dtype=dtype) * eig_L['vectors']
        # V = opt_vg * K + opt_ve * sp.eye(len(K))
        # H_sqrt = cholesky(V).T
        # H_sqrt_inv = H_sqrt.I
        X_t = H_sqrt_inv * X
        Y_t = H_sqrt_inv * self.Y
        (beta_est, mahalanobis_rss, rank, sigma) = linalg.lstsq(X_t, Y_t)
        x_beta = X * beta_est
        residuals = self.Y - x_beta
        rss = residuals.T * residuals
        # x_beta_var = sp.var(x_beta, ddof=1)
        # var_perc = x_beta_var / self.y_var
        res_dict = {'max_ll':opt_ll, 'delta':opt_delta, 'beta':beta_est, 've':opt_ve, 'vg':opt_vg,
            'rss':rss, 'mahalanobis_rss':mahalanobis_rss, 'H_sqrt_inv':H_sqrt_inv,
            'pseudo_heritability':1.0 / (1 + opt_delta)}

        if xs is not None and return_f_stat:
#            rss_ratio = h0_rss / rss_list
#            var_perc = 1 - 1 / rss_ratio
#            f_stats = (rss_ratio - 1) * n_p / float(q)

            h0_X = H_sqrt_inv * self.X
            (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y_t)
            f_stat = (h0_rss / mahalanobis_rss - 1) * p / xs.shape[1]
            res_dict['var_perc'] = 1.0 - mahalanobis_rss / h0_rss
            res_dict['f_stat'] = float(f_stat)
        if return_pvalue:
            p_val = stats.f.sf(f_stat, (xs.shape[1]), p)
            res_dict['p_val'] = float(p_val)
        return res_dict  # , lls, dlls, sp.log(deltas)



    def expedited_REML_t_test(self, snps, ngrids=50, llim= -4, ulim=10, esp=1e-6, verbose=True, eig_L=None):
        """
        Single SNP analysis, i.e. EMMA, but faster than R EMMA.
        """
        assert len(self.random_effects) == 2, "Expedited REMLE only works when we have exactly two random effects."
        K = self.random_effects[1][1]
        if eig_L is None:
            eig_L = self._get_eigen_L_(K)
        num_snps = len(snps)
        f_stats = sp.empty(num_snps)
        vgs = sp.empty(num_snps)
        ves = sp.empty(num_snps)
        max_lls = sp.empty(num_snps)
        var_perc = sp.empty(num_snps)
        rss_list = sp.empty(num_snps)
        betas = []
        p_vals = sp.empty(num_snps)

        # Run null model....
        log.info('Starting EMMA scan for %d SNPs' %len(snps))
        for i, snp in enumerate(snps):
            res = self.get_estimates(eig_L=eig_L, xs=sp.matrix(snp).T, ngrids=ngrids, llim=llim, ulim=ulim,
                               esp=esp, return_pvalue=True, return_f_stat=True)
            f_stats[i] = res['f_stat']
            vgs[i] = res['vg']
            ves[i] = res['ve']
            max_lls[i] = res['max_ll']
            var_perc[i] = res['var_perc']
            betas.append(map(float, list(res['beta'])))
            p_vals[i] = res['p_val']
            rss_list[i] = res['rss']
            if (i+1) % (num_snps / 10) == 0:
                perc = round(100.0 * i /num_snps)
                log.info('Performing EMMA (SNPs completed: %d %%)' % perc)

        return {'ps':p_vals, 'f_stats':f_stats, 'vgs':vgs, 'ves':ves, 'var_perc':var_perc,
            'max_lls':max_lls, 'betas':betas, 'rss':rss_list}



    def emmax_f_test_w_interactions(self, snps, int_af_threshold=15):
        """
        EMMAX implementation (in python)
        Single SNPs
        
        With interactions between SNP and possible cofactors.
        """
        assert len(self.random_effects) == 2, "Expedited REMLE only works when we have exactly two random effects."
        p_0 = len(self.X.T)
        n = self.n


        K = self.random_effects[1][1]
        eig_L = self._get_eigen_L_(K)
        res = self.get_expedited_REMLE(eig_L=eig_L)  # Get the variance estimates..
        delta = res['delta']
        print 'pseudo_heritability:', 1.0 / (1 + delta)
        H_sqr = res['H_sqrt']
        H_sqrt_inv = H_sqr.I
        Y = H_sqrt_inv * self.Y  # The transformed outputs.
        h0_X = H_sqrt_inv * self.X
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        h0_betas = map(float, list(h0_betas))
        f_stats = []
        rss_list = []
        betas_list = []
        p_vals = []
        var_perc = []
        cofactors = sp.array(self.cofactors)
        num_cof = len(cofactors)
        low_int_af_count = 0
        for snp in snps:
            snp_a = sp.array(snp)
            if sp.sum(snp_a) > int_af_threshold:
                interactions = cofactors * snp_a
                interactions = interactions[sp.sum(interactions, axis=1) > int_af_threshold]
                low_int_af_count += num_cof - len(interactions)
                X = sp.hstack([h0_X, H_sqrt_inv * sp.matrix(snp_a).T, H_sqrt_inv * sp.matrix(interactions).T])
            else:
                low_int_af_count += num_cof
                X = sp.hstack([h0_X, H_sqrt_inv * sp.matrix(snp_a).T])
            (betas, rss, p, sigma) = linalg.lstsq(X, Y)
            q = p - p_0
            n_p = n - p
            if not rss:
                if q == 0:
                    print 'No predictability in the marker, moving on...'
                    p_vals.append(1)
                    f_stats.append(0)
                    rss_list.append(h0_rss)
                    betas_list.append(h0_betas)
                    var_perc.append(0)
                    continue
                else:
                    res = (Y - X * betas)
                    rss = res.T * res
            f_stat = ((h0_rss - rss) / q) / (rss / n_p)
            p_val = stats.f.sf(f_stat, q, n_p)
            p_vals.append(p_val[0])
            f_stats.append(f_stat[0])
            rss_list.append(rss[0])
            betas_list.append(map(float, list(betas)))
            var_perc.append(float(1 - rss / h0_rss))

        print 'low_int_af_count:', low_int_af_count

        return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
            'delta':delta, 'pseudo_heritability': 1.0 / (1 + delta), 'var_perc':var_perc}



    def emmax_anova_f_test(self, snps):
        """
        EMMAX implementation (in python)
        Single SNPs
        
        With interactions between SNP and possible cofactors.
        """
        K = self.random_effects[1][1]
        eig_L = self._get_eigen_L_(K)
        res = self.get_expedited_REMLE(eig_L=eig_L)  # Get the variance estimates..
        print 'pseudo_heritability:', res['pseudo_heritability']

        r = self._emmax_anova_f_test_(snps, res['H_sqrt'])
        r['pseudo_heritability'] = res['pseudo_heritability']
        r['max_ll'] = res['max_ll']
        return r



    def _emmax_anova_f_test_(self, snps, H_sqrt, verbose=True):
        """
        EMMAX implementation (in python)
        Single SNPs
        
        With interactions between SNP and possible cofactors.
        """
        n = self.n
        p_0 = len(self.X.T)

        H_sqrt_inv = H_sqrt.I
        Y = H_sqrt_inv * self.Y  # The transformed outputs.
        h0_X = H_sqrt_inv * self.X
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        h0_betas = map(float, list(h0_betas))
        num_snps = len(snps)
        rss_list = sp.repeat(h0_rss, num_snps)
        betas_list = [h0_betas] * num_snps
        var_perc = sp.zeros(num_snps)
        f_stats = sp.zeros(num_snps)
        dfs = sp.zeros(num_snps)
        p_vals = sp.ones(num_snps)

        for i, snp in enumerate(snps):
            groups = sp.unique(snp)
            q = len(groups) - 1  # Null model has 1 df.
            p = p_0 + q
            n_p = n - p
            l = []
            for g in groups:
                l.append(sp.int8(snp == g))
            X = sp.mat(l) * (H_sqrt_inv.T)
            (betas, rss, p, sigma) = linalg.lstsq(X.T, Y, overwrite_a=True)
            if not rss:
                log.debug('No predictability in the marker, moving on...')
                continue
            rss_list[i] = rss[0]
            betas_list[i] = map(float, list(betas))
            rss_ratio = h0_rss / rss
            var_perc[i] = 1 - 1 / rss_ratio
            f_stat = (rss_ratio - 1) * n_p / float(q)
            p_vals[i] = stats.f.sf([f_stat], q, n_p)[0]
            f_stats[i] = f_stat
            dfs[i] = n_p
            if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()

        if num_snps >= 10:
            sys.stdout.write('\n')
        # rss_ratio = h0_rss / rss_list
        # var_perc = 1 - 1 / rss_ratio
        # f_stats = (rss_ratio - 1) * n_p / float(q)
        # p_vals = stats.f.sf(f_stats, q, n_p)


        return {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'betas':betas_list,
            'var_perc':var_perc, 'dfs':dfs}



    def emmax_permutations(self, snps, num_perm, method='REML'):
        """
        EMMAX permutation test
        Single SNPs
        
        Returns the list of max_pvals and max_fstats 
        """
        K = self.random_effects[1][1]
        eig_L = self._get_eigen_L_(K)
        # s = time.time()
        res = self.get_estimates(eig_L=eig_L, method=method)  # Get the variance estimates..
        # print 'Took % .6f secs.' % (time.time() - s)
        # print 'pseudo_heritability:', res['pseudo_heritability']
        q = 1  # Single SNP is being tested
        p = len(self.X.T) + q
        n = self.n
        n_p = n - p
        H_sqrt_inv = res['H_sqrt_inv']
        Y = H_sqrt_inv * self.Y  # The transformed outputs.
        h0_X = H_sqrt_inv * self.X
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        Y = Y - h0_X * h0_betas
        num_snps = len(snps)
        max_fstat_list = []
        min_pval_list = []
        chunk_size = len(Y)
        Ys = sp.mat(sp.zeros((chunk_size, num_perm)))
        for perm_i in range(num_perm):
            # print 'Permutation nr. % d' % perm_i
            sp.random.shuffle(Y)
            Ys[:, perm_i] = Y

        min_rss_list = sp.repeat(h0_rss, num_perm)
        for i in range(0, num_snps, chunk_size):  # Do the dot-product in chuncks!
            snps_chunk = sp.matrix(snps[i:i + chunk_size])
            Xs = snps_chunk * (H_sqrt_inv.T)
            Xs = Xs - sp.mat(sp.mean(Xs, axis=1))
            for j in range(len(Xs)):
                (betas, rss_list, p, sigma) = linalg.lstsq(Xs[j].T, Ys, overwrite_a=True)
                for k, rss in enumerate(rss_list):
                    if not rss:
                        log.debug('No predictability in the marker, moving on...')
                        continue
                    if min_rss_list[k] > rss:
                        min_rss_list[k] = rss
                if num_snps >= 10 and (i + j + 1) % (num_snps / num_perm) == 0:  # Print dots
                    sys.stdout.write('.')
                    sys.stdout.flush()

        if num_snps >= 10:
            sys.stdout.write('\n')
        min_rss = min(rss_list)
        max_f_stats = ((h0_rss / min_rss_list) - 1.0) * n_p / float(q)
        min_pvals = (stats.f.sf(max_f_stats, q, n_p))


        res_d = {'min_ps':min_pvals, 'max_f_stats':max_f_stats}
        return res_d



    def emmax_f_test(self, snps,num_snps, snp_priors=None, Z=None, with_betas=False, method='REML',
            eig_L=None, eig_R=None, emma_num=100, progress_file_writer=None):
        """
        EMMAX implementation (in python)
        Single SNPs
        """
        if not eig_L:
            log.info('Calculating the eigenvalues of K')
            s0 = time.time()
            eig_L = self._get_eigen_L_()
            log.info('Done.')
            log.info('Took %0.2f seconds' % (time.time() - s0))
        if not eig_R:
            log.info("Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'")
            s0 = time.time()
            eig_R = self._get_eigen_R_(X=self.X)
            log.info('Done')
            log.info('Took %0.2f seconds' % (time.time() - s0))

        log.info('Getting variance estimates')
        s0 = time.time()
        res = self.get_estimates(eig_L, method=method, eig_R=eig_R)  # Get the variance estimates..
        log.info('Done.')
        log.info('Took %0.2f seconds' % (time.time() - s0))
        log.info('pseudo_heritability:', res['pseudo_heritability'])

        s0 = time.time()
        r = self._emmax_f_test_(snps, num_snps, res['H_sqrt_inv'], snp_priors=snp_priors, Z=Z, with_betas=with_betas,
                progress_file_writer=progress_file_writer, emma_num=emma_num, eig_L=eig_L)
        log.info('Took %0.2f seconds' % (time.time() - s0))
        r['pseudo_heritability'] = res['pseudo_heritability']
        r['ve'] = res['ve']
        r['vg'] = res['vg']
        r['max_ll'] = res['max_ll']

        return r




    def _emmax_f_test_(self, genotype, num_snps,H_sqrt_inv, snp_priors=None, verbose=True, return_transformed_snps=False,
            Z=None, with_betas=False, progress_file_writer=None, emma_num=100, eig_L=None, **kwargs):
        """
        EMMAX implementation (in python)
        Single SNPs
        
        Methods:
            normal - Normal regression
            qr - Uses QR decomposition to speed up regression with many co-factors.
            
        """
        snps = genotype.get_snps_iterator(is_chunked=True)
        dtype = 'single'
        q = 1  # Single SNP is being tested
        p = len(self.X.T) + q
        n = self.n
        n_p = n - p

        h0_X = sp.mat(H_sqrt_inv * self.X, dtype=dtype)
        Y = H_sqrt_inv * self.Y  # The transformed outputs.
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        assert len(h0_rss) > 0, 'WTF?'
        Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
        h0_betas = map(float, list(h0_betas))

        if Z is not None:
            H_sqrt_inv = H_sqrt_inv * Z

        if not with_betas:
            (Q, R) = qr_decomp(h0_X)  # Do the QR-decomposition for the Gram-Schmidt process.
            Q = sp.mat(Q)
            Q2 = Q * Q.T
            M = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)
        else:
            betas_list = [h0_betas] * num_snps
            M = H_sqrt_inv.T

        rss_list = sp.repeat(h0_rss, num_snps)
        if return_transformed_snps:
            t_snps = []
        if snp_priors is not None:
            snp_priors = sp.array(snp_priors)
            log_h0_rss = sp.log(h0_rss)
            log_bfs = sp.zeros(num_snps)  # Bayes factors
        if not progress_file_writer == None:
            progress_file_writer.update_progress_bar(progress=0.25, task_status='Starting AMM scan')
        log.info('Starting AMM scan',extra={'progress':25})
        i = 0

        for snps_chunk in snps:
            if len(snps_chunk) == 0:
                continue
            Xs = sp.matrix(snps_chunk,dtype=dtype) * M
            for X_j in Xs:
                if with_betas:
                    (betas, rss, p, sigma) = linalg.lstsq(sp.hstack([h0_X, X_j.T]), Y, overwrite_a=True)
                    if rss:
                        betas_list[i] = map(float, list(betas))
                else:
                    (betas, rss, p, sigma) = linalg.lstsq(X_j.T, Y, overwrite_a=True)
                if rss:
                    rss_list[i] = rss[0]
    
                if snp_priors is not None:
                    log_bfs[i] = log_h0_rss - sp.log(rss)  # -(1/2)*log(n)
                if (i+1) % (num_snps / 10) == 0:
                    perc = round(100.0 * i /num_snps)
                    log.info('Performing AMM (SNPs completed: %d %%)' % perc,extra={'progress':25 + 55*perc/100}) 
                i += 1
        rss_ratio = h0_rss / rss_list
        var_perc = 1 - 1 / rss_ratio
        # assert sp.all(var_perc < 1.01), '%f\n%s\n%s' % (h0_rss, str(var_perc[var_perc < 1.01]), str(rss_list[var_perc < 1.01]))
        f_stats = (rss_ratio - 1) * n_p / float(q)
        p_vals = stats.f.sf(f_stats, q, n_p)

        res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
            'h0_rss':h0_rss, 'h0_betas':h0_betas}
        if with_betas:
            res_d['betas'] = betas_list
        if return_transformed_snps:
            res_d['t_snps'] = t_snps
        if snp_priors is not None:
            bfs = sp.exp((log_bfs * n - sp.log(n)) * 1 / 2)  # Bayes factors
            res_d['bfs'] = bfs
            pos = bfs * snp_priors / (1 - snp_priors)  # Posterior odds
            res_d['pos'] = pos
            ppas = pos / (1 + pos)  # Posterior probablities of association
            res_d['ppas'] = ppas

        if emma_num > 0:
            pval_indices = sorted(zip(res_d['ps'], range(num_snps)))[:emma_num]

            if not progress_file_writer == None:
                progress_file_writer.update_progress_bar(task_status='Replacing smallest p-values with exact mixed model p-values')
            log.info('Updating p-values using EMMA for the smallest %d p-values.' % len(pval_indices),extra={'progress':81})
            l = map(list, zip(*pval_indices))
            chr_pos = [genotype.get_chr_pos_from_index(ix) for ix in l[1] ]
            top_snps = genotype.get_snps_from_pos(chr_pos)[1]
            top_emma_res = self.expedited_REML_t_test(top_snps, eig_L=eig_L)
            for pi, p, f, r, v in it.izip(l[1], top_emma_res['ps'], top_emma_res['f_stats'], \
                        top_emma_res['rss'], top_emma_res['var_perc']):
                res_d['ps'][pi] = p
                res_d['f_stats'][pi] = f
                res_d['rss'][pi] = r
                res_d['var_perc'][pi] = v
        return res_d



    def _emmax_GxT_f_test_(self, snps, H_sqrt_inv, T, Z, snp_priors=None, verbose=True, **kwargs):
        """
        EMMAX implementation (in python)
        Single SNPs
        
        Methods:
            normal - Normal regression
            qr - Uses QR decomposition to speed up regression with many co-factors.
            
        """
        dtype = 'single'
        n = self.n
        num_snps = len(snps)

        h0_X = sp.mat(H_sqrt_inv * self.X, dtype=dtype)
        Y = H_sqrt_inv * self.Y  # The transformed outputs.
        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
        Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
        h0_betas = map(float, list(h0_betas))
        T_flat = sp.array(T).flatten()

        betas_g_list = [h0_betas] * num_snps
        betas_gt_list = [h0_betas] * num_snps
        M = H_sqrt_inv.T
        rss_g_list = sp.repeat(h0_rss, num_snps)
        rss_gt_list = sp.repeat(h0_rss, num_snps)
        chunk_size = len(Y)
        for i in range(0, num_snps, chunk_size):  # Do the dot-product in chuncks!
            snps_chunk = sp.matrix(snps[i:i + chunk_size], dtype=dtype) * Z.T
            GT = sp.matrix(sp.array(snps_chunk) * T_flat) * M
            Xs = snps_chunk * M

            for j, X_j, in enumerate(Xs):
                GT_j = GT[j]
                (betas_g, rss_g, p, sigma) = linalg.lstsq(sp.hstack([h0_X, X_j.T]), Y, overwrite_a=True)
                if rss_g:
                    betas_g_list[i + j] = map(float, list(betas_g))
                    rss_g_list[i + j] = rss_g[0]

                    (betas_gt, rss_gt, p, sigma) = linalg.lstsq(sp.hstack([h0_X, X_j.T, GT_j.T]), Y,
                                        overwrite_a=True)
                    if rss_gt:
                        betas_gt_list[i + j] = map(float, list(betas_gt))
                        rss_gt_list[i + j] = rss_gt[0]

                if verbose and num_snps >= 10 and (i + j + 1) % (num_snps / 10) == 0:  # Print dots
                    sys.stdout.write('.')
                    sys.stdout.flush()

        if verbose and num_snps >= 10:
            sys.stdout.write('\n')

        # FIRST THE G
        q = 1  # Single SNP is being tested
        p = len(self.X.T) + q
        n_p = n - p

        rss_g_ratio = h0_rss / rss_g_list
        var_perc_g = 1 - 1 / rss_g_ratio
        f_stats_g = (rss_g_ratio - 1) * n_p / float(q)
        p_vals_g = stats.f.sf(f_stats_g, q, n_p)

        g_res_d = {'ps':p_vals_g, 'f_stats':f_stats_g, 'rss':rss_g_list, 'var_perc':var_perc_g,
            'h0_rss':h0_rss, 'h0_betas':h0_betas, 'betas': betas_g_list}


        # G+GxT
        q = 2
        p = len(self.X.T) + q
        n_p = n - p

        rss_gt_ratio = h0_rss / rss_gt_list
        var_perc_gt = 1 - 1 / rss_gt_ratio
        f_stats_gt = (rss_gt_ratio - 1) * n_p / float(q)
        p_vals_gt = stats.f.sf(f_stats_gt, q, n_p)

        gt_res_d = {'ps':p_vals_gt, 'f_stats':f_stats_gt, 'rss':rss_gt_list, 'var_perc':var_perc_gt,
                'betas': betas_gt_list}


        # G+GXT vs G
        q = 1  # Single interaction is being tested
        p = len(self.X.T) + q
        n_p = n - p

        rss_gt_g_ratio = rss_g_list / rss_gt_list
        var_perc_gt_g = 1 - 1 / rss_gt_g_ratio
        f_stats_gt_g = (rss_gt_g_ratio - 1) * n_p / float(q)
        p_vals_gt_g = stats.f.sf(f_stats_gt_g, q, n_p)

        gt_g_res_d = {'ps':p_vals_gt_g, 'f_stats':f_stats_gt_g, 'var_perc':var_perc_gt_g}
        return {'g_res':g_res_d, 'gt_res':gt_res_d, 'gt_g_res':gt_g_res_d}




    def emmax_two_snps(self, snps, verbose=True):
        """
        Every pair of SNPs, Vincent's results.
        """
        import util
        K = self.random_effects[1][1]
        eig_L = self._get_eigen_L_(K)
        num_snps = len(snps)

        res = self.get_expedited_REMLE(eig_L=eig_L)  # Get the variance estimates..
        delta = res['delta']
        pseudo_her = 1.0 / (1 + delta)
        H_sqrt = res['H_sqrt']
        H_sqrt_inv = H_sqrt.I
        emmax_res = self._emmax_f_test_(snps, H_sqrt, verbose=False, return_transformed_snps=True)
        t_snps = emmax_res['t_snps']
        full_rss = emmax_res['h0_rss']
        h0_X = H_sqrt_inv * self.X
        Y = H_sqrt_inv * self.Y  # The transformed outputs.

        # Contructing result containers        
        p3_f_stat_array = sp.zeros((num_snps, num_snps))
        p3_p_val_array = sp.ones((num_snps, num_snps))
        p4_f_stat_array = sp.zeros((num_snps, num_snps))
        p4_p_val_array = sp.ones((num_snps, num_snps))
        f_stat_array = sp.zeros((num_snps, num_snps))
        p_val_array = sp.ones((num_snps, num_snps))
        rss_array = sp.repeat(full_rss, num_snps * num_snps)
        rss_array.shape = (num_snps, num_snps)
        var_perc_array = sp.zeros((num_snps, num_snps))
        haplotype_counts = [[{} for j in range(i + 1)] for i in range(num_snps)]

        # Fill the diagonals with the single SNP emmax
        for i, snp in enumerate(snps):
            hap_set, hap_counts, haplotypes = genotype.get_haplotypes([snp], self.n,
                                        count_haplotypes=True)
            d = {'num_haps':hap_counts}
            for hap, hap_c in zip(hap_set, hap_counts):
                d[hap] = hap_c
            haplotype_counts[i][i] = d
            p_val_array[i, i] = emmax_res['ps'][i]
            p3_p_val_array[i, i] = p_val_array[i, i]
            p4_p_val_array[i, i] = p_val_array[i, i]
            f_stat_array[i, i] = emmax_res['f_stats'][i]
            p3_f_stat_array[i, i] = f_stat_array[i, i]
            p4_f_stat_array[i, i] = f_stat_array[i, i]
            rss_array[i, i] = emmax_res['rss'][i]
            var_perc_array[i, i] = emmax_res['var_perc'][i]


        identical_snp_count = 0
        no_interaction_count = 0
#        p = len(self.X.T) + q + 1 #Adding one SNP as a cofactor.
#        n_p = self.n - p
        for i, snp1 in enumerate(snps):
            snp1 = snps[i]
            for j in range(i):
                snp2 = snps[j]
                if i == j: continue  # Skip diagonals..

                # Haplotype counts 
                hap_set, hap_counts, haplotypes = genotype.get_haplotypes([snp1, snp2], self.n,
                                            count_haplotypes=True)
                groups = set(haplotypes)
                d = {'num_haps':len(hap_counts)}
                for hap, hap_c in zip(hap_set, hap_counts):
                    d[hap] = hap_c
                haplotype_counts[i][j] = d

                # Fill the upper right part with more powerful of two tests.

                if emmax_res['ps'][i] < emmax_res['ps'][j]:
                    rss_array[j, i] = emmax_res['rss'][i]
                    max_i = i
                else:
                    rss_array[j, i] = emmax_res['rss'][j]
                    max_i = j

                if d['num_haps'] == 2:
                    identical_snp_count += 1
                    continue
                elif d['num_haps'] == 3:
                    n_p = self.n - 3
                    no_interaction_count += 1
                    # Do ANOVA
                    l = []
                    for g in groups:
                        l.append(sp.int8(haplotypes == g))
                    X = sp.mat(l) * (H_sqrt_inv.T)
                    (betas, rss, p, sigma) = linalg.lstsq(X.T, Y)
                    rss_array[i, j] = rss[0]
                    var_perc_array[i, j] = 1 - rss / full_rss
                    f_stat = (rss_array[j, i] / rss - 1) * n_p  # FINISH
                    p_val = stats.f.sf([f_stat], 1, n_p)[0]
                    p3_f_stat_array[j, i] = f_stat
                    p3_f_stat_array[i, j] = f_stat
                    p3_p_val_array[j, i] = p_val
                    p3_p_val_array[i, j] = p_val

                    f_stat = ((full_rss - rss) / 2) / (rss / n_p)
                    f_stat_array[j, i] = f_stat
                    p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]


                elif d['num_haps'] == 4:  # With interaction
                    n_p = self.n - 3
                    # Do additive model
                    # snp_mat = H_sqrt_inv * sp.mat([snp1, snp2]).T #Transformed inputs
                    snp_mat = sp.mat([t_snps[i], t_snps[j]]).T  # Transformed inputs
                    X = sp.hstack([h0_X, snp_mat])
                    (betas, rss, rank, s) = linalg.lstsq(X, Y)
                    f_stat = (rss_array[j, i] / rss - 1) * n_p  # Compared to only one SNP
                    p_val = stats.f.sf([f_stat], 1, n_p)[0]
                    rss_array[i, j] = rss
                    p3_f_stat_array[j, i] = f_stat
                    p3_p_val_array[j, i] = p_val

#                    v_f_stat_array[j, i] = f_stat
#                    v_p_val_array[j, i] = stats.f.sf([f_stat], 1, n_p)[0]

                    f_stat = ((full_rss - rss) / 2) / (rss / n_p)  # Compared to only the intercept
                    f_stat_array[j, i] = f_stat
                    p_val_array[j, i] = stats.f.sf([f_stat], 2, n_p)[0]

                    # Generate the interaction, and test it.
                    i_snp = snp1 & snp2
                    snp_mat = H_sqrt_inv * sp.mat([i_snp]).T
                    X = sp.hstack([h0_X, sp.mat([t_snps[max_i]]).T, snp_mat])
                    (betas, rss, rank, s) = linalg.lstsq(X, Y)
                    f_stat = (rss_array[j, i] / rss - 1) * n_p  # Compared to only one SNP
                    p_val = stats.f.sf([f_stat], 1, n_p)[0]
                    p3_f_stat_array[i, j] = f_stat
                    p3_p_val_array[i, j] = p_val


                    # full model p-value
                    n_p = self.n - 4
                    l = []
                    for g in groups:
                        l.append(sp.int8(haplotypes == g))
                    X = sp.mat(l) * (H_sqrt_inv.T)
                    (betas, rss, p, sigma) = linalg.lstsq(X.T, Y)

                    f_stat = ((rss_array[j, i] - rss) / 2) / (rss / n_p)  # Compared to only one SNP
                    p_val = stats.f.sf([f_stat], 2, n_p)[0]
                    p4_f_stat_array[j, i] = f_stat
                    p4_p_val_array[j, i] = p_val

                    f_stat = (rss_array[i, j] / rss - 1) * n_p  # Compared to two SNPs
                    p4_f_stat_array[i, j] = f_stat
                    p4_p_val_array[i, j] = stats.f.sf([f_stat], 1, n_p)[0]

                    f_stat = ((full_rss - rss) / 3) / (rss / n_p)  # Compared to only intercept
                    f_stat_array[i, j] = f_stat
                    p_val_array[i, j] = stats.f.sf([f_stat], 3, n_p)[0]
                    rss_array[j, i] = rss

            if num_snps >= 10 and (i + 1) % (num_snps / 10) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()

        print no_interaction_count, identical_snp_count

        #FINISH res dict!!!
        res_dict = {'p3_ps':p3_p_val_array, 'p3_f_stats':p3_f_stat_array, 'p4_ps':p4_p_val_array,
            'p4_f_stats':p4_f_stat_array, 'rss':rss_array, 'var_perc':var_perc_array,
            'pseudo_heritability': pseudo_her, 'haplotype_counts':haplotype_counts,
            'f_stats':f_stat_array, 'ps':p_val_array}
        return res_dict


    def calc_statistics(self):
        """
        Returns all sorts of statistics used in stepwise regression
        
        Log Likelihood, BIC, modified BIC, extended BIC, RSS, mahalnobis RSS 
        """
        pass





