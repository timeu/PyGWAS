"""
Contains functions to perform Genome Wide Association Mapping (GWAS)

Author: Bjarni J. Vilhjalmsson
Email: bjarni.vilhjalmsson@gmail.com

"""
import logging
from linear_models import LinearModel,LinearMixedModel
import numpy as np
import scipy as sp
from timer import Timer

log = logging.getLogger(__name__)


def _calc_bic_(ll, num_snps, num_par, n):
    bic = -2 * (ll) + num_par * sp.log(n)
    extended_bic = bic + \
        2 * _log_choose_(num_snps, num_par - 2)
    modified_bic = bic + \
        2 * (num_par) * sp.log(num_snps / 2.2 - 1)
    return (bic, extended_bic, modified_bic)

def _log_choose_(n, k):
     # log_fact = lambda n: n * sp.log(n) - n + (sp.log(n * (1 + 4 * n * (1 + 2 * n))) / 6) + sp.log(sp.pi) / 2
     # Srinivasa Ramanujan approximation of log(n!)
    if k == 0 or n == k:
        return 0
    if n < k:
        raise Exception('Out of range.')
    return sp.sum(map(sp.log, range(n, n - k, -1))) - sp.sum(map(sp.log, range(k, 0, -1)))


def mm_lrt_test(y, K):
    """
    Likelihood ratio test for whether the data (y) fits a mixed model with 
    two random terms fits significantly better.
    """
    lm = LinearModel(y)
    lmm = LinearMixedModel(y)
    lmm.add_random_effect(K)
    lmm_res = lmm.get_ML()
    ll0 = lm.get_ll()
    ll1 = lmm_res['max_ll']
    D = 2 * (ll1 - ll0)
    pval = sp.chi2.sf(D, 1)
    return {'pval':pval, 'lrt_stat':D}


def emma(snps, phenotypes, K, cofactors=None):
    """
    Run EMMAX
    """
    lmm = LinearMixedModel(phenotypes)
    lmm.add_random_effect(K)
    if cofactors:
        for cofactor in cofactors:
            lmm.add_factor(cofactor)

    log.info("Running EMMA (python)")
    t = Timer()
    res = lmm.expedited_REML_t_test(snps)
    log.info('Took: %s' % t.stop(True))
    return res



def emmax(snps, phenotypes, K, cofactors=None, Z=None, with_betas=False, progress_file_writer=None, emma_num=0):
    """
    Run EMMAX
    """
    lmm = LinearMixedModel(phenotypes)
    if Z != None:
        lmm.add_random_effect(Z * K * Z.T)
        if cofactors:
            for cofactor in cofactors:
                lmm.add_factor(Z * cofactor)
    else:
        lmm.add_random_effect(K)
        if cofactors:
            for cofactor in cofactors:
                lmm.add_factor(cofactor)

    if not progress_file_writer == None:
        progress_file_writer.update_progress_bar(0.4, 'Performing EMMAX')
    log.info("Running EMMAX")
    t = Timer()
    res = lmm.emmax_f_test(snps, Z=Z, with_betas=with_betas, progress_file_writer=progress_file_writer,
            emma_num=emma_num)
    log.info('Took:%s' % t.stop(True))
    return res



def emmax_perm_test(snps, phenotypes, K, num_perm=100):
    """
    Run EMMAX
    """
    lmm = LinearMixedModel(phenotypes)
    lmm.add_random_effect(K)
    log.info("Running %d EMMAX-permutation (writes %d dots)" % (num_perm, num_perm))
    t = Timer()
    res = lmm.emmax_permutations(snps, num_perm)
    p_f_list = zip(res['min_ps'], res['max_f_stats'])
    p_f_list.sort()
    threshold = p_f_list[len(p_f_list) / 20]
    res['threshold_05'] = threshold
    log.info('Tresholds should be:%s' % threshold)
    log.info('Took: %s' % t.stop(True))
    return res


def anova(snps, phenotypes):
    """
    Run EMMAX
    """
    lmm = LinearModel(phenotypes)

    log.info("Running ANOVA")
    t = Timer()
    res = lmm.anova_f_test(snps)
    log.info('Took: %s' % t.stop(True))
    return res



def emmax_anova(snps, phenotypes, K):
    """
    Run EMMAX
    """
    lmm = LinearMixedModel(phenotypes)
    lmm.add_random_effect(K)

    log.info("Running EMMAX-ANOVA")
    t = Timer()
    res = lmm.emmax_anova_f_test(snps)
    log.info('Took:%s' % t.stop(True))
    return res


def lin_reg_step(phen_vals, sd, cof_chr_pos_list, progress_file_writer=None, plot_prefix=None):
    """
    Standard linear regression single SNPs..
    
    Returns various stats useful for stepwise regression.
    """
    import bisect
    t = Timer()
    lm = LinearModel(phen_vals)

    h0_rss = lm.get_rss()
    step_dict = {}
    step_dict['h0_rss'] = h0_rss

    
    log.info('Looking up cofactors')
    cof_indices,cof_snps  = sd.get_snps_from_pos(cof_chr_pos_list)
    lm.set_factors(cof_snps)

    if not progress_file_writer == None:
        progress_file_writer.update_progress_bar(progress=0.20, task_status='Performing linear regression')
        progress_file_writer.set_step(0.05)


    log.info('Performing linear regression.')
    r = lm.fast_f_test(sd, progress_file_writer=progress_file_writer)
    min_pval_i = sp.argmin(r['ps'])
    step_dict['min_pval_i'] = min_pval_i
    step_dict['min_pval'] = r['ps'][min_pval_i]
    step_dict['mahalanobis_rss'] = r['rss'][min_pval_i]
    step_dict['min_pval_chr_pos'] = sd.get_chr_pos_from_index(min_pval_i)

    num_snps = sd.num_snps
    num_par = lm.X.shape[1] + 1
    step_dict['num_snps'] = num_snps
    step_dict['num_par'] = num_par
    rss = lm.get_rss()
    ll = lm.get_ll(rss)
    (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n)  # Calculate the BICs
    step_dict['ebic'] = extended_bic
    step_dict['mbic'] = modified_bic
    step_dict['bic'] = bic

    step_dict['rss'] = rss
    perc_var_expl = 1.0 - (rss / h0_rss)
    step_dict['perc_var_expl'] = perc_var_expl

    # Calculate maximum cofactor p-value\
    log.info('Updating the cofactor p-values')
    cof_pvals = []
    cof_chrom_pos_pval_list = []
    for i, snp in enumerate(cof_snps):
        t_cofactors = cof_snps[:]
        del t_cofactors[i]
        lm.set_factors(t_cofactors)
        res = lm.fast_f_test([snp])
        cof_pval = res['ps'][0]
        cof_pvals.append(cof_pval)
        cof_chrom_pos_pval_list.append((cof_chr_pos_list[i][0], cof_chr_pos_list[i][1], -math.log10(cof_pval)))
    for i, pval in zip(cof_indices, cof_pvals):
        r['ps'][i] = pval
    if len(cof_pvals):
        step_dict['max_cof_pval'] = max(cof_pvals)
    else:
        step_dict['max_cof_pval'] = 0.0
    log.info('Took:%s' % t.stop(True))

    return {'stats':step_dict, 'res':r,'cof_chrom_pos_pval_list':cof_chrom_pos_pval_list}


def calculate_pseudo_heritability(phen_vals,K):
    """
    Retrieve pseudo heritability for phenotype

    """
    lmm = LinearMixedModel(phen_vals)
    lmm.add_random_effect(K)
    reml_dict = lmm.get_REML()
    return reml_dict['pseudo_heritability']


def emmax_step(phen_vals, genotype, K, cof_chr_pos_list, eig_L=None, eig_R=None, progress_file_writer=None, plot_prefix=None,
        emma_num=100):
    """
    EMMAX single SNPs
    
    Returns various stats useful for stepwise regression.
    """
    import bisect
    t = Timer()
    lmm = LinearMixedModel(phen_vals)
    lmm.add_random_effect(K)

    if not eig_L:
        log.debug('Calculating the eigenvalues of K')
        eig_L = lmm._get_eigen_L_()
        log.debug('Done.')
    reml_dict = lmm.get_REML(eig_L=eig_L)
    h0_rss = float(reml_dict['rss'])
    cof_indices = []
    """
    chrom_pos_list = sd.get_chr_pos_list()
    log.debug('Looking up cofactors')

    for chrom_pos in cof_chr_pos_list:
        i = bisect.bisect(chrom_pos_list, chrom_pos) - 1
        assert chrom_pos_list[i] == chrom_pos, 'Cofactor missing??'
        cof_indices.append(i)"""
    num_snps = genotype.num_snps
    cof_snps_ix,cof_snps = genotype.get_snps_from_pos(cof_chr_pos_list)
    lmm.set_factors(cof_snps)
    if not eig_R:
        log.debug("Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'")
        eig_R = lmm._get_eigen_R_()
        log.debug('Done')


    log.debug('Getting variance components estimates')
    reml_dict = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
    ml_dict = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
    log.debug('Done.')
    log.info('pseudo_heritability: %s' % reml_dict['pseudo_heritability'])
    H_sqrt_inv = reml_dict['H_sqrt_inv']

    if not progress_file_writer == None:
        progress_file_writer.update_progress_bar(progress=0.20, task_status='Performing AMM')
        progress_file_writer.set_step(0.05)


    r = lmm._emmax_f_test_(genotype, num_snps,H_sqrt_inv, progress_file_writer=progress_file_writer, emma_num=emma_num)
    min_pval_i = sp.argmin(r['ps'])
    step_dict = {}
    step_dict['min_pval_i'] = min_pval_i
    step_dict['min_pval'] = r['ps'][min_pval_i]
    step_dict['mahalanobis_rss'] = r['rss'][min_pval_i]
    step_dict['min_pval_chr_pos'] = genotype.get_chr_pos_from_index(min_pval_i)
    step_dict['h0_rss'] = h0_rss

    ll = ml_dict['max_ll']
    step_dict['max_ll'] = ll

    num_par = lmm.X.shape[1] + 1
    step_dict['num_snps'] = num_snps
    step_dict['num_par'] = num_par
    (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n)  # Calculate the BICs
    step_dict['ebic'] = extended_bic
    step_dict['mbic'] = modified_bic
    step_dict['bic'] = bic

    step_dict['ve'] = reml_dict['ve']
    step_dict['vg'] = reml_dict['vg']
    p_her = reml_dict['pseudo_heritability']
    step_dict['pseudo_heritability'] = p_her

    step_dict['rss'] = float(reml_dict['rss'])
    step_dict['reml_mahalanobis_rss'] = reml_dict['mahalanobis_rss']
    perc_var_expl = 1.0 - (step_dict['rss'] / h0_rss)
    step_dict['perc_var_expl'] = perc_var_expl
    step_dict['remain_perc_gen_var'] = (1 - perc_var_expl) * p_her
    step_dict['remain_perc_err_var'] = (1 - perc_var_expl) * (1 - p_her)

    # Calculate maximum cofactor p-value\
    log.debug('Updating the cofactor p-values')
    cof_pvals = []
    cof_chrom_pos_pval_list = []
    for i, snp in enumerate(cof_snps):
        t_cofactors = cof_snps[:]
        del t_cofactors[i]
        lmm.set_factors(t_cofactors)
        res = lmm._emmax_f_test_([snp], H_sqrt_inv, emma_num=0)
        cof_pval = res['ps'][0]
        cof_pvals.append(cof_pval)
        cof_chrom_pos_pval_list.append((cof_chr_pos_list[i][0], cof_chr_pos_list[i][1], -math.log10(cof_pval)))
    for i, pval in zip(cof_indices, cof_pvals):
        r['ps'][i] = pval
    if len(cof_pvals):
        step_dict['max_cof_pval'] = max(cof_pvals)
    else:
        step_dict['max_cof_pval'] = 0.0
    # step_dict['cofactor_snps'] = cof_snps
    log.info('Took %s' % t.stop(True))

    return {'stats':step_dict, 'res':r, 'upd_H_sqrt_inv':H_sqrt_inv,'cof_chrom_pos_pval_list':cof_chrom_pos_pval_list}




def emmax_step_wise(phenotypes, K, sd=None, num_steps=10, file_prefix=None, forward_backwards=True,
        local=False, cand_gene_list=None, plot_xaxis=True, with_qq_plots=True, sign_threshold=None,
        log_qq_max_val=5, highlight_loci=None, save_pvals=False, pval_file_prefix=None, snp_priors=None,
        K2=None, snp_choose_criteria='pval', emma_num=0, markersize=3, chrom_col_map=None, **kwargs):
    """
    Run step-wise EMMAX forward-backward.
    """
    import gwaResults as gr
    if local:
        with_qq_plots = False


    if sd:
        kwargs['snps'] = sd.getSnps()
        kwargs['positions'] = sd.getPositions()
        kwargs['chromosomes'] = sd.get_chr_list()
        d = sd.get_mafs()
        kwargs['macs'] = d['mafs']
        kwargs['mafs'] = d['marfs']
    if snp_priors:
        print 'Using provided SNP priors'
        kwargs['snp_priors'] = snp_priors[:]

    snps = kwargs['snps'][:]
    positions = kwargs['positions'][:]
    chromosomes = kwargs['chromosomes'][:]
    mafs = kwargs['mafs'][:]
    macs = kwargs['macs'][:]

    chr_pos_list = zip(chromosomes, positions)
    lmm = LinearMixedModel(phenotypes)
    lmm.add_random_effect(K)
    if K2 != None:
        lmm.add_random_effect(K2)
    num_snps = len(snps)

    if snp_priors == None:
        print "Using default SNP priors"
        snp_priors = [1.0 / num_snps] * num_snps

    if not sign_threshold:  # Then use Bonferroni threshold
        sign_threshold = 1.0 / (num_snps * 20.0)

    print "Running EMMAX stepwise"
    s1 = time.time()
    step_info_list = []
    cofactors = []  # A list of the loci found, together with their statistics.
    cofactor_snps = []
    step_i = 0
    num_par = 2  # mean and variance scalar
    num_pher_0 = 0

    if K2 != None:  # Then first estimate K
        res = lmm.get_estimates_3()
        pherit = res['perc_var1']
        print res['perc_var1'], res['perc_var2']
        K = res['opt_k']
    eig_L = lmm._get_eigen_L_(K=K)
    eig_R = lmm._get_eigen_R_(K=K)

    reml_res = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
    ml_res = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
    H_sqrt_inv = reml_res['H_sqrt_inv']
    ll = ml_res['max_ll']
    rss = float(reml_res['rss'])
    reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
    criterias = {'ebics':[], 'mbics':[], 'bonf':[], 'mbonf':[]}
    (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n)  # Calculate the BICs
    criterias['ebics'].append(extended_bic)
    criterias['mbics'].append(modified_bic)
    max_cofactor_pval = 0  # 5e-324 #min float, a hack to fix an annoying bug
    criterias['mbonf'].append(max_cofactor_pval)
    criterias['bonf'].append(0)

    criterias['min_cof_ppa'] = [1]  # Posterior probability of association
    cof_snp_priors = []
    ppa_cofactors = []

    action = 'None'
    if K2 == None:
        pherit = reml_res['pseudo_heritability']
    print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f, num_snps=%d' % \
        (step_i, action, num_par, pherit, ll, rss, reml_mahalanobis_rss, \
         bic, extended_bic, modified_bic, num_snps)
    print 'Cofactors:', _cofactors_to_string_(cofactors)
    quantiles_dict = {'log':[], 'norm':[], 'labels':[]}

    for step_i in range(1, num_steps + 1):
        emmax_res = lmm._emmax_f_test_(snps, H_sqrt_inv, snp_priors=snp_priors, emma_num=emma_num)
        if step_i == 1:
            first_emmax_res = emmax_res
        min_pval_i = sp.argmin(emmax_res['ps'])
        min_pval = emmax_res['ps'][min_pval_i]
        mahalnobis_rss = emmax_res['rss'][min_pval_i]
        min_pval_chr_pos = chr_pos_list[min_pval_i]
        max_ppa_i = sp.argmax(emmax_res['ppas'])
        max_ppa = emmax_res['ppas'][max_ppa_i]
        print 'Min p-value:', min_pval
        criterias['bonf'].append(min_pval)
        print 'Min Mahalanobis RSS:', mahalnobis_rss
        step_info = {'pseudo_heritability':pherit, 'rss':rss, \
            'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'mahalanobis_rss':mahalnobis_rss,
            'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic, 'mbonf':max_cofactor_pval,
            'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:],
            'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos,
            'max_ppa':max_ppa, 'max_ppa_pval':emmax_res['ps'][max_ppa_i],
            'max_ppa_chr_pos':chr_pos_list[max_ppa_i], 'ppa_cofactors':map(tuple, ppa_cofactors[:])}
        ppas = emmax_res['ppas'].tolist()

        if snp_choose_criteria == 'pval':
            snp_i = min_pval_i
        elif snp_choose_criteria == 'ppas':
            snp_i = max_ppa_i


        ex_pvals = emmax_res['ps'].tolist()
        ex_perc_var_expl = emmax_res['var_perc'].tolist()

        # Plot gwas results per step 
        if file_prefix:
            _plot_manhattan_and_qq_(file_prefix, step_i - 1, ex_pvals, quantiles_dict, positions=positions,
                    chromosomes=chromosomes, mafs=mafs, macs=macs, perc_var_expl=ex_perc_var_expl,
                    plot_bonferroni=True,
                    highlight_markers=cofactors, cand_genes=cand_gene_list, plot_xaxis=plot_xaxis,
                    log_qq_max_val=log_qq_max_val, with_qq_plots=with_qq_plots,
                    highlight_loci=highlight_loci, write_pvals=save_pvals, ppas=ppas,
                    highlight_ppa_markers=ppa_cofactors, markersize=markersize,
                    chrom_col_map=chrom_col_map)
        if save_pvals or pval_file_prefix:
            res = gr.Result(scores=ex_pvals, perc_var_expl=ex_perc_var_expl, **kwargs)
            res.filter_percentile(0.02, reversed=True)            
            if save_pvals:
                step_info['res'] = res
            if pval_file_prefix:
                pval_file_name = '%s_step%d.pvals' % (pval_file_prefix, step_i)
                res.write_to_file(pval_file_name, only_pickled=True, additional_columns='perc_var_expl')

        step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_pvals)
        step_info['pval_median'] = agr.calc_median(ex_pvals)
        print step_info['kolmogorov_smirnov'], step_info['pval_median']
        step_info_list.append(step_info)


        # Adding the new SNP as a cofactor
        lmm.add_factor(snps[snp_i])
        cofactor_snps.append(snps[snp_i])

        if K2 != None:  # Again first estimate K
            res = lmm.get_estimates_3()
            pherit = res['perc_var1']
            K = res['opt_k']
            eig_L = lmm._get_eigen_L_(K=K)


        eig_R = lmm._get_eigen_R_(X=lmm.X, K=K)
        reml_res = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
        ml_res = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
        H_sqrt_inv = reml_res['H_sqrt_inv']
        ll = ml_res['max_ll']
        rss = float(reml_res['rss'])
        reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
        num_par += 1
        action = '+'
        cof_snp_priors.append(snp_priors[snp_i])
        ppa_cofactors.append([chromosomes[snp_i], positions[snp_i], max_ppa])
        cofactors.append([chromosomes[snp_i], positions[snp_i], min_pval])


        # Re-estimate the p-value of the cofactors... with the smallest in the list.
        cofactor_pvals = []
        if snp_priors != None:
            cofactor_ppas = []  # Posterior probabilities of association
        for i, snp in enumerate(cofactor_snps):
            t_cofactors = cofactor_snps[:]
            del t_cofactors[i]
            lmm.set_factors(t_cofactors)
            res = lmm._emmax_f_test_([snp], H_sqrt_inv, snp_priors=[cof_snp_priors[i]], emma_num=0)
            cofactor_ppas.append(res['ppas'][0])
            ppa_cofactors[i][2] = res['ppas'][0]
            pval = res['ps'][0]
            cofactor_pvals.append(pval)
            cofactors[i][2] = -math.log10(pval)
        lmm.set_factors(cofactor_snps)
        max_cofactor_pval = max(cofactor_pvals)
        criterias['mbonf'].append(max_cofactor_pval)
        if snp_priors != None:
            criterias['min_cof_ppa'].append(min(cofactor_ppas))


        # Remove the found SNP from considered SNPs
        del snps[snp_i]
        del snp_priors[snp_i]
        del positions[snp_i]
        del chromosomes[snp_i]
        del chr_pos_list[snp_i]
        del mafs[snp_i]
        del macs[snp_i]
        num_snps -= 1


        (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n)  # Calculate the BICs
        criterias['ebics'].append(extended_bic)
        criterias['mbics'].append(modified_bic)

        if K2 == None:
            pherit = reml_res['pseudo_heritability']
        print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f, num_snps=%d' % \
            (step_i, action, num_par, pherit, ll, rss, reml_mahalanobis_rss, \
             bic, extended_bic, modified_bic, num_snps)

        print 'Cofactors:', _cofactors_to_string_(cofactors)
        print ppa_cofactors
#        if reml_res['pseudo_heritability'] < 0.01 and num_pher_0 < 1:
#            num_pher_0 += 1
#        elif reml_res['pseudo_heritability'] < 0.01:
        # if reml_res['pseudo_heritability'] < 0.01:
        if pherit < 0.001:
            if num_pher_0 < 1:
                num_pher_0 += 1
            else:
                print 'Breaking early, since pseudoheritability is close to 0.'
                break

    emmax_res = lmm._emmax_f_test_(snps, H_sqrt_inv, snp_priors=snp_priors, emma_num=emma_num)  #FINISH!!!
    min_pval_i = sp.argmin(emmax_res['ps'])
    min_pval = emmax_res['ps'][min_pval_i]
    mahalnobis_rss = emmax_res['rss'][min_pval_i]
    min_pval_chr_pos = chr_pos_list[min_pval_i]
    max_ppa_i = sp.argmax(emmax_res['ppas'])
    ppas = emmax_res['ppas'].tolist()
    print 'Min p-value:', min_pval
    print 'Min Mahalanobis RSS:', mahalnobis_rss
    step_info = {'pseudo_heritability':pherit, 'rss':rss, 'reml_mahalanobis_rss': reml_res['mahalanobis_rss'],
        'mahalanobis_rss':mahalnobis_rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
        'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:],
        'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos,
        'max_ppa':emmax_res['ppas'][max_ppa_i], 'max_ppa_pval':emmax_res['ps'][max_ppa_i],
        'max_ppa_chr_pos':chr_pos_list[max_ppa_i], 'ppa_cofactors':map(tuple, ppa_cofactors[:])}


    ex_pvals = emmax_res['ps'].tolist()
    ex_perc_var_expl = emmax_res['var_perc'].tolist()
    if save_pvals:
        step_info['ps'] = ex_pvals
    print "Generating plots"
    # Now plotting!
    if file_prefix:
        _plot_manhattan_and_qq_(file_prefix, step_i, ex_pvals, quantiles_dict, positions=positions,
                    chromosomes=chromosomes, mafs=mafs, macs=macs, perc_var_expl=ex_perc_var_expl,
                    plot_bonferroni=True,
                    highlight_markers=cofactors, cand_genes=cand_gene_list, plot_xaxis=plot_xaxis,
                    log_qq_max_val=log_qq_max_val, with_qq_plots=with_qq_plots,
                    highlight_loci=highlight_loci, write_pvals=save_pvals, ppas=ppas,
                    highlight_ppa_markers=ppa_cofactors, markersize=markersize,
                    chrom_col_map=chrom_col_map)
        # Plot posterior probabilities of association
    if save_pvals or pval_file_prefix:
        res = gr.Result(scores=ex_pvals, perc_var_expl=ex_perc_var_expl, **kwargs)
        res.filter_percentile(0.02, reversed=True)            
        if save_pvals:
            step_info['res'] = res
        if pval_file_prefix:
            pval_file_name = '%s_step%d.pvals' % (pval_file_prefix, step_i)
            res.write_to_file(pval_file_name, only_pickled=True, additional_columns='perc_var_expl')


    step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(ex_pvals)
    step_info['pval_median'] = agr.calc_median(ex_pvals)
    print step_info['kolmogorov_smirnov'], step_info['pval_median']
    step_info_list.append(step_info)

    if pval_file_prefix != None or save_pvals:
        res = gr.Result(scores=ex_pvals, perc_var_expl=ex_perc_var_expl, **kwargs)
        res.filter_percentile(0.02, reversed=True)
        if pval_file_prefix != None:
            pval_file_name = '%s_step%d.pvals' % (pval_file_prefix, step_i)
            res.write_to_file(pval_file_name, only_pickled=True, additional_columns='perc_var_expl')
        if save_pvals:
            step_info['res'] = res

    max_num_cofactors = len(cofactors)

    # Now backward stepwise.
    if forward_backwards:
        print 'Starting backwards..'
        while len(cofactor_snps) > 1:
            step_i += 1
            f_stats = sp.zeros(len(cofactor_snps))
            ppas = sp.zeros(len(cofactor_snps))
            for i, snp in enumerate(cofactor_snps):
                t_cofactors = cofactor_snps[:]
                del t_cofactors[i]
                lmm.set_factors(t_cofactors)
                res = lmm._emmax_f_test_([snp], H_sqrt_inv, snp_priors=[cof_snp_priors[i]], emma_num=0)
                ppas[i] = res['ppas'][0]
                cofactors[i][2] = -math.log10(res['ps'][0])
                f_stats[i] = res['f_stats'][0]
            if snp_choose_criteria == 'pval':
                i_to_remove = f_stats.argmin()
            elif snp_choose_criteria == 'ppas':
                i_to_remove = ppas.argmin()
            del ppa_cofactors[i_to_remove]
            del cofactor_snps[i_to_remove]
            del cofactors[i_to_remove]
            lmm.set_factors(cofactor_snps)
            num_snps += 1


            # Re-estimating the REML and ML.
            if K2 != None:  # Again first estimate K
                res = lmm.get_estimates_3()
                pherit = res['perc_var1']
                K = res['opt_k']
                eig_L = lmm._get_eigen_L_(K=K)
            eig_R = lmm._get_eigen_R_(X=lmm.X, K=K)
            reml_res = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
            ml_res = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
            ll = ml_res['max_ll']
            H_sqrt_inv = reml_res['H_sqrt_inv']
            rss = float(reml_res['rss'])
            reml_mahalanobis_rss = float(reml_res['mahalanobis_rss'])
            num_par -= 1
            action = '-'


            # Update the p-values
            cofactor_pvals = []
            if snp_priors != None:
                cofactor_ppas = []  # Posterior probabilities of association
            for i, snp in enumerate(cofactor_snps):
                t_cofactors = cofactor_snps[:]
                del t_cofactors[i]
                lmm.set_factors(t_cofactors)
                res = lmm._emmax_f_test_([snp], H_sqrt_inv, snp_priors=[cof_snp_priors[i]], emma_num=0)
                cofactor_ppas.append(res['ppas'][0])
                ppa_cofactors[i][2] = res['ppas'][0]
                pval = res['ps'][0]
                cofactor_pvals.append(pval)
                cofactors[i][2] = -math.log10(pval)
            max_cofactor_pval = max(cofactor_pvals)
            criterias['mbonf'].append(max_cofactor_pval)
            criterias['min_cof_ppa'].append(min(cofactor_ppas))

            # Calculate the BICs
            (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lmm.n)
            criterias['ebics'].append(extended_bic)
            criterias['mbics'].append(modified_bic)
            if K2 == None:
                pherit = reml_res['pseudo_heritability']
            print '\nStep %d: action=%s, num_par=%d, p_her=%0.4f, ll=%0.2f, rss=%0.2f, reml_m_rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f, num_snps=%d' % \
                (step_i, action, num_par, pherit, ll, rss,
                reml_mahalanobis_rss, bic, extended_bic, modified_bic, num_snps)

            print 'Cofactors:', _cofactors_to_string_(cofactors)
            print ppa_cofactors

            step_info = {'pseudo_heritability':pherit, 'rss':rss, \
                'reml_mahalanobis_rss': reml_res['mahalanobis_rss'], 'll':ll, 'bic':bic,
                'e_bic':extended_bic, 'm_bic':modified_bic, 'mbonf':max_cofactor_pval,
                'cofactors':map(tuple, cofactors[:]), 'cofactor_snps':cofactor_snps[:],
                'mahalanobis_rss':None, 'min_pval':None, 'min_pval_chr_pos':None,
                'kolmogorov_smirnov':None, 'pval_median':None,
                'ppa_cofactors': map(tuple, ppa_cofactors[:])}
            step_info_list.append(step_info)
            print step_info['kolmogorov_smirnov'], step_info['pval_median']

    opt_dict, opt_indices = _analyze_opt_criterias_(criterias, sign_threshold, max_num_cofactors, file_prefix,
                        with_qq_plots, lmm, step_info_list, quantiles_dict,
                        plot_bonferroni=True, cand_genes=cand_gene_list, plot_xaxis=plot_xaxis,
                        log_qq_max_val=log_qq_max_val, eig_L=eig_L, type='emmax',
                        highlight_loci=highlight_loci, write_pvals=save_pvals,
                        markersize=markersize, chrom_col_map=chrom_col_map, emma_num=emma_num,
                        **kwargs)

    for step_i in opt_indices:
        for h in ['mahalanobis_rss', 'min_pval', 'min_pval_chr_pos', 'kolmogorov_smirnov', 'pval_median', 'res']:
            step_info_list[step_i][h] = opt_indices[step_i][h]
            if save_pvals:
                step_info_list[step_i]['res'] = opt_indices[step_i]['res']

    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)

    if file_prefix:
        _plot_stepwise_stats_(file_prefix, step_info_list, sign_threshold, type='emmax')

    res_dict = {'step_info_list':step_info_list, 'first_emmax_res':first_emmax_res, 'opt_dict':opt_dict}

    return res_dict





def lm_step_wise(phenotypes, sd=None, num_steps=10, file_prefix=None, forward_backwards=True, local=False,
        cand_gene_list=None, plot_xaxis=True, with_qq_plots=True, sign_threshold=None, log_qq_max_val=5,
        highlight_loci=None, save_pvals=False, markersize=3, chrom_col_map=None, **kwargs):
    """
    Run simple step-wise linear model forward-backward.
    """
    # import plotResults as pr
    if local:
        with_qq_plots = False

    if sd:
        kwargs['snps'] = sd.getSnps()
        kwargs['positions'] = sd.getPositions()
        kwargs['chromosomes'] = sd.get_chr_list()
        d = sd.get_mafs()
        kwargs['macs'] = d['mafs']
        kwargs['mafs'] = d['marfs']

    snps = kwargs['snps'][:]
    positions = kwargs['positions'][:]
    chromosomes = kwargs['chromosomes'][:]
    mafs = kwargs['mafs'][:]
    macs = kwargs['macs'][:]

    chr_pos_list = zip(chromosomes, positions)
    lm = LinearModel(phenotypes)
    num_snps = len(snps)

    if not sign_threshold:  # Then use Bonferroni threshold
        sign_threshold = 1.0 / (num_snps * 20.0)

    print "Running step-wise LM"
    s1 = time.time()
    step_info_list = []
    cofactors = []  # A list of the loci found, together with their statistics.
    cofactor_snps = []
    step_i = 0
    num_par = 2  # mean and variance scalar

    rss = lm.get_rss()
    ll = lm.get_ll(rss)
    criterias = {'ebics':[], 'mbics':[], 'bonf':[], 'mbonf':[]}
    (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n)  # Calculate the BICs
    criterias['ebics'].append(extended_bic)
    criterias['mbics'].append(modified_bic)
    max_cofactor_pval = 0
    criterias['mbonf'].append(max_cofactor_pval)
    criterias['bonf'].append(0)
    action = 'None'
    print '\nStep %d: action=%s, num_par=%d, ll=%0.2f, rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
        (step_i, action, num_par, ll, rss, bic, extended_bic, modified_bic)
    print 'Cofactors:', _cofactors_to_string_(cofactors)
    quantiles_dict = {'log':[], 'norm':[], 'labels':[]}

    for step_i in range(1, num_steps + 1):
        lm_res = lm.fast_f_test(snps)
        if step_i == 1:
            first_lm_res = lm_res
        min_pval_i = sp.argmin(lm_res['ps'])
        min_pval = lm_res['ps'][min_pval_i]
        min_pval_chr_pos = chr_pos_list[min_pval_i]
        print 'Min p-value:', min_pval
        criterias['bonf'].append(min_pval)
        step_info = {'rss':rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
                'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]),
                'cofactor_snps':cofactor_snps[:], 'min_pval':min_pval,
                'min_pval_chr_pos': min_pval_chr_pos}

        lm_pvals = lm_res['ps'].tolist()
        # Plot gwas results per step 
        if file_prefix:
            _plot_manhattan_and_qq_(file_prefix, step_i - 1, lm_pvals, quantiles_dict, positions=positions,
                chromosomes=chromosomes, mafs=mafs, macs=macs, plot_bonferroni=True, highlight_markers=cofactors,
                cand_genes=cand_gene_list, plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val,
                with_qq_plots=with_qq_plots, highlight_loci=highlight_loci, write_pvals=save_pvals,
                markersize=markersize, chrom_col_map=chrom_col_map)
        if save_pvals:
            step_info['ps'] = lm_pvals


        if cand_gene_list:
            # Calculate candidate gene enrichments.
            pass
        step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(lm_pvals)
        step_info['pval_median'] = agr.calc_median(lm_pvals)
        print step_info['kolmogorov_smirnov'], step_info['pval_median']
        step_info_list.append(step_info)


        # Adding the new SNP as a cofactor
        lm.add_factor(snps[min_pval_i])
        cofactor_snps.append(snps[min_pval_i])
        rss = lm.get_rss()
        ll = lm.get_ll(rss)
        num_par += 1
        action = '+'

        cofactors.append([min_pval_chr_pos[0], min_pval_chr_pos[1], min_pval])


        # Re-estimate the p-value of the cofactors... with the smallest in the list.
        cofactor_pvals = []
        for i, snp in enumerate(cofactor_snps):
            t_cofactors = cofactor_snps[:]
            del t_cofactors[i]
            lm.set_factors(t_cofactors)
            pval = lm.fast_f_test([snp])['ps'][0]
            cofactor_pvals.append(pval)
            cofactors[i][2] = -math.log10(pval)
        lm.set_factors(cofactor_snps)
        max_cofactor_pval = max(cofactor_pvals)
        criterias['mbonf'].append(max_cofactor_pval)

        # Remove the found SNP from considered SNPs
        del snps[min_pval_i]
        del positions[min_pval_i]
        del chromosomes[min_pval_i]
        del chr_pos_list[min_pval_i]
        del mafs[min_pval_i]
        del macs[min_pval_i]
        num_snps -= 1


        (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n)  # Calculate the BICs
        criterias['ebics'].append(extended_bic)
        criterias['mbics'].append(modified_bic)

        print '\nStep %d: action=%s, num_par=%d, ll=%0.2f, rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
            (step_i, action, num_par, ll, rss, bic, extended_bic, modified_bic)
        print 'Cofactors:', _cofactors_to_string_(cofactors)

    lm_res = lm.fast_f_test(snps)
    min_pval_i = sp.argmin(lm_res['ps'])
    min_pval = lm_res['ps'][min_pval_i]
    min_pval_chr_pos = chr_pos_list[min_pval_i]
    print 'Min p-value:', min_pval
    step_info = {'rss':rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
        'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]),
        'cofactor_snps':cofactor_snps[:], 'min_pval':min_pval, 'min_pval_chr_pos': min_pval_chr_pos}
    lm_pvals = lm_res['ps'].tolist()
    if save_pvals:
        step_info['ps'] = lm_pvals

    # Now plotting!
    print "Generating plots"
    if file_prefix:
        _plot_manhattan_and_qq_(file_prefix, step_i, lm_pvals, quantiles_dict, positions=positions,
                    chromosomes=chromosomes, mafs=mafs, macs=macs, plot_bonferroni=True, highlight_markers=cofactors,
                    cand_genes=cand_gene_list, plot_xaxis=plot_xaxis, log_qq_max_val=log_qq_max_val,
                    with_qq_plots=with_qq_plots, highlight_loci=highlight_loci, write_pvals=save_pvals,
                    markersize=markersize, chrom_col_map=chrom_col_map)

    max_num_cofactors = len(cofactors)
    step_info['kolmogorov_smirnov'] = agr.calc_ks_stats(lm_pvals)
    step_info['pval_median'] = agr.calc_median(lm_pvals)
    print step_info['kolmogorov_smirnov'], step_info['pval_median']
    step_info_list.append(step_info)

    # Now backward stepwise.
    if forward_backwards:
        print 'Starting backwards..'
        while len(cofactor_snps) > 1:
            step_i += 1
            f_stats = sp.zeros(len(cofactor_snps))
            for i, snp in enumerate(cofactor_snps):
                t_cofactors = cofactor_snps[:]
                del t_cofactors[i]
                lm.set_factors(t_cofactors)
                res = lm.fast_f_test([snp])
                cofactors[i][2] = -math.log10(res['ps'][0])
                f_stats[i] = res['f_stats'][0]
            i_to_remove = f_stats.argmin()
            del cofactor_snps[i_to_remove]
            del cofactors[i_to_remove]
            lm.set_factors(cofactor_snps)
            num_snps += 1

            # Re-estimating the REML and ML.
            rss = lm.get_rss()
            ll = lm.get_ll(rss)
            num_par -= 1
            action = '-'

            # Update the p-values
            cofactor_pvals = []
            for i, snp in enumerate(cofactor_snps):
                t_cofactors = cofactor_snps[:]
                del t_cofactors[i]
                lm.set_factors(t_cofactors)
                res = lm.fast_f_test([snp])
                pval = res['ps'][0]
                cofactor_pvals.append(pval)
                cofactors[i][2] = -math.log10(pval)
            max_cofactor_pval = max(cofactor_pvals)
            criterias['mbonf'].append(max_cofactor_pval)

            # Calculate the BICs
            (bic, extended_bic, modified_bic) = _calc_bic_(ll, num_snps, num_par, lm.n)
            criterias['ebics'].append(extended_bic)
            criterias['mbics'].append(modified_bic)
            print '\nStep %d: action=%s, num_par=%d, ll=%0.2f, rss=%0.2f, bic=%0.2f, extended_bic=%0.2f, modified_bic=%0.2f' % \
                (step_i, action, num_par, ll, rss, bic, extended_bic, modified_bic)
            print 'Cofactors:', _cofactors_to_string_(cofactors)

            step_info = {'rss':rss, 'll':ll, 'bic':bic, 'e_bic':extended_bic, 'm_bic':modified_bic,
                'mbonf':max_cofactor_pval, 'cofactors':map(tuple, cofactors[:]),
                'cofactor_snps':cofactor_snps[:], 'min_pval':None, 'min_pval_chr_pos':None,
                'kolmogorov_smirnov':None, 'pval_median':None}
            step_info_list.append(step_info)
            print cofactors

    opt_dict, opt_indices = _analyze_opt_criterias_(criterias, sign_threshold, max_num_cofactors, file_prefix,
                        with_qq_plots, lm, step_info_list, quantiles_dict,
                        plot_bonferroni=True, cand_genes=cand_gene_list, plot_xaxis=plot_xaxis,
                        log_qq_max_val=log_qq_max_val, type='lm', highlight_loci=highlight_loci,
                        write_pvals=save_pvals, markersize=markersize,
                        chrom_col_map=chrom_col_map, **kwargs)


    for step_i in opt_indices:
        for h in ['min_pval', 'min_pval_chr_pos', 'kolmogorov_smirnov', 'pval_median']:
            step_info_list[step_i][h] = opt_indices[step_i][h]

    if file_prefix:
        _plot_stepwise_stats_(file_prefix, step_info_list, sign_threshold, type == 'lm')

    res_dict = {'step_info_list':step_info_list, 'first_lm_res':first_lm_res, 'opt_dict':opt_dict}

    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    return res_dict




def linear_model(snps, phenotypes, cofactors=None):
    lm = LinearModel(phenotypes)
    if cofactors:
        for cofactor in cofactors:
            lm.add_factor(cofactor)
    log.info("Running a standard linear model")
    t = Timer()
    res = lm.fast_f_test(snps)
    log.info('Took: %s' % t.stop(True))
    return res


def kruskal_wallis(genotype, phenVals, useTieCorrection=True):
    t = Timer()
    num_accessions=len(genotype.accessions)
    num_snps = genotype.num_snps
    assert num_accessions== len(phenVals), "SNPs and phenotypes are not equal length."
    ranks, group_counts = _kw_get_ranks_(phenVals)
    assert len(group_counts) == len(set(ranks)), 'Somethings wrong..'
    tieCorrection = 1.0
    if useTieCorrection:
        n_total = len(ranks)
        ones_array = np.repeat(1.0, len(group_counts))
        s = np.sum(group_counts * (group_counts * group_counts - ones_array))
        s = s / (n_total * (n_total * n_total - 1.0))
        tieCorrection = 1.0 / (1 - s)

    n = len(phenVals)
    ds = np.zeros(num_snps)
    c = 12.0 / (n * (n + 1))
    snps = genotype.get_snps_iterator()
    for i, snp in enumerate(snps):
        ns = np.bincount(snp)
        rs = np.array([0.0, np.sum(snp * ranks)])
        rs[0] = np.sum(ranks) - rs[1]
        nominator = 0.0
        mean_r = (n + 1) / 2.0
        for j in range(2):
            v = (rs[j] / ns[j]) - mean_r
            nominator += ns[j] * v * v
        ds[i] = (c * nominator) * tieCorrection

    ps = sp.stats.chi2.sf(ds, 1)

    log.info('Took %s' % t.stop(True))
    return {"ps":ps, "ds":ds}

def _kw_get_ranks_(values, withTies=True):
    """
    returns ranks (large values w. large ranks)
    """
    srt_vals = zip(values[:], range(len(values)))
    srt_vals.sort()
    srt_ranks = np.zeros(len(srt_vals))
    group_counts = []
    if withTies:
        i = 0
        while i < len(srt_vals):
            curVal = srt_vals[i][0]
            j = 1
            while i + j < len(srt_vals) and srt_vals[i + j][0] == curVal:
                j += 1
            rank = i + (j + 1) / 2.0
            max_i = i + j
            group_counts.append(j)
            while i < max_i:
                srt_ranks[i] = rank
                i += 1

    else:
        srt_ranks = np.arange(len(srt_vals))
        group_counts = [1] * len(srt_vals)

    ranks = np.zeros(len(srt_vals))

    for i, (val, val_index) in enumerate(srt_vals):
        ranks[val_index] = srt_ranks[i]
    return ranks, np.array(group_counts)

