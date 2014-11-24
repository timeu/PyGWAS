#
# class MultipleTraitLMM(LinearMixedModel):
#    """
#    A class to encompass multiple traits/environments mixed model.
#    """
#    def __init__(self, joint_phen_vals=None, phenotype_list=None, joint_ets=None, ecotype_list=None,
#            unique_ets=None, K=None, trait_vec=None, Z=None, ecotype_maps=None, phen_corr_mat=None,
#            p_herits=None):
#        """
#        The fixed effects should be a list of fixed effect lists (SNPs)
#        """
#        if joint_phen_vals == None: # Joint phenotype vector.
#            joint_phen_vals = []
#            if phenotype_list != None:
#                for phen in phenotype_list:
#                    joint_phen_vals.extend(phen)
#        self.n = len(joint_phen_vals)
#        self.Y = sp.matrix(joint_phen_vals).T
#
#        if joint_ets == None: # Joint ecotype vector.
#            joint_ets = []
#            if ecotype_list != None:
#                for ets in ecotype_list:
#                    joint_ets.append(ets)
#                #self.et_list = ecotype_list
#
#        self.joint_ets = joint_ets
#        if unique_ets == None:
#            unique_ets = list(set(self.joint_ets))
#        self.unique_ets = unique_ets
#
#        self.ys = phenotype_list
#        if ecotype_list:
#            self.num_ets = map(len, ecotype_list)
#            phen_corr_mat = self.get_phen_corr(phenotype_list, ecotype_list)
#        if phen_corr_mat != None:
#            self.phen_corr_mat = phen_corr_mat
#        self.X = sp.mat(sp.repeat(1, self.n)).T #The intercept
#        self.num_unique = len(unique_ets)
#
#        if K != None:
#            self.K = kinship.scale_k(K)
#
#        #E in the model.
#        if trait_vec == None:
#            trait_vec = []
#            shift = self.num_ets[0]
#            for num_ets in self.num_ets[1:]: #Skip the first one..
#                trait_vec.append([0] * shift + [1] * num_ets + [0] * (self.n - num_ets - shift))
#                shift += num_ets
#            self.E = sp.mat(E).T
#        self.trait_vec = trait_vec
#
#        if p_herits != None:
#            self.p_herits = p_herits
#
#
#
#
#    def get_estimates_5(self, xs=None, ngrids=[5, 5, 5, 5, 5], llim= -3, ulim=3, method='REML',
#                verbose=False, dtype='single', debug=False):
#        #1. Get marginal estimates, use them as guides.  Assumes that the trait has been standardized.
#        sqrt_h = sp.sqrt(self.p_herits[0] * self.p_herits[1])
#        pcorr = self.phen_corr_mat[1, 0]
#        if pcorr == 0.0:
#            pcorr = 0.001
#        elif pcorr == 1.0:
#            pcorr = 0.999
#        if sqrt_h != 0:
#            start_gen_corr = pcorr / sqrt_h
#            if start_gen_corr > 0.999:
#                start_gen_corr = 0.999
#        else:
#            start_gen_corr = 0.999
#        start_g0g1_ratio = sp.log(self.p_herits[0] / self.p_herits[1])
#        start_g0e0_ratio = sp.log(self.p_herits[0] / (1.0 - self.p_herits[0]))
#        start_g1e1_ratio = sp.log(self.p_herits[1] / (1.0 - self.p_herits[1]))
#
#        #2. Set up grid, and estimate parameters.  3 ratios and one corr parameter to search over.
#        #for g0g1_ratio in
#
#
#        #3. Calculate the log-likelihoods, differentiate them, etc.
#        pass
#        # return LinearMixedModel.get_estimates_3(self, xs=xs, ngrids=ngrids, llim=llim, ulim=ulim, method=method, verbose=verbose, dtype=dtype, debug=debug)
#
#
#    def get_variance_matrix(self):
#        """
#        Assumes kinship is there...
#        """
#        gen_var_list = []
#        err_var_list = []
#        her_list = []
#        r_list = self.get_marginal_estimates()
#        for res in r_list:
#            her_list.append(res['pseudo_heritability'])
#            err_var_list.append(res['ve'])
#            gen_var_list.append(res['vg'])
#        print err_var_list, gen_var_list, her_list
#        V = sp.zeros((self.n, self.n))
#        m_i = 0
#        for i, ets_map1 in enumerate(self.ecotype_maps):
#            num_ets1 = len(ets_map1)
#            m_j = 0
#            for j, ets_map2 in enumerate(self.ecotype_maps):
#                num_ets2 = len(ets_map2)
#                if i == j:
#                    var_sclice = gen_var_list[i] * self.K[ets_map1, :][:, ets_map2] \
#                            + err_var_list[i] * sp.eye(num_ets1)
#                else:
#                    if  her_list[i] > 0 and her_list[j] > 0:
#                        rho = self.corr_mat[i, j] / sp.sqrt(her_list[i] * her_list[j])
#                        if rho > 1.0: rho = 1.0
#                        if rho < -1.0: rho = -1.0
#                    else:
#                        rho = 0
#                    print i, j, rho
#                    var_sclice = rho * math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
#                            self.K[ets_map1, :][:, ets_map2]
#                V[m_i:m_i + num_ets1, m_j:m_j + num_ets2] = var_sclice
#                m_j += num_ets2
#            m_i += num_ets1
#        H_sqrt = cholesky(V)
#        H_sqrt_inv = sp.mat(H_sqrt.T).I
#        return {'V':V, 'H_sqrt':H_sqrt, 'H_sqrt_inv':H_sqrt_inv}
#
#
#
#
#    def get_phen_corr(self, phenotype_list, ecotype_list):
#        print 'Calculating phenotype correlation'
#        num_traits = len(phenotype_list)
#        corr_mat = sp.ones((num_traits, num_traits))
#        for i, (pvs1, ets1) in enumerate(zip(phenotype_list, ecotype_list)):
#            for j in range(i):
#                ets2 = ecotype_list[j]
#                pvs2 = phenotype_list[j]
#                common_ets = set(ets1).intersection(set(ets2))
#                ets_ix1 = map(ets1.index, common_ets)
#                ets_ix2 = map(ets2.index, common_ets)
#                vs1 = [pvs1[et_i] for et_i in ets_ix1]
#                vs2 = [pvs2[et_i] for et_i in ets_ix2]
#                corr_mat[i, j] = sp.corrcoef(vs1, vs2)[0, 1]
#                corr_mat[j, i] = corr_mat[i, j]
#        print corr_mat
#        return corr_mat
#
#
#
#    def emmax_f_test_fm(self, snps, with_betas=False, type='common'):
#        """
#        EMMAX Full model f-test 
#        returns three sets of p-values
#        
#        """
#        s1 = time.time()
#        print 'Getting variance estimates'
#        res = self.get_variance_matrix()
#        print 'Done.'
#
#        X = sp.hstack([self.X, self.E])
#        print 'Perfoming GLS scan'
#        r = self._emmax_f_test_(snps, X, res['H_sqrt_inv'], E=self.E, with_betas=with_betas)
#        secs = time.time() - s1
#        if secs > 60:
#            mins = int(secs) / 60
#            secs = secs - mins * 60
#            print 'Took %d mins and %f seconds to run EMMAX on the traits...' % (mins, secs)
#        else:
#            print 'Took %f seconds to run EMMAX on the traits..' % (secs)
#        return r
#
#
#
#    def emmax_f_test(self, snps, with_betas=False, type='common'):
#        """
#        EMMAX implementation
#        Single SNPs, multiple traits.
#        
#        """
#        s1 = time.time()
#        print 'Getting variance estimates'
#        res = self.get_variance_matrix()
#        print 'Done.'
#
#        X = sp.hstack([self.X, self.E])
#        print 'Perfoming GLS scan'
#        r = self._emmax_f_test_(snps, X, res['H_sqrt_inv'], with_betas=with_betas)
#        secs = time.time() - s1
#        if secs > 60:
#            mins = int(secs) / 60
#            secs = secs - mins * 60
#            print 'Took %d mins and %f seconds to run EMMAX on the traits...' % (mins, secs)
#        else:
#            print 'Took %f seconds to run EMMAX on the traits..' % (secs)
#        return r
#
#
#
#    def _emmax_f_test_(self, snps, X, H_sqrt_inv, E=None, verbose=True, with_betas=False, return_transformed_snps=False):
#        """
#        EMMAX implementation (in python)
#        Single SNPs, multiple traits
#        
#        Methods:
#            normal - Normal regression
#            qr - Uses QR decomposition to speed up regression with many co-factors.
#            
#        If E (environment) is provided then GxE factors are considered as well.
#            
#        """
#        dtype = 'single'
#        q = 1  # Single SNP is being tested
#        p = len(X.T) + q
#        n = self.n
#        n_p = n - p
#        num_snps = len(snps)
#
#        h0_X = sp.mat(H_sqrt_inv * X, dtype=dtype)
#        Y = H_sqrt_inv * self.Y    #The transformed outputs.
#        (h0_betas, h0_rss, h0_rank, h0_s) = linalg.lstsq(h0_X, Y)
#        Y = sp.mat(Y - h0_X * h0_betas, dtype=dtype)
#        h0_betas = map(float, list(h0_betas))
#
#        #H_sqrt_inv = H_sqrt_inv * self.Z
#        (Q, R) = qr_decomp(h0_X)  #Do the QR-decomposition for the Gram-Schmidt process.
#        Q = sp.mat(Q)
#        Q2 = Q * Q.T
#        M_0 = sp.mat(H_sqrt_inv.T * (sp.eye(n) - Q2), dtype=dtype)
#        M = self.Z.T * M_0
#
#        if E != None:
#            q_e = q + E.shape[1]
#            p_e = p + E.shape[1]
#            n_p_e = n - p_e
#            M_Es = []  #A list of matrices to genereate each GxE interaction
#            for i in range(E.shape[1]):
#                M_Es.append(sp.mat(self.Z * sp.array(E[:, i:i + 1]), dtype=dtype).T * M_0)
#            rss_full_list = sp.repeat(h0_rss, num_snps)
#
#
#        rss_list = sp.repeat(h0_rss, num_snps)
#        chunk_size = len(Y)
#        for i in range(0, num_snps, chunk_size): #Do the dot-product in chuncks!
#            snps_chunk = sp.matrix(snps[i:i + chunk_size], dtype=dtype)
#            Xs = snps_chunk * M
#            if E != None:
#                X_Es = []
#                for M_E in M_Es:
#                    X_Es.append(snps_chunk * M_E)
#            for j, X_j in enumerate(Xs):
#                (betas, rss, p, sigma) = linalg.lstsq(X_j.T, Y)
#                if rss:
#                    rss_list[i + j] = rss[0]
#                    if E != None:
#                        l = [X_E[j] for X_E in X_Es]
#                        l.append(X_j)
#                        X_mat = sp.mat(sp.vstack(l), dtype=dtype).T
#                        (betas, rss_full, p, sigma) = linalg.lstsq(X_mat, Y, overwrite_a=True)
#                        if rss_full and (rss - rss_full) / h0_rss > -0.0001:
#                            rss_full_list[i + j] = rss_full[0] if rss_full < rss else rss[0]
#                        elif rss_full and (rss - rss_full) / h0_rss <= -0.0001:
#                            print "Funny business"
#                            pdb.set_trace()
#
#                if verbose and num_snps >= 10 and (i + j + 1) % (num_snps / 10) == 0: #Print dots
#                    sys.stdout.write('.')
#                    sys.stdout.flush()
#
#        if verbose and num_snps >= 10:
#            sys.stdout.write('\n')
#        rss_ratio = h0_rss / rss_list
#        var_perc = 1 - 1 / rss_ratio
#        f_stats = (rss_ratio - 1) * n_p / float(q)
#        p_vals = stats.f.sf(f_stats, q, n_p)
#
#        res_d = {'ps':p_vals, 'f_stats':f_stats, 'rss':rss_list, 'var_perc':var_perc,
#            'h0_rss':h0_rss, 'h0_betas':h0_betas}
#
#        if E != None:
#            rss_ratio = h0_rss / rss_full_list
#            res_d['var_perc_full'] = 1 - 1 / rss_ratio
#            f_stats_full = (rss_ratio - 1) * n_p_e / float(q_e)
#            res_d['ps_full'] = stats.f.sf(f_stats_full, q_e, n_p_e)
#            res_d['f_stats_full'] = f_stats_full
#
#
#            rss_ratio = rss_list / rss_full_list
#            res_d['var_perc_ge'] = 1 - 1 / rss_ratio
#            f_stats_ge = (rss_ratio - 1) * (n_p_e + q) / float(q_e - q)
#            res_d['ps_ge'] = stats.f.sf(f_stats_ge, q_e - q, n_p_e + q)
#            res_d['f_stats_ge'] = f_stats_ge
#        return res_d



# def _test_joint_analysis_(run_id='FT_suzi_log'):
#    import dataParsers as dp
#    import phenotypeData as pd
#    import env
#    #filename = "/Users/bjarnivilhjalmsson/Projects/Data/phenotypes/seed_size.csv"
#    filename = env.env['phen_dir'] + 'phen_raw_112210.csv'
#    phed = pd.parse_phenotype_file(filename)
#    phed.convert_to_averages()
#
#    mac_threshold = 15
#    #pids = [187, 188, 189, 190]
#    pids = [5, 6, 7]
#    #pids = [164, 165, 166]
#    #pids = [1, 2, 3, 4]
#
#    joint_phen_vals = []
#    joint_ecotypes = []
#    gen_var_list = []
#    err_var_list = []
#    her_list = []
#    k_list = []
#    for pid in pids:
#        print phed.get_name(pid)
#        sd = dp.load_250K_snps(debug_filter=0.05)
#        sd.coordinate_w_phenotype_data(phed, pid)
#        phed.normalize_values(pid)
#        phenotypes = phed.get_values(pid)
#        phed.log_transform(pid)
#        #Create a new phenotype.
#        joint_phen_vals.extend(phenotypes)
#        ecotypes = phed.get_ecotypes(pid)
#        #Figure out all ecotypes which are in 250K and for which we have phenotypic values.
#        joint_ecotypes.extend(ecotypes)
#        K = load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', ecotypes)
#        #Get pseudoheritabilities for all phenotypes
#        lmm = LinearMixedModel(phenotypes)
#        lmm.add_random_effect(K)
#        res = lmm.get_REML()
#        print 'Pseudo-heritatbility:', res['pseudo_heritability']
#        her_list.append(res['pseudo_heritability'])
#        err_var_list.append(res['ve'])
#        gen_var_list.append(res['vg'])
#
#
#    #E in the model.
#    num_vals = len(joint_phen_vals)
#    joint_X = [[1] * num_vals]
#    shift = phed.num_vals(pids[0])
#    for pid in pids[1:]:
#        n = phed.num_vals(pid)
#        joint_X.append([0] * shift + [1] * n + [0] * (num_vals - n - shift))
#        shift += n
#
#
#    #Get correlations between traits..
#    corr_mat, pids = phed.get_correlations(pids)
#    #corr_mat = corr_mat * corr_mat
#    print gen_var_list, err_var_list, corr_mat, pids
#
#    unique_ecotypes = list(set(joint_ecotypes))
#    sd = dp.load_250K_snps()
#    sd.filter_accessions(unique_ecotypes)
#    if mac_threshold:
#        sd.filter_mac_snps(mac_threshold) #Filter MAF SNPs!
#    #unique_ecotypes = sd.accessions
#    K = load_kinship_from_file(env.env['data_dir'] + 'kinship_matrix_cm72.pickled', unique_ecotypes)
#
#    #Create ecotype maps
#    ecotype_maps = {}
#    for pid in pids:
#        l = []
#        ets = phed.get_ecotypes(pid)
#        for et in ets:
#            l.append(unique_ecotypes.index(et))
#        ecotype_maps[pid] = l
#
#    #Construct new variance matrix
#    V = sp.zeros((len(joint_ecotypes), len(joint_ecotypes)))
#    m_i = 0
#    for i, pid1 in enumerate(pids):
#        ets_map1 = ecotype_maps[pid1]
#        num_ets1 = len(ets_map1)
#        m_j = 0
#        for j, pid2 in enumerate(pids):
#            ets_map2 = ecotype_maps[pid2]
#            num_ets2 = len(ets_map2)
#            if i == j:
#                var_sclice = math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
#                        K[ets_map1, :][:, ets_map2] + err_var_list[i] * sp.eye(num_ets1)
#            else:
#                if  her_list[i] > 0 and her_list[j] > 0:
#                    rho = corr_mat[i, j] / sp.sqrt(her_list[i] * her_list[j])
#                    if rho > 1.0: rho = 1.0
#                    if rho < -1.0: rho = 1.0
#                else:
#                    rho = 0
#                print i, j, rho
#                var_sclice = rho * math.sqrt(gen_var_list[i] * gen_var_list[j]) * \
#                        K[ets_map1, :][:, ets_map2]
#            V[m_i:m_i + num_ets1, m_j:m_j + num_ets2] = var_sclice
#            m_j += num_ets2
#        m_i += num_ets1
#
#    print'Performing Cholesky decomposition'
#    H_sqrt = cholesky(V)
#    H_sqrt_inv = sp.mat(H_sqrt.T).I
#
#    #pdb.set_trace()
#
#    #5. Run analysis..
#    lmm = LinearMixedModel(joint_phen_vals)
#    lmm.set_factors(joint_X, False)
#    #print lmm.X
#    #Need to fix the SNPs!!!
#    Z = sp.int16(sp.mat(joint_ecotypes).T == sp.mat(unique_ecotypes))
#    snps = sd.getSnps()
#     positions = sd.getPositions()
#     chromosomes = sd.get_chr_list()
#    #pdb.set_trace()
#    res = lmm.emmax_full_model_gxe(snps, joint_X[1:], H_sqrt_inv, Z=Z)
#    #pdb.set_trace()
#    import gwaResults as gr
#    res_h01 = gr.Result(scores=res['ps_h01'].tolist(), positions=positions, chromosomes=chromosomes)
#    res_h02 = gr.Result(scores=res['ps_h02'].tolist(), positions=positions, chromosomes=chromosomes)
#    res_h12 = gr.Result(scores=res['ps_h12'].tolist(), positions=positions, chromosomes=chromosomes)
#    agr.qq_plot({'EMMAX_G':res_h01, 'EMMAX_FULL':res_h02, 'EMMAX_GxE':res_h12}, 1000, method_types=['emma', 'emma', 'emma'],
#        mapping_labels=['EMMAX_G', 'EMMAX_FULL', 'EMMAX_GxE'], phenName='joint', pngFile=env.env['results_dir'] + 'qq_plot_%s.png' % run_id)
#    agr.log_qq_plot({'EMMAX_G':res_h01, 'EMMAX_FULL':res_h02, 'EMMAX_GxE':res_h12}, 1000, 7, method_types=['emma', 'emma', 'emma'],
#        mapping_labels=['EMMAX_G', 'EMMAX_FULL', 'EMMAX_GxE'], phenName='joint', pngFile=env.env['results_dir'] + 'log_qq_plot%s.png' % run_id)
#    png_file_name = env.env['results_dir'] + 'manhattan_h01%s.png' % run_id
#    res_h01.neg_log_trans()
#    res_h01.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
#    png_file_name = env.env['results_dir'] + 'manhattan_h02%s.png' % run_id
#    res_h02.neg_log_trans()
#    res_h02.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)
#    png_file_name = env.env['results_dir'] + 'manhattan_h12%s.png' % run_id
#    res_h12.neg_log_trans()
#    res_h12.plot_manhattan(png_file=png_file_name, plot_bonferroni=True)



def get_emma_reml_estimates(y, K, K2=None, cofactors=None, include_intercept=True):
    lmm = LinearMixedModel(y)
    lmm.add_random_effect(K)
    if cofactors != None:
        lmm.set_factors(cofactors, include_intercept=include_intercept)
    if K2 != None:
        lmm.add_random_effect(K2)
        r = lmm.get_estimates_3()
        lmm.set_random_effect([r['opt_k']])
    res = lmm.get_REML()
    res['Y_t'] = res['H_sqrt_inv'] * lmm.Y  # The transformed outputs.
    res['X_t'] = res['H_sqrt_inv'] * lmm.X
    res['lmm'] = lmm
    if K2 != None:
        res['perc_var1'] = r['perc_var1']
        res['perc_var2'] = r['perc_var2']
    return res


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
    pval = stats.chi2.sf(D, 1)
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

    print "Running EMMA (python)"
    s1 = time.time()
    res = lmm.expedited_REML_t_test(snps)
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
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
    print "Running EMMAX"
    s1 = time.time()
    res = lmm.emmax_f_test(snps, Z=Z, with_betas=with_betas, progress_file_writer=progress_file_writer,
            emma_num=emma_num)
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    return res



def emmax_perm_test(snps, phenotypes, K, num_perm=100):
    """
    Run EMMAX
    """
    lmm = LinearMixedModel(phenotypes)
    lmm.add_random_effect(K)
    print "Running %d EMMAX-permutation (writes %d dots)" % (num_perm, num_perm)
    s1 = time.time()
    res = lmm.emmax_permutations(snps, num_perm)
    p_f_list = zip(res['min_ps'], res['max_f_stats'])
    p_f_list.sort()
#    print p_f_list[:10]
    threshold = p_f_list[len(p_f_list) / 20]
    res['threshold_05'] = threshold
    print 'Tresholds should be:', threshold
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    return res


def anova(snps, phenotypes):
    """
    Run EMMAX
    """
    lmm = LinearModel(phenotypes)

    print "Running ANOVA"
    s1 = time.time()
    res = lmm.anova_f_test(snps)
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    return res



def emmax_anova(snps, phenotypes, K):
    """
    Run EMMAX
    """
    lmm = LinearMixedModel(phenotypes)
    lmm.add_random_effect(K)

    print "Running EMMAX-ANOVA"
    s1 = time.time()
    res = lmm.emmax_anova_f_test(snps)
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    return res


#
# def _get_interactions_(isnp, snps, mac_threshold=15):
#    isnps = []
#    cofactor_indices = []
#    anti_isnp = sp.vectorize(lambda x: 1 if x == 0 else 0)(isnp)
#    for i, snp in enumerate(snps):
#        min_count = min(min(sp.bincount(isnp & snp)), min(sp.bincount(isnp | snp)), \
#            min(sp.bincount(snp & anti_isnp)), min(sp.bincount(snp | anti_isnp)))
#        if min_count > mac_threshold and min:
#            isnps.append(isnp & snp)
#            cofactor_indices.append(i)
#    return isnps, cofactor_indices



def _log_choose_(n, k):
     # log_fact = lambda n: n * sp.log(n) - n + (sp.log(n * (1 + 4 * n * (1 + 2 * n))) / 6) + sp.log(sp.pi) / 2
     # Srinivasa Ramanujan approximation of log(n!)
    if k == 0 or n == k:
        return 0
    if n < k:
        raise Exception('Out of range.')
    return sp.sum(map(sp.log, range(n, n - k, -1))) - sp.sum(map(sp.log, range(k, 0, -1)))


def _calc_bic_(ll, num_snps, num_par, n):
    bic = -2 * (ll) + num_par * sp.log(n)
    extended_bic = bic + \
        2 * _log_choose_(num_snps, num_par - 2)
    modified_bic = bic + \
        2 * (num_par) * sp.log(num_snps / 2.2 - 1)
    return (bic, extended_bic, modified_bic)


def _plot_manhattan_and_qq_(file_prefix, step_i, pvals, quantiles_dict=None, plot_bonferroni=True,
            highlight_markers=None, cand_genes=None, plot_xaxis=True, log_qq_max_val=5, with_qq_plots=True,
            num_dots=1000, simple_qq=False, highlight_loci=None, write_pvals=False,
            highlight_ppa_markers=None, ppas=None, markersize=3, chrom_col_map=None, perc_var_expl=None,
            save_pvals=False, **kwargs):
    import pylab
    cm = pylab.get_cmap('hsv')

    png_file_name = '%s_step%d.png' % (file_prefix, step_i)
    res = gr.Result(scores=pvals, perc_var_expl=perc_var_expl, **kwargs)


    res.plot_manhattan(png_file=png_file_name, plot_bonferroni=True, highlight_markers=highlight_markers,
                cand_genes=cand_genes, plot_xaxis=plot_xaxis, highlight_loci=highlight_loci,
                markersize=markersize, chrom_col_map=chrom_col_map, neg_log_transform=True)
    ret_dict = {'manhattan':png_file_name}


    if ppas != None:
        png_file_name = '%s_ppa_step%d.png' % (file_prefix, step_i)

        ppas_res = gr.Result(scores=ppas, perc_var_expl=perc_var_expl, **kwargs)
        ppas_res.plot_manhattan(png_file=png_file_name, highlight_markers=highlight_ppa_markers,
                    cand_genes=cand_genes, plot_xaxis=plot_xaxis,
                    highlight_loci=highlight_loci, ylab='Post. prob. of assoc.',
                    min_score=0.0, max_score=1.0, markersize=3)
        ret_dict['ppa_manhattan'] = png_file_name

    # Plot QQ plot
    if with_qq_plots and quantiles_dict != None:
        # calculate qq-plot line..
        log_quantiles = agr.get_log_quantiles(pvals, num_dots=num_dots, max_val=log_qq_max_val)
        quantiles = agr.get_quantiles(pvals, num_dots=num_dots)
        quantiles_dict['log'].append(log_quantiles)
        quantiles_dict['norm'].append(quantiles)
        quantiles_dict['labels'].append('Step %d' % (step_i))
        # plot all lines
        log_png_file_name = '%s_step%d_log_qqplot.png' % (file_prefix, step_i)
        png_file_name = '%s_step%d_qqplot.png' % (file_prefix, step_i)
        if step_i == 0:
            color_list = [cm(0.0)]
        else:
            color_list = [cm(i / float(step_i) * 0.7) for i in range(step_i + 1)]
        if simple_qq:
            log_qs = [quantiles_dict['log'][0], log_quantiles]
            qs = [quantiles_dict['norm'][0], quantiles]
            q_labs = [quantiles_dict['labels'][0], quantiles_dict['labels'][-1]]
            lcols = [cm(0), cm(0.7)]
            agr.simple_log_qqplot(log_qs, log_png_file_name, quantile_labels=q_labs, line_colors=lcols,
                    num_dots=num_dots, max_val=log_qq_max_val)
            agr.simple_qqplot(qs, png_file_name, quantile_labels=q_labs, line_colors=color_list, num_dots=num_dots)

        elif step_i > 6:
            log_qs = quantiles_dict['log'][0:6] + [log_quantiles]
            qs = quantiles_dict['norm'][0:6] + [quantiles]
            q_labs = quantiles_dict['labels'][0:6] + [quantiles_dict['labels'][-1]]
            lcols = [cm(i / 6.0 * 0.7) for i in range(7)]
            agr.simple_log_qqplot(log_qs, log_png_file_name, quantile_labels=q_labs, line_colors=lcols, num_dots=num_dots,
                    max_val=log_qq_max_val)
            agr.simple_qqplot(qs, png_file_name, quantile_labels=q_labs, line_colors=lcols, num_dots=num_dots)
        else:
            agr.simple_log_qqplot(quantiles_dict['log'], log_png_file_name, quantile_labels=quantiles_dict['labels'],
                    line_colors=color_list, num_dots=num_dots, max_val=log_qq_max_val)
            agr.simple_qqplot(quantiles_dict['norm'], png_file_name, quantile_labels=quantiles_dict['labels'],
                    line_colors=color_list, num_dots=num_dots)
        ret_dict['log_qq'] = log_png_file_name
        ret_dict['qq'] = png_file_name
    if write_pvals:
        pval_file_name = '%s_step%d.pvals' % (file_prefix, step_i)
        res.write_to_file(pval_file_name, additional_columns='perc_var_expl')
        if ppas != None:
        #    ppas_res.filter_percentile(0.9)
            pval_file_name = '%s_step%d.ppas' % (file_prefix, step_i)
            ppas_res.write_to_file(pval_file_name, additional_columns='perc_var_expl')
        if save_pvals:
            res.filter_percentile(0.02, reversed=True)
            ret_dict['res'] = res
    return ret_dict


def _analyze_opt_criterias_(criterias, sign_threshold, max_num_cofactors, file_prefix, with_qq_plots,
                            lm, step_info_list, quantiles_dict, plot_bonferroni=True, cand_genes=None,
                            plot_xaxis=True, log_qq_max_val=5, eig_L=None, type='emmax', highlight_loci=None,
                            write_pvals=False, snp_priors=None, ppa_threshold=0.5, emma_num=None,
                            save_pvals=False , **kwargs):
    """
    Copies or plots optimal criterias
    """
    ret_dict = {}
    opt_indices = {}
    opt_file_dict = {}
    for c in criterias:
        print 'GWAs for optimal %s criteria:' % c
        if c == 'bonf':
            opt_list = sp.arange(max_num_cofactors + 1)
            for i, pval in enumerate(criterias['bonf']):
                if pval > sign_threshold:
                    opt_list[i] = -1
            i_opt = opt_list.argmax()
        elif c == 'mbonf':
            fw_opt_list = sp.arange(max_num_cofactors + 1)
            for i in range(max_num_cofactors + 1):
                pval = criterias[c][i]
                if pval > sign_threshold:
                    fw_opt_list[i] = -1
            fw_i_opt = fw_opt_list.argmax()
            fw_max = fw_opt_list[fw_i_opt]

            if max_num_cofactors > 1:
                shift = max_num_cofactors + 1
                bw_opt_list = sp.arange(max_num_cofactors - 1, 0, -1)
                for i in range(len(bw_opt_list)):
                    pval = criterias[c][i + shift]
                    if pval > sign_threshold:
                        bw_opt_list[i] = -1
                bw_i_opt = bw_opt_list.argmax()
                bw_max = bw_opt_list[bw_i_opt]
                bw_i_opt = bw_opt_list.argmax() + shift
                if bw_max == fw_max:
                    i_opt = bw_i_opt if criterias[c][fw_i_opt] > criterias[c][bw_i_opt] else fw_i_opt
                else:
                    i_opt = bw_i_opt if bw_max > fw_max else fw_i_opt
            else:
                i_opt = fw_i_opt
        elif c == 'min_cof_ppa':
            fw_opt_list = sp.arange(max_num_cofactors + 1)
            for i in range(max_num_cofactors + 1):
                ppa = criterias[c][i]
                if ppa < ppa_threshold:
                    fw_opt_list[i] = -1
            fw_i_opt = fw_opt_list.argmax()
            fw_max = fw_opt_list[fw_i_opt]

            if max_num_cofactors > 1:
                shift = max_num_cofactors + 1
                bw_opt_list = sp.arange(max_num_cofactors - 1, 0, -1)
                for i in range(len(bw_opt_list)):
                    ppa = criterias[c][i + shift]
                    if ppa < ppa_threshold:
                        bw_opt_list[i] = -1
                bw_i_opt = bw_opt_list.argmax()
                bw_max = bw_opt_list[bw_i_opt]
                bw_i_opt = bw_opt_list.argmax() + shift
                if bw_max == fw_max:
                    i_opt = bw_i_opt if criterias[c][fw_i_opt] > criterias[c][bw_i_opt] else fw_i_opt
                else:
                    i_opt = bw_i_opt if bw_max > fw_max else fw_i_opt
            else:
                i_opt = fw_i_opt

        else:
            cur_min_val = criterias[c][0]
            min_indices = [0]
            for i in range(1, len(criterias[c])):
                v = criterias[c][i]
                if v < cur_min_val:
                    cur_min_val = v
                    min_indices = [i]
                if v == cur_min_val:
                    min_indices.append(i)
            i_opt = min(min_indices)
            # i_opt = sp.argmin(criterias[c])
        print "    %d'th step was optimal." % i_opt
        ret_dict[c] = i_opt
        if i_opt <= max_num_cofactors:
            # Copy the pngs...
            if file_prefix:
                png_file_name = '%s_step%d.png' % (file_prefix, i_opt)
                opt_png_file_name = '%s_step%d_opt_%s.png' % (file_prefix, i_opt, c)
                if platform.system() == 'Linux' or platform.system() == 'Darwin':
                    # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)
                    if snp_priors != None:
                        png_file_name = '%s_ppa_step%d.png' % (file_prefix, i_opt)
                        opt_png_file_name = '%s_ppa_step%d_opt_%s.png' % (file_prefix, i_opt, c)
                        # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)
                    if with_qq_plots:
                        qq_png_file_name = '%s_step%d_qqplot.png' % (file_prefix, i_opt)
                        opt_qq_png_file_name = '%s_step%d_opt_%s_qqplot.png' % (file_prefix, i_opt, c)
                        # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', qq_png_file_name, opt_qq_png_file_name)
                        log_qq_png_file_name = '%s_step%d_log_qqplot.png' % (file_prefix, i_opt)
                        opt_log_qq_png_file_name = '%s_step%d_opt_%s_log_qqplot.png' % (file_prefix, i_opt, c)
                        # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', log_qq_png_file_name, opt_log_qq_png_file_name)
        elif i_opt in opt_file_dict:
            if file_prefix:
                png_file_name = opt_file_dict[i_opt]['manhattan']
                opt_png_file_name = '%s_step%d_opt_%s.png' % (file_prefix, i_opt, c)
                if platform.system() == 'Linux' or platform.system() == 'Darwin':
                    # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)
                    if snp_priors != None:
                        png_file_name = opt_file_dict[i_opt]['ppa_manhattan']
                        opt_png_file_name = '%s_ppa_step%d_opt_%s.png' % (file_prefix, i_opt, c)
                        # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', png_file_name, opt_png_file_name)

                    if with_qq_plots:
                        qq_png_file_name = opt_file_dict[i_opt]['qq']
                        opt_qq_png_file_name = '%s_step%d_opt_%s_qqplot.png' % (file_prefix, i_opt, c)
                        # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', qq_png_file_name, opt_qq_png_file_name)
                        log_qq_png_file_name = opt_file_dict[i_opt]['log_qq']
                        opt_log_qq_png_file_name = '%s_step%d_opt_%s_log_qqplot.png' % (file_prefix, i_opt, c)
                        # os.spawnlp(os.P_NOWAIT, 'cp', 'cp', log_qq_png_file_name, opt_log_qq_png_file_name)

        elif not i_opt in opt_indices:
            # Perfom GWAS witht the optimal cofactors
            cofactor_snps = step_info_list[i_opt]['cofactor_snps']
            cofactors = step_info_list[i_opt]['cofactors']
            print cofactors
            lm.set_factors(cofactor_snps)
            if type == 'emmax':
                eig_R = lm._get_eigen_R_(X=lm.X)
                reml_res = lm.get_REML(eig_L=eig_L, eig_R=eig_R)
                H_sqrt_inv = reml_res['H_sqrt_inv']
                l_res = lm._emmax_f_test_(kwargs['snps'], H_sqrt_inv, snp_priors=snp_priors,
                                emma_num=emma_num)
                min_pval_i = l_res['ps'].argmin()
                mahalnobis_rss = l_res['rss'][min_pval_i]
                print 'Min Mahalanobis RSS:', mahalnobis_rss
            elif type == 'lm':
                l_res = lm.fast_f_test(kwargs['snps'])
                min_pval_i = l_res['ps'].argmin()

            min_pval = l_res['ps'][min_pval_i]

            min_pval_chr_pos = (kwargs['chromosomes'][min_pval_i], kwargs['positions'][min_pval_i])
            print 'Min p-value:', min_pval
            l_pvals = l_res['ps'].tolist()
            l_perc_var_expl = l_res['var_perc'].tolist()
            opt_indices[i_opt] = {'min_pval':min_pval, 'min_pval_chr_pos':min_pval_chr_pos,
                        'kolmogorov_smirnov':agr.calc_ks_stats(l_pvals),
                        'pval_median':agr.calc_median(l_pvals)}
            if file_prefix:
                opt_file_prefix = '%s_opt_%s' % (file_prefix, c)
                if snp_priors:
                    ppa_cofactors = step_info_list[i_opt]['ppa_cofactors']
                    ppas = l_res['ppas'].tolist()
                else:
                    ppas = None
                    ppa_cofactors = None
                opt_file_dict[i_opt], res = _plot_manhattan_and_qq_(opt_file_prefix, i_opt, l_pvals, quantiles_dict,
                                plot_bonferroni=True, highlight_markers=cofactors,
                                cand_genes=cand_genes, plot_xaxis=plot_xaxis,
                                log_qq_max_val=log_qq_max_val, with_qq_plots=with_qq_plots,
                                simple_qq=True, highlight_loci=highlight_loci,
                                write_pvals=write_pvals, highlight_ppa_markers=ppa_cofactors,
                                ppas=ppas, perc_var_expl=l_perc_var_expl, save_pvals=save_pvals,
                                **kwargs)
                if save_pvals:
                    opt_indices['res'] = opt_file_dict[i_opt]['res']

            if type == 'emmax':
                opt_indices[i_opt]['mahalanobis_rss'] = mahalnobis_rss
            opt_indices[i_opt]['res'] = res
    return ret_dict, opt_indices


def _cofactors_to_string_(cofactors):
    st = ''
    if len(cofactors) > 0:
        for tup in cofactors[:-1]:
            st += '%d_%d_%f,' % (tup[0], tup[1], tup[2])
        st += '%d_%d_%f' % (cofactors[-1][0], cofactors[-1][1], cofactors[-1][2])
    return st


def _plot_stepwise_stats_(file_prefix, step_info_list, sign_threshold, type='emmax'):
    import pylab
    rss_list = []
    ll_list = []
    bic_list = []
    e_bic_list = []
    m_bic_list = []
    mbonf_list = []
    min_pval_list = []
    ks_stat_list = []
    pval_median_list = []
    full_steps = []
    d_keys = [ 'rss', 'll', 'bic', 'e_bic', 'm_bic', 'min_pval', 'mbonf', 'pval_median']
    if type == 'emmax':
        p_her_list = []
        reml_mahalanobis_rss_list = []
        mahalanobis_rss_list = []
        d_keys += ['pseudo_heritability', 'reml_mahalanobis_rss', 'mahalanobis_rss']

    f = open(file_prefix + "_stats.csv", 'w')
    f.write(','.join(['step_nr'] + d_keys + ['kolmogorov_smirnov', 'min_pval_pos_chr', 'cofactors']) + '\n')
    for i, si in enumerate(step_info_list):
        st = ','.join(map(str, [i] + [si[k] for k in d_keys]))
        if si['kolmogorov_smirnov']:
            st += ',%d' % si['kolmogorov_smirnov']['D']
        else:
            st += ','
        if si['min_pval_chr_pos']:
            st += ',%d_%d,' % si['min_pval_chr_pos']
        else:
            st += ',,'
        st += _cofactors_to_string_(si['cofactors'])
        st += '\n'
        f.write(st)
        rss_list.append(float(si['rss']))
        if type == 'emmax':
            p_her_list.append(float(si['pseudo_heritability']))
            reml_mahalanobis_rss_list.append(float(si['reml_mahalanobis_rss']))
            if si['mahalanobis_rss']:
                mahalanobis_rss_list.append(float(si['mahalanobis_rss']))
        ll_list.append(si['ll'])
        bic_list.append(si['bic'])
        e_bic_list.append(si['e_bic'])
        m_bic_list.append(si['m_bic'])
        mbonf_list.append(si['mbonf'])
        if si['min_pval']:
            min_pval_list.append(float(si['min_pval']))
            ks_stat_list.append(si['kolmogorov_smirnov']['D'])
            pval_median_list.append(si['pval_median'])
            full_steps.append(i)
    f.close()
    pylab.figure(figsize=(6, 4))
    pylab.axes([0.15, 0.1, 0.83, 0.85])
    num_steps = len(step_info_list)
    rss_list.append(rss_list[0])
    y_max = max(rss_list)
    pylab.axis([-0.025 * num_steps, 1.025 * num_steps, 0, 1.025 * y_max])
    pylab.plot(range(len(rss_list)), rss_list, 'o-', alpha=0.7)
    pylab.ylabel('RSS')
    pylab.xlabel('Number of steps')
    pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
    pylab.savefig(file_prefix + '_stats_rss.png')
    pylab.clf()
    pylab.axes([0.15, 0.1, 0.83, 0.85])
    pylab.axis([-0.025 * num_steps / 2, 1.025 * num_steps / 2, 0, 1.025 * y_max])
    pylab.plot(range((num_steps / 2) + 1), rss_list[:(num_steps / 2) + 1], 'o-', alpha=0.7)
    pylab.ylabel('RSS')
    pylab.xlabel('Number of steps')
    pylab.savefig(file_prefix + '_stats_rss2.png')
    pylab.clf()
    ll_list.append(ll_list[0])
    pylab.plot(range(len(ll_list)), ll_list, 'o-', alpha=0.7)
    pylab.ylabel('Log likelihood')
    pylab.xlabel('Number of steps')
    pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
    pylab.savefig(file_prefix + '_stats_ll.png')
    pylab.clf()
    bic_list.append(bic_list[0])
    pylab.plot(range(len(bic_list)), bic_list, 'o-', alpha=0.7)
    pylab.ylabel('BIC')
    pylab.xlabel('Number of steps')
    pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
    pylab.savefig(file_prefix + '_stats_bic.png')
    pylab.clf()
    e_bic_list.append(e_bic_list[0])
    pylab.plot(range(len(e_bic_list)), e_bic_list, 'o-', alpha=0.7)
    pylab.ylabel('Extended BIC')
    pylab.xlabel('Number of steps')
    pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
    pylab.savefig(file_prefix + '_stats_ebic.png')
    pylab.clf()
    m_bic_list.append(m_bic_list[0])
    pylab.plot(range(len(m_bic_list)), m_bic_list, 'o-', alpha=0.7)
    pylab.ylabel('Modified BIC')
    pylab.xlabel('Number of steps')
    pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
    pylab.savefig(file_prefix + '_stats_mbic.png')
    pylab.clf()
    pylab.plot(full_steps, -sp.log10(min_pval_list), 'o-', alpha=0.7)
    pylab.axhline(y= -sp.log10(sign_threshold), c='k', linestyle='--')
    pylab.ylabel('Min. p-value')
    pylab.xlabel('Number of steps')
    pylab.savefig(file_prefix + '_stats_pval.png')
    pylab.clf()
    pylab.plot(full_steps, ks_stat_list, 'o-', alpha=0.7)
    pylab.ylabel('Kolmogorov-Smirnov statistic')
    pylab.xlabel('Number of steps')
    pylab.savefig(file_prefix + '_stats_ks.png')
    pylab.clf()
    pylab.plot(full_steps, pval_median_list, 'o-', alpha=0.7)
    pylab.axhline(y=0, c='k', linestyle='--')
    pylab.ylabel('0.5 - median p-value')
    pylab.xlabel('Number of steps')
    pylab.savefig(file_prefix + '_stats_mpval.png')
    pylab.clf()
    pylab.plot(range(1, len(mbonf_list)), -sp.log10(mbonf_list[1:]), 'o-', alpha=0.7)
    pylab.axhline(y= -sp.log10(sign_threshold), c='k', linestyle='--')
    pylab.ylabel('Max cofactor pvalue (-log10[pval])')
    pylab.xlabel('Number of steps')
    pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
    pylab.savefig(file_prefix + '_stats_mbonf.png')
    pylab.clf()

    # Plotting variance partition plots
    if type == 'emmax':
        p_her_list.append(p_her_list[0])
        pylab.plot(range(len(p_her_list)), p_her_list, 'o-', alpha=0.7)
        pylab.ylabel('Pseudo-heritability')
        pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
        pylab.savefig(file_prefix + '_stats_p_her.png')
        pylab.clf()
        reml_mahalanobis_rss_list.append(reml_mahalanobis_rss_list[0])
        pylab.plot(range(len(reml_mahalanobis_rss_list)), reml_mahalanobis_rss_list, 'o-', alpha=0.7)
        pylab.ylabel('REML Mahalanobis RSS')
        pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
        pylab.savefig(file_prefix + '_stats_reml_mahalanobis_rss.png')
        pylab.clf()
        pylab.plot(full_steps, mahalanobis_rss_list, 'o-', alpha=0.7)
        pylab.ylabel('Mahalanobis RSS')
        pylab.savefig(file_prefix + '_stats_mahalanobis_rss.png')
        pylab.clf()

        max_rss = max(rss_list)
        rss_array = sp.array(rss_list) / max_rss
        p_her_array = rss_array * sp.array(p_her_list)
        genetic_variance = p_her_array + (1 - rss_array)
        variance_explained = (1 - rss_array)
        pylab.figure(figsize=(10, 2.5))
        pylab.axes([0.07, 0.18, 0.92, 0.78])
        pylab.fill_between([0, num_steps], 0, 1.1, color='#DD3333', alpha=0.8, label='Error')
        pylab.fill_between(sp.arange(num_steps + 1), 0, genetic_variance, color='#22CC44', alpha=0.8, label='Genetic variance')
        pylab.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained, color='#2255AA', alpha=0.8, label='Variance explained')
        pylab.ylabel('Partition of variance')
        pylab.xlabel('Step number')
        pylab.axvline(x=num_steps / 2, c='k', linestyle=':')
        pylab.legend(loc=1, ncol=3, shadow=True)
        pylab.axis([0, num_steps, 0, 1])
        pylab.savefig(file_prefix + '_stats_variances.png', format='png')

        pylab.figure(figsize=(6, 2.5))
        pylab.axes([0.12, 0.18, 0.86, 0.78])
        num_steps = num_steps / 2
        pylab.fill_between([0, num_steps ], 0, 1.1, color='#DD3333', alpha=0.8, label='Error')
        pylab.fill_between(sp.arange(num_steps + 1), 0, genetic_variance[:num_steps + 1], color='#22CC44', alpha=0.8, label='Genetic variance')
        pylab.fill_between(sp.arange(num_steps + 1), -0.1, variance_explained[:num_steps + 1], color='#2255AA', alpha=0.8, label='Variance explained')
        pylab.ylabel('Partition of variance')
        pylab.xlabel('Step number')
        pylab.legend(loc=1, ncol=3, shadow=True)
        pylab.axis([0, num_steps, 0, 1])
        pylab.savefig(file_prefix + '_stats_variances_forward.png', format='png')


def lin_reg_step(phen_vals, sd, cof_chr_pos_list, progress_file_writer=None, plot_prefix=None):
    """
    Standard linear regression single SNPs..
    
    Returns various stats useful for stepwise regression.
    """
    import gwaResults as gr
    import bisect
    s1 = time.time()
    lm = LinearModel(phen_vals)

    h0_rss = lm.get_rss()
    step_dict = {}
    step_dict['h0_rss'] = h0_rss

    chrom_pos_list = sd.get_chr_pos_list()
    print 'Looking up cofactors'
    cof_indices = []
    for chrom_pos in cof_chr_pos_list:
        i = bisect.bisect(chrom_pos_list, chrom_pos) - 1
        assert chrom_pos_list[i] == chrom_pos, 'Cofactor missing??'
        cof_indices.append(i)
    snps = sd.get_snps()
    cof_snps = [snps[i] for i in cof_indices]
    lm.set_factors(cof_snps)

    if not progress_file_writer == None:
        progress_file_writer.update_progress_bar(progress=0.20, task_status='Performing linear regression')
        progress_file_writer.set_step(0.05)


    print 'Performing linear regression.'
    r = lm.fast_f_test(snps, progress_file_writer=progress_file_writer)
    min_pval_i = sp.argmin(r['ps'])
    step_dict['min_pval_i'] = min_pval_i
    step_dict['min_pval'] = r['ps'][min_pval_i]
    step_dict['mahalanobis_rss'] = r['rss'][min_pval_i]
    step_dict['min_pval_chr_pos'] = chrom_pos_list[min_pval_i]

    num_snps = len(snps)
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
    print 'Updating the cofactor p-values'
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
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)

    if plot_prefix != None:
        res = gr.Result(scores=list(r['ps']), snps_data=sd)
        png_file_name = '%s_%s_manhattan.png' % (plot_prefix, _cofactors_to_string_(cof_chrom_pos_pval_list))
        res.neg_log_trans()
        res.filter_percentile(0.1)
        res.plot_manhattan(png_file=png_file_name, plot_bonferroni=True, highlight_markers=cof_chrom_pos_pval_list)

    return {'stats':step_dict, 'res':r}





def emmax_step(phen_vals, sd, K, cof_chr_pos_list, eig_L=None, eig_R=None, progress_file_writer=None, plot_prefix=None,
        emma_num=100):
    """
    EMMAX single SNPs
    
    Returns various stats useful for stepwise regression.
    """

    import gwaResults as gr
    import bisect
    s1 = time.time()
    lmm = LinearMixedModel(phen_vals)
    lmm.add_random_effect(K)

    if not eig_L:
        print 'Calculating the eigenvalues of K'
        eig_L = lmm._get_eigen_L_()
        print 'Done.'
    reml_dict = lmm.get_REML(eig_L=eig_L)
    h0_rss = float(reml_dict['rss'])

    chrom_pos_list = sd.get_chr_pos_list()
    print 'Looking up cofactors'
    cof_indices = []
    for chrom_pos in cof_chr_pos_list:
        i = bisect.bisect(chrom_pos_list, chrom_pos) - 1
        assert chrom_pos_list[i] == chrom_pos, 'Cofactor missing??'
        cof_indices.append(i)
    snps = sd.get_snps()
    cof_snps = [snps[i] for i in cof_indices]
    lmm.set_factors(cof_snps)
    if not eig_R:
        print "Calculating the eigenvalues of S(K+I)S where S = I-X(X'X)^-1X'"
        eig_R = lmm._get_eigen_R_()
        print 'Done'


    print 'Getting variance components estimates'
    reml_dict = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
    ml_dict = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
    print 'Done.'
    print 'pseudo_heritability:', reml_dict['pseudo_heritability']
    H_sqrt_inv = reml_dict['H_sqrt_inv']

    if not progress_file_writer == None:
        progress_file_writer.update_progress_bar(progress=0.20, task_status='Performing AMM')
        progress_file_writer.set_step(0.05)


    r = lmm._emmax_f_test_(snps, H_sqrt_inv, progress_file_writer=progress_file_writer, emma_num=emma_num)
    min_pval_i = sp.argmin(r['ps'])
    step_dict = {}
    step_dict['min_pval_i'] = min_pval_i
    step_dict['min_pval'] = r['ps'][min_pval_i]
    step_dict['mahalanobis_rss'] = r['rss'][min_pval_i]
    step_dict['min_pval_chr_pos'] = chrom_pos_list[min_pval_i]
    step_dict['h0_rss'] = h0_rss

    # Update the variance components
#    lmm.add_factor(snps[min_pval_i])
#    eig_R = lmm._get_eigen_R_()
#    print 'Updating the variance components estimates'
#    reml_dict = lmm.get_REML(eig_L=eig_L, eig_R=eig_R)
#    ml_dict = lmm.get_ML(eig_L=eig_L, eig_R=eig_R)
#    print 'Done.'
    ll = ml_dict['max_ll']
    step_dict['max_ll'] = ll
    num_snps = len(snps)
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
#    rss = lmm.get_rss()[0]

    step_dict['rss'] = float(reml_dict['rss'])
    step_dict['reml_mahalanobis_rss'] = reml_dict['mahalanobis_rss']
    perc_var_expl = 1.0 - (step_dict['rss'] / h0_rss)
    step_dict['perc_var_expl'] = perc_var_expl
    step_dict['remain_perc_gen_var'] = (1 - perc_var_expl) * p_her
    step_dict['remain_perc_err_var'] = (1 - perc_var_expl) * (1 - p_her)

    # Calculate maximum cofactor p-value\
    print 'Updating the cofactor p-values'
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
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)

    if plot_prefix != None:
        res = gr.Result(scores=list(r['ps']), snps_data=sd)
        png_file_name = '%s_%s_manhattan.png' % (plot_prefix, _cofactors_to_string_(cof_chrom_pos_pval_list))
        res.neg_log_trans()
        res.filter_percentile(0.1)
        res.plot_manhattan(png_file=png_file_name, plot_bonferroni=True, highlight_markers=cof_chrom_pos_pval_list)

    return {'stats':step_dict, 'res':r, 'upd_H_sqrt_inv':H_sqrt_inv}




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
    print "Running a standard linear model"
    s1 = time.time()
    res = lm.fast_f_test(snps)
    secs = time.time() - s1
    if secs > 60:
        mins = int(secs) / 60
        secs = secs - mins * 60
        print 'Took %d mins and %f seconds.' % (mins, secs)
    else:
        print 'Took %f seconds.' % (secs)
    return res




def local_vs_global_mm(y, local_k, global_k, K, h0_res=None):
    """
    Local vs. global kinship mixed model.
    """
    if h0_res == None:
        lmm0 = LinearMixedModel(Y=y)
        lmm0.add_random_effect(K)
        eig_L = lmm0._get_eigen_L_()
        h0_res = lmm0.get_estimates(eig_L)

    lmm = LinearMixedModel(Y=y)
    lmm.add_random_effect(global_k)
    lmm.add_random_effect(local_k)
    h1_res = lmm.get_estimates_3()
    lrt_stat = 2 * (h1_res['max_ll'] - h0_res['max_ll'])
    pval = stats.chi2.sf(lrt_stat, 1)
    # print 'p-value:', pval
    res_dict = {'pval':pval, 'perc_var1':h1_res['perc_var1'], 'perc_var2':h1_res['perc_var2'],
            'pseudo_heritability0':h0_res['pseudo_heritability'],
            'pseudo_heritability1':h1_res['pseudo_heritability']}
    return res_dict


def chrom_vs_rest_mm(y, sd, kinship_method='ibd', global_k=None):
    if global_k == None:
        if kinship_method == 'ibd':
            K = sd.get_ibd_kinship_matrix()
        elif kinship_method == 'ibs':
            K = sd.get_ibs_kinship_matrix()
        else:
            raise NotImplementedError
    else:
        K = global_k
    lmm0 = LinearMixedModel(Y=y)
    lmm0.add_random_effect(K)
    eig_L = lmm0._get_eigen_L_()
    h0_res = lmm0.get_estimates(eig_L)

    chromosomes = []
    pvals = []
    perc_variances1 = []
    perc_variances2 = []
    h1_heritabilities = []
    for chrom in sd.chromosomes:
        d = sd.get_chrom_vs_rest_kinships(chrom, global_kinship=K, kinship_method=kinship_method)
        if d['local_k'] != None and d['global_k'] != None:
            local_k = kinship.scale_k(d['local_k'])
            global_k = kinship.scale_k(d['global_k'])
            res_dict = local_vs_global_mm(y, local_k, global_k, K, h0_res=h0_res)
            chromosomes.append(chrom)
            perc_variances1.append(res_dict['perc_var1'])
            perc_variances2.append(res_dict['perc_var2'])
            h1_heritabilities.append(res_dict['pseudo_heritability1'])
            pvals.append(res_dict['pval'])
    return {'pvals':pvals, 'perc_variances2':perc_variances2, 'perc_variances1':perc_variances1,
        'h0_heritability':h0_res['pseudo_heritability'], 'h1_heritabilities':h1_heritabilities,
        'chromosomes':chromosomes}


def local_vs_global_mm_scan(y, sd, file_prefix='/tmp/temp', window_size=1000000, jump_size=500000, kinship_method='ibd', global_k=None):
    """
    Local vs. global kinship mixed model.
    """
    print 'Starting Mixed model, local vs. global kinship scan...'
    print 'window size is %d, and jump size is %d' % (window_size, jump_size)
    import gwaResults as gr
    if global_k == None:
        if kinship_method == 'ibd':
            K = sd.get_ibd_kinship_matrix()
        elif kinship_method == 'ibs':
            K = sd.get_ibs_kinship_matrix()
        else:
            raise NotImplementedError
    else:
        K = global_k
    lmm0 = LinearMixedModel(Y=y)
    lmm0.add_random_effect(K)
    eig_L = lmm0._get_eigen_L_()
    h0_res = lmm0.get_estimates(eig_L)

    genome_length = sd.get_genome_length()
    est_num_chunks = genome_length / jump_size
    chromosomes = []
    positions = []
    pvals = []
    perc_variances1 = []
    perc_variances2 = []
    h1_heritabilities = []
    chunk_i = 0
    for sdl, chrom in zip(sd.snpsDataList, sd.chromosomes):
        for focal_pos in range(sdl.positions[0], sdl.positions[-1], jump_size):
            chunk_i += 1
            d = sd.get_local_n_global_kinships((chrom, focal_pos), window_size,
                                        global_kinship=K,
                                        kinship_method=kinship_method)
            if d['local_k'] != None and d['global_k'] != None:
                local_k = kinship.scale_k(d['local_k'])
                global_k = kinship.scale_k(d['global_k'])
                # print "Chromosome=%d, position=%d" % (chrom, focal_pos)
                res_dict = local_vs_global_mm(y, local_k, global_k, K, h0_res=h0_res)
                chromosomes.append(chrom)
                positions.append(focal_pos)
                perc_variances1.append(res_dict['perc_var1'])
                perc_variances2.append(res_dict['perc_var2'])
                h1_heritabilities.append(res_dict['pseudo_heritability1'])
                pvals.append(res_dict['pval'])

                # print 'H0: pseudo_heritability=%0.2f' % (res_dict['h0_res']['pseudo_heritability'])
                # print 'H1: pseudo_heritability=%0.2f, perc_var1=%0.2f, perc_var2=%0.2f' % \
                #        (res_dict['h1_res']['pseudo_heritability'],
                #        res_dict['h1_res']['perc_var1'],
                #        res_dict['h1_res']['perc_var2'])
            if est_num_chunks >= 100 and (chunk_i + 1) % int(est_num_chunks / 100) == 0:  # Print dots
                sys.stdout.write('.')
                sys.stdout.flush()
            elif est_num_chunks < 100:
                print chunk_i

    pval_res = gr.Result(scores=pvals, positions=positions, chromosomes=chromosomes)
    pval_res.neg_log_trans()
    pval_res.plot_manhattan(png_file=file_prefix + '_lrt_pvals.png', percentile=0, plot_bonferroni=True)
    perc_var_res = gr.Result(scores=perc_variances2, positions=positions, chromosomes=chromosomes)
    perc_var_res.plot_manhattan(png_file=file_prefix + '_perc_var_explained.png', percentile=0,
                ylab='% of variance explained')
    return {'pvals':pvals, 'perc_variances2':perc_variances2, 'perc_variances1':perc_variances1,
        'h0_heritability':h0_res['pseudo_heritability'], 'h1_heritabilities':h1_heritabilities,
        'chromosomes':chromosomes, 'positions':positions}



def local_vs_global_gene_mm_scan(y, sd, file_prefix='/tmp/temp', radius=20000, kinship_method='ibd',
                global_k=None, tair_ids=None, plot_gene_trees=False, ets=None):
    """
    Local vs. global kinship mixed model.
    """
    print 'Starting Mixed model, local vs. global kinship scan...'
    import gwaResults as gr
    import dataParsers as dp
    if global_k == None:
        if kinship_method == 'ibd':
            K = sd.get_ibd_kinship_matrix()
        elif kinship_method == 'ibs':
            K = sd.get_ibs_kinship_matrix()
        else:
            raise NotImplementedError
    else:
        K = global_k
    lmm0 = LinearMixedModel(Y=y)
    lmm0.add_random_effect(K)
    eig_L = lmm0._get_eigen_L_()
    h0_res = lmm0.get_estimates(eig_L)

    gene_dict = dp.parse_tair_gff_file()
    if tair_ids == None:
        tair_ids = gene_dict.keys()
    tair_ids.sort()
    chromosomes = []
    positions = []
    pvals = []
    perc_variances1 = []
    perc_variances2 = []
    h1_heritabilities = []
    mapped_tair_ids = []
#    chunk_i = 0
    for i, tair_id in enumerate(tair_ids):
        gd = gene_dict[tair_id]
        chrom = gd['chromosome']
        if chrom not in ['1', '2', '3', '4', '5']:
            continue
        chrom = int(chrom)
        start_pos = gd['start_pos'] - radius
        stop_pos = gd['end_pos'] + radius
        mean_pos = (start_pos + stop_pos) / 2
        d = sd.get_local_n_global_kinships(chrom=chrom, start_pos=start_pos, stop_pos=stop_pos,
                            global_kinship=K, kinship_method=kinship_method)
        if d['local_k'] != None and d['global_k'] != None:
            local_k = kinship.scale_k(d['local_k'])
            global_k = kinship.scale_k(d['global_k'])
            # print "Chromosome=%d, position=%d" % (chrom, focal_pos)
            res_dict = local_vs_global_mm(y, local_k, global_k, K, h0_res=h0_res)
            chromosomes.append(chrom)
            positions.append(mean_pos)
            perc_variances1.append(res_dict['perc_var1'])
            perc_variances2.append(res_dict['perc_var2'])
            h1_heritabilities.append(res_dict['pseudo_heritability1'])
            pvals.append(res_dict['pval'])
            mapped_tair_ids.append(tair_id)
            if plot_gene_trees and ets != None:
                tree_file = file_prefix + '_%s_%d_tree.pdf' % (tair_id, radius)
                y_strs = map(lambda x: '%0.2f' % x, y)
                snpsdata.plot_tree(local_k, tree_file, ets, verbose=True, label_values=y_strs)
                continue

            # print 'H0: pseudo_heritability=%0.2f' % (res_dict['h0_res']['pseudo_heritability'])
            # print 'H1: pseudo_heritability=%0.2f, perc_var1=%0.2f, perc_var2=%0.2f' % \
            #        (res_dict['h1_res']['pseudo_heritability'],
            #        res_dict['h1_res']['perc_var1'],
            #        res_dict['h1_res']['perc_var2'])
        if (i + 1) % int(len(tair_ids) / 100) == 0:  # Print dots
            sys.stdout.write('.')
            sys.stdout.flush()

    pval_res = gr.Result(scores=pvals, positions=positions, chromosomes=chromosomes)
    pval_res.neg_log_trans()
    pval_res.plot_manhattan(png_file=file_prefix + '_lrt_pvals.png', percentile=0, plot_bonferroni=True)
    perc_var_res = gr.Result(scores=perc_variances2, positions=positions, chromosomes=chromosomes)
    perc_var_res.plot_manhattan(png_file=file_prefix + '_perc_var_explained.png', percentile=0,
                ylab='% of variance explained')
    return {'pvals':pvals, 'perc_variances2':perc_variances2, 'perc_variances1':perc_variances1,
        'h0_heritability':h0_res['pseudo_heritability'], 'h1_heritabilities':h1_heritabilities,
        'chromosomes':chromosomes, 'positions':positions, 'tair_ids':mapped_tair_ids}



