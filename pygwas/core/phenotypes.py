   def get_correlations(self, pids=None):
        """
        Returns correlation matrix between traits
        
        All traits are used if pids is left empty.
        """
        import bisect
        if not pids:
            pids = sorted(self.phen_dict.keys())

        num_traits = len(pids)
        corr_mat = sp.ones((num_traits, num_traits))
        for i, pid1 in enumerate(pids):
            pd = self.get_avg_value_dict(pid1)
            ets1 = pd['ecotypes']
            pvs1 = pd['values']
            for j, pid2 in enumerate(pids[:i]):
                pd = self.get_avg_value_dict(pid2)
                ets2 = pd['ecotypes']
                pvs2 = pd['values']
                common_ets = set(ets1).intersection(set(ets2))
                ets_ix1 = map(ets1.index, common_ets)
                ets_ix2 = map(ets2.index, common_ets)
                vs1 = [pvs1[et_i] for et_i in ets_ix1]
                vs2 = [pvs2[et_i] for et_i in ets_ix2]
                corr_mat[i, j] = sp.corrcoef(vs1, vs2)[0, 1]
                corr_mat[j, i] = corr_mat[i, j]
        return corr_mat, pids

    def get_correlation(self, pid1, phed, pid2):
        """
        Returns the correlation with another phenotype_data object 
        """
        assert pid1 in self.phen_dict, 'phenotype ID %d missing in the self phed??' % pid1
        assert pid2 in phed.phen_dict, 'phenotype ID %d missing in the self phed??' % pid2
        pd = self.get_avg_value_dict(pid1)
        ets1 = pd['ecotypes']
        pvs1 = pd['values']
        pd = phed.get_avg_value_dict(pid2)
        ets2 = pd['ecotypes']
        pvs2 = pd['values']
        common_ets = set(ets1).intersection(set(ets2))
        ets_ix1 = map(ets1.index, common_ets)
        ets_ix2 = map(ets2.index, common_ets)
        vs1 = [pvs1[et_i] for et_i in ets_ix1]
        vs2 = [pvs2[et_i] for et_i in ets_ix2]
        return sp.corrcoef(vs1, vs2)[0, 1]

