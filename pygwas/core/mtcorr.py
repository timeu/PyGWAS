"""
Multiple test corrections
"""
import scipy as sp


def get_bh_thres(pvals, fdr_thres=0.05):
    """
    Implements Benjamini-Hochberg FDR threshold (1995)
    """
    m = len(pvals)   
    s_pvals = sorted(pvals) 
    for i, p in enumerate(s_pvals):
        thes_pval = ((i + 1.0) / float(m)) * fdr_thres
        if p > thes_pval:
            break
        
    return {'thes_pval':thes_pval, 'thres_i':i}
        
        
def get_bhy_thres(pvals, fdr_thres=0.05):
    """
    Implements the Benjamini-Hochberg-Yekutieli procedure (2001).
    Assumes arbitrary dependence between variables.
    """
    #Euler-Mascheroni constant 
    gamma = 0.57721566490153286060651209008    
    m = len(pvals)   
    m_float = float(m) 
    s_pvals = sorted(pvals) 
    s = 1.0
    for i, p in enumerate(s_pvals):
        if i > 2:
            s = s + 1.0/(i-1)
        thes_pval = ((i + 1.0) / m_float) * fdr_thres / s
        if p > thes_pval:
            break        
    return {'thes_pval':thes_pval, 'thres_i':i}

        
def bh_test():
    import random
    pvals = [random.random() for i in range(250000)]
    print 'BHY', get_bhy_thres(pvals, fdr_thres=0.05)
    print 'BH', get_bh_thres(pvals, fdr_thres=0.05)


if __name__ == '__main__':
    for i in range(100):
        bh_test()
