"""
A container for functions which aim to analyze or process gwas results, for some aim.
"""

import scipy as sp
from scipy import stats
import logging
import math

log = logging.getLogger(__name__)






def calc_median(scores,exp_median=None):
    s = sp.copy(scores)
    s.sort()
    median = s[len(s) / 2]
    del s
    return median

def _estAreaBetweenCurves_(quantiles, expQuantiles):
    area = 0
    for i in range(0, len(quantiles) - 1):
        area += (expQuantiles[i + 1] - expQuantiles[i]) * (abs(quantiles[i + 1] - expQuantiles[i + 1] + quantiles[i] - expQuantiles[i])) / 2.0
    return area

def calc_ks_stats(scores, exp_scores=None):

    if exp_scores:
        (D, p_val) = stats.ks_2samp(scores, exp_scores)
    else:
        (D, p_val) = stats.kstest(scores, stats.uniform.cdf)
    return {'D':D, 'p_val':p_val}

def _getExpectedPvalueQuantiles_(numQuantiles):
    quantiles = []
    for i in range(numQuantiles):
        quantiles.append(float(i) + 0.5 / (numQuantiles + 1))
    return quantiles


def calculate_sp_pval(phen_vals):
    r = stats.shapiro(phen_vals)
    if sp.isfinite(r[0]):
        sp_pval = r[1]
    else:
        sp_pval = 0.0
    return sp_pval


def get_log_quantiles(scores, num_dots=1000, max_val=5):
    """
    Uses scipy
    """
    scores = sp.copy(sp.array(scores))
    scores.sort()
    indices = sp.array(10 ** ((-sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1)) * max_val) \
                * len(scores), dtype='int')
    return -sp.log10(scores[indices])




def get_quantiles(scores, num_dots=1000):
    """
    Uses scipy
    """
    scores = sp.copy(sp.array(scores))
    scores.sort()
    indices = [int(len(scores) * i / (num_dots + 2)) for i in range(1, num_dots + 1)]
    return scores[indices]



def calculate_qqplot_data(pvals,num_dots=1000):
    max_val = -math.log10(min(pvals))
    quantiles = get_quantiles(pvals, num_dots=num_dots)
    exp_quantiles = _getExpectedPvalueQuantiles_(num_dots)
    log_quantiles = get_log_quantiles(pvals, num_dots=num_dots, max_val=max_val)
    exp_log_quantiles = sp.arange(1, num_dots + 1, dtype='single') / (num_dots + 1) * max_val

    quantiles_dict = {'quantiles':quantiles, 'exp_quantiles':exp_quantiles,
            'log_quantiles':log_quantiles, 'exp_log_quantiles':exp_log_quantiles}
    return quantiles_dict






