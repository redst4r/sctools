import pandas as pd
import numpy as np
"""
some code around John Storeys Q-value concept
which seems alot more intuitive (but more mathematical/statistical) then
the BH-correction
"""

def _storey_estimate_pi0(pvals, lambda_):
    """
    estimate the fraction of null-hyoptheses (pi0) in the data.

    Relies on at "spike-and-slab"-like histogram of pvalues, where the high pvalues (e.g. >0.5)
    mostly come from the Null-Hyopthesis. 
    We pick all pvalues > lambda_, threat these as all being from H0 and extrapolate to the overall fraction of null-hypthoseses

    """
    m = len(pvals)
    W = np.sum(pvals > lambda_)
    pi0_hat = W/((1-lambda_) * m)
    return pi0_hat

def _storey_fdr(pvals, gamma, lambda_):
    """
    gamma: siginficant level, (usually alpha)
    lambda_: where to cutt of pvalues to be sure to have H0 tests (0.5 is a good value, anthingun >0.5 is probably a TN)
    """
    m = len(pvals)
    R = np.sum(pvals<= gamma)

    pi0_hat = _storey_estimate_pi0(pvals, lambda_)

    pr_pg = (R if R >0 else 1) / m

    pFDR = (pi0_hat*gamma) / ((pr_pg) * (1-(1-gamma)**m))

    FDR = (pi0_hat * gamma ) / pr_pg
    # print(pFDR, FDR)
    return pFDR, FDR

def storey_fdr(pvals, gamma, lambda_, n_bootstraps):
    """
    Estimate an upper bound of  FDR (and pFDR) of the pvalues IF thresholded at gamma.
    Lambda_ is needed to estimate the fraction of null-hypthoeseis

    # Parameters
    :param gamma: q-value threshold, e.g. 0.1 or 0.05
    :param labmda_: for estimating pi0, treat all p-values > lambda_ as beining from H0
    :param n_bootstraps: To get a handle on the FDR variability, perform bootstraps (and use the "biggest" FDR as an upper bound)
    """
    bootstrapped_pFDR = []
    bootstrapped_FDR = []
    for b in range(n_bootstraps):
        pvals_boot = np.random.choice(pvals, size=len(pvals), replace=True,)
        pFDR, FDR = _storey_fdr(pvals_boot, gamma, lambda_)
        bootstrapped_pFDR.append(pFDR)
        bootstrapped_FDR.append(FDR)

    pFDR_upper_bound = np.percentile(bootstrapped_pFDR, 95)
    FDR_upper_bound = np.percentile(bootstrapped_FDR, 95)

    return pFDR_upper_bound, FDR_upper_bound


def storey_qvalue(pvals, lambda_):
    """
    estimate the Q-values using Storey's technique (i.e. thresholding the q-values at \alpha guarantees a FDR <=
    alpha)

    # Parameters:
    :param pvals: list of UNCORRECTED pvalues. Should viusally follow a spike-and-lab-histogram!
    :param lambda_: to estimate pi0 (fraction of null-hypotheses) treat all p>lambda_ as coming from the H0
    """

    ix_sort = np.argsort(pvals)
    pvals_sorted=pvals[ix_sort]
    ix_undo_argsort = np.argsort(ix_sort)

    qvals = np.zeros_like(pvals)

    for (i, p) in reversed(list(enumerate(pvals_sorted))):
        # print(i)
        if i == len(qvals)-1:
            q, _ = _storey_fdr(pvals, p, lambda_)
            qvals[i] = q
        else:
            q, _ =  _storey_fdr(pvals, p, lambda_)
            qvals[i] = np.minimum(q, qvals[i+1])

    return qvals[ix_undo_argsort]

