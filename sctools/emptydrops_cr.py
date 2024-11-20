# from sctools.emptydrops_cr import find_nonambient_barcodes, estimate_profile_sgt, est_background_profile_sgt
#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

""" Functions for calling cell-associated barcodes """
from collections import namedtuple
import numpy as np
import numpy.ma as ma
import scipy.stats as sp_stats
import tqdm

# from cellranger.analysis.diffexp import adjust_pvalue_bh
# import cellranger.stats as cr_stats
# import cellranger.sgt as cr_sgt

# Number of additional barcodes to consider after the initial cell calling
N_CANDIDATE_BARCODES=20000

# Number of partitions (max number of barcodes to consider for ambient estimation)
N_PARTITIONS=90000

# Drop this top fraction of the barcodes when estimating ambient.
MAX_OCCUPIED_PARTITIONS_FRAC = 0.5

# Minimum number of UMIS per barcode to consider after the initial cell calling
MIN_UMIS = 500

# Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling
MIN_UMI_FRAC_OF_MEDIAN = 0.01

# Maximum adjusted p-value to call a barcode as non-ambient
MAX_ADJ_PVALUE = 0.01

def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
    """ Estimate a gene expression profile by Simple Good Turing.
    Args:
      raw_mat (sparse matrix): Sparse matrix of all counts
      barcode_indices (np.array(int)): Barcode indices to use
      nz_feat (np.array(int)): Indices of features that are non-zero at least once
    Returns:
      profile (np.array(float)): Estimated probabilities of length len(nz_feat).
    """
    # Initial profile estimate
    prof_mat = matrix[:,barcode_indices]

    profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
    zero_feat = np.flatnonzero(profile == 0)

    # Simple Good Turing estimate
    p_smoothed, p0 = sgt_proportions(profile[np.flatnonzero(profile)])

    # Distribute p0 equally among the zero elements.
    p0_i = p0/len(zero_feat)

    profile_p = np.repeat(p0_i, len(nz_feat))
    profile_p[np.flatnonzero(profile)] = p_smoothed

    assert np.isclose(profile_p.sum(), 1.0)
    return profile_p


# Construct a background expression profile from barcodes with <= T UMIs
def est_background_profile_sgt(matrix, use_bcs):
    """ Estimate a gene expression profile on a given subset of barcodes.
         Use Good-Turing to smooth the estimated profile.
    Args:
      matrix (scipy.sparse.csc_matrix): Sparse matrix of all counts
      use_bcs (np.array(int)): Indices of barcodes to use (col indices into matrix)
    Returns:
      profile (use_features, np.array(float)): Estimated probabilities of length use_features.
    """
    # Use features that are nonzero anywhere in the data
    use_feats = np.flatnonzero(np.asarray(matrix.sum(1)))

    # Estimate background profile
    bg_profile_p = estimate_profile_sgt(matrix, use_bcs, use_feats)

    return (use_feats, bg_profile_p)

NonAmbientBarcodeResult = namedtuple('NonAmbientBarcodeResult',
                                     ['eval_bcs',      # Candidate barcode indices (n)
                                      'log_likelihood',# Ambient log likelihoods (n)
                                      'pvalues',       # pvalues (n)
                                      'pvalues_adj',   # B-H adjusted pvalues (n)
                                      'is_nonambient', # Boolean nonambient calls (n)
                                      ])

def find_nonambient_barcodes(adata, orig_cell_bcs,
                             min_umi_frac_of_median=MIN_UMI_FRAC_OF_MEDIAN,
                             min_umis_nonambient=MIN_UMIS,
                             max_adj_pvalue=MAX_ADJ_PVALUE, n_cores=0):
    """ Call barcodes as being sufficiently distinct from the ambient profile

    Args:
      adata: Full expression matrix.
      orig_cell_bcs (iterable of str): Strings of initially-called cell barcodes.
    Returns:
    TBD
    """
    # Estimate an ambient RNA profile
    umis_per_bc = np.array(adata.X.sum(1)).flatten()
    bc_order = np.argsort(umis_per_bc)

    # Take what we expect to be the barcodes associated w/ empty partitions.
    if False:
        empty_bcs = bc_order[::-1][int(N_PARTITIONS/2):N_PARTITIONS]  # this is weird, they simply consider barcodes 45k-end for this?!
        empty_bcs.sort()
    else:
        import warnings
        warnings.warn('changed the code here!!')
        empty_bcs = np.where(umis_per_bc < min_umis_nonambient)[0]
        empty_bcs.sort()
    # Require non-zero barcodes (i.e. empty droplets, but >0 counts)
    nz_bcs = np.flatnonzero(umis_per_bc)
    nz_bcs.sort()

    use_bcs = np.intersect1d(empty_bcs, nz_bcs, assume_unique=True)

    if len(use_bcs) > 0:
        try:
            print(f'Estimating background profile with {len(use_bcs)} barcodes')
            # need to convert to int!! otherwise downstream exception complaining that we cant cast float32 to int64
            eval_features, ambient_profile_p = est_background_profile_sgt(adata.X.T.astype(int), use_bcs)  #  transpose,, it expects BCs per column
        except SimpleGoodTuringError as e:
            print(str(e))
            return None
    else:
        eval_features = np.zeros(0, dtype=int)
        ambient_profile_p = np.zeros(0)

    # Choose candidate cell barcodes
    orig_cell_bc_set = set(orig_cell_bcs)
    orig_cells = np.flatnonzero(np.fromiter((bc in orig_cell_bc_set for bc in adata.obs.index),
                                            count=adata.shape[0], dtype=bool))

    # No good incoming cell calls
    if orig_cells.sum() == 0:
        print('orig_cells.sum()==0, returning None')
        return None

    # Look at non-cell barcodes above a minimum UMI count
    eval_bcs = np.ma.array(np.arange(adata.shape[0]))
    eval_bcs[orig_cells] = ma.masked  # masking everthing thats clearly a cell

    median_initial_umis = np.median(umis_per_bc[orig_cells])
    min_umis = int(max(min_umis_nonambient, round(np.ceil(median_initial_umis * min_umi_frac_of_median))))
    print('Median UMIs of initial cell calls: {}'.format(median_initial_umis))
    print('Min UMIs: {}'.format(min_umis))

    eval_bcs[umis_per_bc < min_umis] = ma.masked
    n_unmasked_bcs = len(eval_bcs) - eval_bcs.mask.sum()

    # Take the top N_CANDIDATE_BARCODES by UMI count, of barcodes that pass the above criteria
    eval_bcs = np.argsort(ma.masked_array(umis_per_bc, mask=eval_bcs.mask))[0:n_unmasked_bcs][-N_CANDIDATE_BARCODES:]

    if len(eval_bcs) == 0:
        print('len(eval_bcs)==0, returning None')
        return None

    assert not np.any(np.isin(eval_bcs, orig_cells))
    print('Calculating the difference from ambience for:')
    print('Number of candidate bcs: {}'.format(len(eval_bcs)))
    print('Range candidate bc umis: {}, {}'.format(umis_per_bc[eval_bcs].min(), umis_per_bc[eval_bcs].max()))

    eval_mat = adata.X.T[eval_features, :][:, eval_bcs]

    if len(ambient_profile_p) == 0:
        obs_loglk = np.repeat(np.nan, len(eval_bcs))
        pvalues = np.repeat(1, len(eval_bcs))
        sim_loglk = np.repeat(np.nan, len(eval_bcs))
        print('ambient profile all zero, retunrin None')
        return None

    # Compute observed log-likelihood of barcodes being generated from ambient RNA
    obs_loglk = eval_multinomial_loglikelihoods(eval_mat, ambient_profile_p)

    # Simulate log likelihoods
    # this needs to be paralellized!
    # Addendum: For each umi_count (e.g. 500,...,1000), this simulates the logp that you'd get with this number of UMIs
    # again ->int cast
    distinct_ns, sim_loglk = simulate_multinomial_loglikelihoods(ambient_profile_p, umis_per_bc[eval_bcs].astype(int), num_sims=10000, verbose=True, n_cores=n_cores)
    # print(distinct_ns.shape, sim_loglk.shape)

    # Compute p-values
    # basically checks how far out the observed loglike is from the simulated null-distribution
    pvalues = compute_ambient_pvalues(umis_per_bc[eval_bcs], obs_loglk, distinct_ns, sim_loglk)

    pvalues_adj = adjust_pvalue_bh(pvalues)

    is_nonambient = pvalues_adj <= max_adj_pvalue

    return NonAmbientBarcodeResult(
        eval_bcs=eval_bcs,
        log_likelihood=obs_loglk,
        pvalues=pvalues,
        pvalues_adj=pvalues_adj,
        is_nonambient=is_nonambient,
    )

def eval_multinomial_loglikelihoods(matrix, profile_p, max_mem_gb=0.1):
    """Compute the multinomial log PMF for many barcodes
    Args:
      matrix (scipy.sparse.csc_matrix): Matrix of UMI counts (feature x barcode)
      profile_p (np.ndarray(float)): Multinomial probability vector
      max_mem_gb (float): Try to bound memory usage.
    Returns:
      log_likelihoods (np.ndarray(float)): Log-likelihood for each barcode
    """
    gb_per_bc = float(matrix.shape[0] * matrix.dtype.itemsize) / (1024**3)
    bcs_per_chunk = max(1, int(round(max_mem_gb/gb_per_bc)))
    num_bcs = matrix.shape[1]

    loglk = np.zeros(num_bcs)

    for chunk_start in range(0, num_bcs, bcs_per_chunk):
        chunk = slice(chunk_start, chunk_start+bcs_per_chunk)
        matrix_chunk = matrix[:,chunk].transpose().toarray()
        n = matrix_chunk.sum(1)
        loglk[chunk] = sp_stats.multinomial.logpmf(matrix_chunk, n, p=profile_p)
    return loglk



def _sim(profile_p, distinct_n, jump, n_sample_feature_block=1000000):
    log_profile_p = np.log(profile_p)
    loglk = np.zeros((len(distinct_n)), dtype=float)

    curr_counts = np.ravel(sp_stats.multinomial.rvs(distinct_n[0], profile_p, size=1))
    curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[0], p=profile_p)

    loglk[0] = curr_loglk

    # this is a reservoir of random genes that get incremeted
    sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
    k = 0

    for i in range(1, len(distinct_n)):
        step = distinct_n[i] - distinct_n[i-1]
        if step >= jump:
            # Instead of iterating for each n, sample the intermediate ns all at once
            curr_counts += np.ravel(sp_stats.multinomial.rvs(step, profile_p, size=1))
            curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[i], p=profile_p)
            assert not np.isnan(curr_loglk)
        else:
            # Iteratively sample between the two distinct values of n
            for n in range(distinct_n[i-1]+1, distinct_n[i]+1):
                j = sampled_features[k]
                k += 1
                if k >= n_sample_feature_block:
                    # Amortize this operation
                    sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
                    k = 0
                curr_counts[j] += 1
                curr_loglk += log_profile_p[j] + np.log(float(n)/curr_counts[j])

        loglk[i] = curr_loglk

    return loglk

import multiprocessing as mp

def simulate_multinomial_loglikelihoods(profile_p, umis_per_bc,
                                        num_sims=1000, jump=1000,
                                        n_sample_feature_block=1000000, verbose=False, n_cores=0):
    """Simulate draws from a multinomial distribution for various values of N.
       Uses the approximation from Lun et al. ( https://www.biorxiv.org/content/biorxiv/early/2018/04/04/234872.full.pdf )
    Args:
      profile_p (np.ndarray(float)): Probability of observing each feature.
      umis_per_bc (np.ndarray(int)): UMI counts per barcode (multinomial N).
      num_sims (int): Number of simulations per distinct N value.
      jump (int): Vectorize the sampling if the gap between two distinct Ns exceeds this.
      n_sample_feature_block (int): Vectorize this many feature samplings at a time.
    Returns:
      (distinct_ns (np.ndarray(int)), log_likelihoods (np.ndarray(float)):
      distinct_ns is an array containing the distinct N values that were simulated.
      log_likelihoods is a len(distinct_ns) x num_sims matrix containing the
        simulated log likelihoods.
    """
    distinct_n = np.flatnonzero(np.bincount(umis_per_bc))

    loglk = np.zeros((len(distinct_n), num_sims), dtype=float)
    num_all_n = np.max(distinct_n) - np.min(distinct_n)


    sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
    k = 0

    import tqdm
    if n_cores==0:
        for sim_idx in tqdm.trange(num_sims):
            loglk[:, sim_idx] = _sim(profile_p, distinct_n, jump, n_sample_feature_block)
    else:
        print('parallel')
        with mp.Pool(n_cores) as pool:
            res = pool.starmap(_sim, ((profile_p, distinct_n, jump, n_sample_feature_block) for _ in range(num_sims)))
            loglk = np.stack(res).T

    return distinct_n, loglk




def simulate_multinomial_loglikelihoods_old(profile_p, umis_per_bc,
                                        num_sims=1000, jump=1000,
                                        n_sample_feature_block=1000000, verbose=False):
    """Simulate draws from a multinomial distribution for various values of N.

       Uses the approximation from Lun et al. ( https://www.biorxiv.org/content/biorxiv/early/2018/04/04/234872.full.pdf )

    Args:
      profile_p (np.ndarray(float)): Probability of observing each feature.
      umis_per_bc (np.ndarray(int)): UMI counts per barcode (multinomial N).
      num_sims (int): Number of simulations per distinct N value.
      jump (int): Vectorize the sampling if the gap between two distinct Ns exceeds this.
      n_sample_feature_block (int): Vectorize this many feature samplings at a time.
    Returns:
      (distinct_ns (np.ndarray(int)), log_likelihoods (np.ndarray(float)):
      distinct_ns is an array containing the distinct N values that were simulated.
      log_likelihoods is a len(distinct_ns) x num_sims matrix containing the
        simulated log likelihoods.
    """
    distinct_n = np.flatnonzero(np.bincount(umis_per_bc))

    loglk = np.zeros((len(distinct_n), num_sims), dtype=float)
    num_all_n = np.max(distinct_n) - np.min(distinct_n)
    if verbose:
        print( 'Number of distinct N supplied: %d' % len(distinct_n))
        print( 'Range of N: %d' % num_all_n)
        print( 'Number of features: %d' % len(profile_p))

    sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
    k = 0

    log_profile_p = np.log(profile_p)

    for sim_idx in tqdm.trange(num_sims):
#         if verbose and sim_idx % 100 == 99:
#             sys.stdout.write('.')
#             sys.stdout.flush()
        curr_counts = np.ravel(sp_stats.multinomial.rvs(distinct_n[0], profile_p, size=1))

        curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[0], p=profile_p)

        loglk[0, sim_idx] = curr_loglk

        for i in range(1, len(distinct_n)):
            step = distinct_n[i] - distinct_n[i-1]
            if step >= jump:
                # Instead of iterating for each n, sample the intermediate ns all at once
                curr_counts += np.ravel(sp_stats.multinomial.rvs(step, profile_p, size=1))
                curr_loglk = sp_stats.multinomial.logpmf(curr_counts, distinct_n[i], p=profile_p)
                assert not np.isnan(curr_loglk)
            else:
                # Iteratively sample between the two distinct values of n
                for n in range(distinct_n[i-1]+1, distinct_n[i]+1):
                    j = sampled_features[k]
                    k += 1
                    if k >= n_sample_feature_block:
                        # Amortize this operation
                        sampled_features = np.random.choice(len(profile_p), size=n_sample_feature_block, p=profile_p, replace=True)
                        k = 0
                    curr_counts[j] += 1
                    curr_loglk += log_profile_p[j] + np.log(float(n)/curr_counts[j])

            loglk[i, sim_idx] = curr_loglk

#     if verbose:
#         sys.stdout.write('\n')

    return distinct_n, loglk

# alternatively
# from statsmodels.stats.multitest import multipletests
# multipletests(p, method='fdr_bh')
def adjust_pvalue_bh(p):
    """ Multiple testing correction of p-values using the Benjamini-Hochberg procedure """
    descending = np.argsort(p)[::-1]
    # q = p * N / k where p = p-value, N = # tests, k = p-value rank
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

    # Return to original order
    return q[np.argsort(descending)]

def compute_ambient_pvalues(umis_per_bc, obs_loglk, sim_n, sim_loglk):
    """Compute p-values for observed multinomial log-likelihoods
    Args:
      umis_per_bc (nd.array(int)): UMI counts per barcode
      obs_loglk (nd.array(float)): Observed log-likelihoods of each barcode deriving from an ambient profile
      sim_n (nd.array(int)): Multinomial N for simulated log-likelihoods
      sim_loglk (nd.array(float)): Simulated log-likelihoods of shape (len(sim_n), num_simulations)
    Returns:
      pvalues (nd.array(float)): p-values
    """
    assert len(umis_per_bc) == len(obs_loglk)
    assert sim_loglk.shape[0] == len(sim_n)

    # Find the index of the simulated N for each barcode
    sim_n_idx = np.searchsorted(sim_n, umis_per_bc)
    num_sims = sim_loglk.shape[1]

    num_barcodes = len(umis_per_bc)

    pvalues = np.zeros(num_barcodes)

    for i in range(num_barcodes):
        num_lower_loglk = np.sum(sim_loglk[sim_n_idx[i],:] < obs_loglk[i])
        pvalues[i] = float(1 + num_lower_loglk) / (1 + num_sims)
    return pvalues



#!/usr/bin/env python
#
# Copyright (c) 2018 10X Genomics, Inc. All rights reserved.
#

"""
Simple Good-Turing estimator.
Based on S implementation in
  William A. Gale & Geoffrey Sampson (1995) Good-turing frequency estimation without tears,
  Journal of Quantitative Linguistics, 2:3, 217-237, DOI: 10.1080/09296179508590051
"""



class SimpleGoodTuringError(Exception):
    pass

def _averaging_transform(r, nr):
    d = np.concatenate((np.ones(1, dtype=int), np.diff(r)))
    dr = np.concatenate((
        0.5 * (d[1:] + d[0:-1]),
        np.array((d[-1],), dtype=float),
        ))
    return nr.astype(float)/dr

def _rstest(r, coef):
    return r * np.power(1 + 1/r, 1 + coef)

def simple_good_turing(xr, xnr):
    """Make a Simple Good-Turing estimate of the frequencies.

    Args:
      xr (np.array(int)): Non-zero item frequencies
      xnr (np.array(int)): Non-zero frequencies of frequencies
    Returns:
      (rstar (np.array(float)), p0 (float)):
        rstar: The adjusted non-zero frequencies
        p0: The total probability of unobserved items
    """

    xr = xr.astype(float)
    xnr = xnr.astype(float)

    xN = np.sum(xr*xnr)

    # Get Linear Good-Turing estimate
    xnrz = _averaging_transform(xr, xnr)
    slope, intercept, _, _, _ = sp_stats.linregress(np.log(xr), np.log(xnrz))

    if slope > -1:
        raise SimpleGoodTuringError("The log-log slope is > -1 (%d); the SGT estimator is not applicable to these data." % slope)

    xrst = _rstest(xr,slope)
    xrstrel = xrst/xr

    # Get traditional Good-Turing estimate
    xrtry = xr == np.concatenate((xr[1:]-1, np.zeros(1)))
    xrstarel = np.zeros(len(xr))
    xrstarel[xrtry] = (xr[xrtry]+1) / xr[xrtry] * \
                      np.concatenate((xnr[1:], np.zeros(1)))[xrtry] / xnr[xrtry]

    # Determine when to switch from GT to LGT estimates
    tursd = np.ones(len(xr))
    for i in range(len(xr)):
        if xrtry[i]:
            tursd[i] = float(i+2) / xnr[i] * np.sqrt(xnr[i+1] * (1 + xnr[i+1]/xnr[i]))

    xrstcmbrel = np.zeros(len(xr))
    useturing = True
    for r in range(len(xr)):
        if not useturing:
            xrstcmbrel[r]  = xrstrel[r]
        else:
            if np.abs(xrstrel[r]-xrstarel[r]) * (1+r)/tursd[r] > 1.65:
                xrstcmbrel[r] = xrstarel[r]
            else:
                useturing = False
                xrstcmbrel[r] = xrstrel[r]

    # Renormalize the probabilities for observed objects
    sumpraw = np.sum(xrstcmbrel * xr * xnr / xN)

    xrstcmbrel = xrstcmbrel * (1 - xnr[0] / xN) / sumpraw
    p0 = xnr[0]/xN

    return (xr * xrstcmbrel, p0)

def sgt_proportions(frequencies):
    """Use Simple Good-Turing estimate to adjust for unobserved items

    Args:
      frequencies (np.array(int)): Nonzero frequencies of items
    Returns:
      (pstar (np.array(float)), p0 (float)):
        pstar: The adjusted non-zero proportions
        p0: The total probability of unobserved items
    """
    if len(frequencies) == 0:
        raise ValueError("Input frequency vector is empty")
    if np.count_nonzero(frequencies) != len(frequencies):
        raise ValueError("Frequencies must be greater than zero")

    freqfreqs = np.bincount(frequencies)
    assert freqfreqs[0] == 0
    use_freqs = np.flatnonzero(freqfreqs)

    if len(use_freqs) < 10:
        raise SimpleGoodTuringError("Too few non-zero frequency items (%d). Aborting SGT." % len(use_freqs))


    rstar, p0 = simple_good_turing(use_freqs, freqfreqs[use_freqs])

    # rstar contains the smoothed frequencies.
    # Map each original frequency r to its smoothed rstar.
    rstar_dict = dict(zip(use_freqs, rstar))

    rstar_sum = np.sum(freqfreqs[use_freqs] * rstar)
    rstar_i = np.fromiter((rstar_dict[f] for f in frequencies),
                          dtype=float, count=len(frequencies))
    pstar = (1 - p0) * (rstar_i / rstar_sum)

    assert np.isclose(p0 + np.sum(pstar), 1)
    return (pstar, p0)
