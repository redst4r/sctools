"""
COPIED from scanpy.tl for now to fix the memory consumption.
also submitted as githb issue.

Calculate scores based on the expression of gene lists.
"""
from typing import Sequence, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.sparse import issparse

import logging as logg


def _sparse_nanmean(X, axis):
    """
    np.nanmean equivalent for sparse matrices

    # testing:
    >>> import numpy as np
    >>> from scipy.sparse import csr_matrix
    >>> R, C = 60, 50
    >>> A = np.random.rand(R,C) * (np.random.rand(R,C) < 0.3)
    >>> S = csr_matrix(A)

    # col sum
    >>> np.testing.assert_allclose(A.mean(0),np.array(_sparse_nanmean(S,0)).flatten())

    # rowsum
    >>> np.testing.assert_allclose(A.mean(1), np.array(_sparse_nanmean(S,1)).flatten())

    # now with nan
    >>> A = np.random.rand(R,C) * (np.random.rand(R,C) < 0.3)
    >>> masknan = (np.random.rand(R,C) < 0.3)
    >>> A[masknan] = np.nan

    >>> np.testing.assert_allclose(np.nanmean(A,1), np.array(_sparse_nanmean(csr_matrix(A),1)).flatten())
    >>> np.testing.assert_allclose(np.nanmean(A,0), np.array(_sparse_nanmean(csr_matrix(A),0)).flatten())

    # edge case of only NaNs per row
    >>> A = np.full((10,1), np.nan)
    >>> np.testing.assert_allclose(np.nanmean(A,0), np.array(_sparse_nanmean(csr_matrix(A),0)).flatten())
    """

    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    s = Y.sum(axis)
    m = s / n_elements

    return m


def score_genes(
    adata: AnnData,
    gene_list: Sequence[str],
    ctrl_size: int = 50,
    gene_pool: Optional[Sequence[str]] = None,
    n_bins: int = 25,
    score_name: str = 'score',
    random_state = 0,
    copy: bool = False,
    use_raw: bool = None,
) -> Optional[AnnData]:
    """\
    Score a set of genes [Satija15]_.

    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.

    This reproduces the approach in Seurat [Satija15]_ and has been implemented
    for Scanpy by Davide Cittaro.

    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
    ctrl_size
        Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    gene_pool
        Genes for sampling the reference set. Default is all genes.
    n_bins
        Number of expression level bins for sampling.
    score_name
        Name of the field to be added in `.obs`.
    random_state
        The random seed for sampling.
    copy
        Copy `adata` or modify it inplace.
    use_raw
        Use `raw` attribute of `adata` if present.

        .. versionchanged:: 1.4.5
           Default value changed from `False` to `None`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with an additional field
    `score_name`.

    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    start = logg.info(f'computing score {score_name!r}')
    adata = adata.copy() if copy else adata

    if random_state is not None:
        np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.raw.var_names if use_raw else adata.var_names
    genes_to_ignore = []
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
        else:
            genes_to_ignore.append(gene)
    if len(genes_to_ignore) > 0:
        logg.warning(f'genes are not in var_names and ignored: {genes_to_ignore}')
    gene_list = set(gene_list_in_var[:])

    if len(gene_list) == 0:
        logg.warning('provided gene list has length 0, scores as 0')
        adata.obs[score_name] = 0
        return adata if copy else None

    if gene_pool is None:
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    if use_raw is None:
        use_raw = True if adata.raw is not None else False
    _adata = adata.raw if use_raw else adata

    _adata_subset = _adata[:, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata
    if issparse(_adata_subset.X):
        obs_avg = pd.Series(
            np.array(_adata_subset.X.mean(axis=0)).flatten(), index=gene_pool)  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata_subset.X, axis=0), index=gene_pool)  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)] # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[list(gene_list)]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        # uses full r_genes if ctrl_size > len(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))

    # To index, we need a list – indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)

    X_list = _adata[:, gene_list].X
    if issparse(X_list):
        X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
    else:
        X_list = np.nanmean(X_list, axis=1)

    X_control = _adata[:, control_genes].X
    if issparse(X_control):
        X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        X_control = np.nanmean(X_control, axis=1)


    if len(gene_list) == 0:
        # We shouldn't even get here, but just in case
        logg.info(
            f'could not add \n'
            f'    {score_name!r}, score of gene set (adata.obs)'
        )
        return adata if copy else None
    elif len(gene_list) == 1:
        if _adata[:, gene_list].X.ndim == 2:
            vector = _adata[:, gene_list].X.toarray()[:, 0] # new anndata
        else:
            vector =  _adata[:, gene_list].X  # old anndata
        score = vector - X_control
    else:
        score = X_list - X_control

    adata.obs[score_name] = pd.Series(np.array(score).ravel(), index=adata.obs_names)

    # logg.info(
    #     '    finished',
    #     time=start,
    #     deep=(
    #         'added\n'
    #         f'    {score_name!r}, score of gene set (adata.obs).\n'
    #         f'    {len(control_genes)} total control genes are used.'
    #     ),
    # )
    return adata if copy else None


def score_genes_cell_cycle(
    adata: AnnData,
    s_genes: Sequence[str],
    g2m_genes: Sequence[str],
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """\
    Score cell cycle genes [Satija15]_.

    Given two lists of genes associated to S phase and G2M phase, calculates
    scores and assigns a cell cycle phase (G1, S or G2M). See
    :func:`~scanpy.tl.score_genes` for more explanation.

    Parameters
    ----------
    adata
        The annotated data matrix.
    s_genes
        List of genes associated with S phase.
    g2m_genes
        List of genes associated with G2M phase.
    copy
        Copy `adata` or modify it inplace.
    **kwargs
        Are passed to :func:`~scanpy.tl.score_genes`. `ctrl_size` is not
        possible, as it's set as `min(len(s_genes), len(g2m_genes))`.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **S_score** : `adata.obs`, dtype `object`
        The score for S phase for each cell.
    **G2M_score** : `adata.obs`, dtype `object`
        The score for G2M phase for each cell.
    **phase** : `adata.obs`, dtype `object`
        The cell cycle phase (`S`, `G2M` or `G1`) for each cell.

    See also
    --------
    score_genes

    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    logg.info('calculating cell cycle phase')

    adata = adata.copy() if copy else adata
    ctrl_size = min(len(s_genes), len(g2m_genes))
    # add s-score
    score_genes(adata, gene_list=s_genes, score_name='S_score', ctrl_size=ctrl_size, **kwargs)

    # add g2m-score
    score_genes(adata, gene_list=g2m_genes, score_name='G2M_score', ctrl_size=ctrl_size, **kwargs)
    scores = adata.obs[['S_score', 'G2M_score']]

    # default phase is S
    phase = pd.Series('S', index=scores.index)

    # if G2M is higher than S, it's G2M
    phase[scores.G2M_score > scores.S_score] = 'G2M'

    # if all scores are negative, it's G1...
    phase[np.all(scores < 0, axis=1)] = 'G1'

    adata.obs['phase'] = phase
    logg.info('    \'phase\', cell cycle phase (adata.obs)')
    return adata if copy else None



def my_ssgsea(x: pd.Series, omega: float, gs):
    """
    GSEA gene set scoring. See https://www.pathwaycommons.org/guide/primers/data_analysis/gsea/ for some details
    In brief:
    1. Rank genes
    2. go over the ranked genes, starting with the higest. Calculate a running sum, incrementing the sum
       each time we encounter a gene in the set, decrement the counter when we encouter a gene outside the set
    3. Aggregate the running sum: Maximum deviation from zero

    :param x: Rank-order statistic. Note that this IS NOT THE RANK itself. This is the quantity thats ranked. It is ALSO DIRECTLY ENTERNG the test ()!
    :param omega: exponent of the increments, see above link (its called alpha there)
    :param gs: Gene set, some iterable of gene names
    df
    """

    assert isinstance(x, pd.Series)
    gene_set = set(gs)

    # actually, get rid of genes that are not even in the data
    # those will screw up our Fi_not vector/values
    gene_set = gene_set & set(x.index)

    #first sort by absolute expression value, starting with the highest expressed genes first
    xsorted = x.sort_values(axis=0, ascending=False, inplace=False)
    keys_sorted = xsorted.index.tolist()

    # xsorted contains the "effect sizes", which should be a positive contributoin to the running sum.
    # Notice the |r| notation in the publications!
    xsorted = np.abs(xsorted)

    if False:
        ## seems slow, maybe because we dont iterate directly over keys_sorted
        nom_vector = np.zeros(len(x))
        for t in range(len(x)):
            # the_rank = len(x)-t   # actually the highest rank (the first one) is N, the lowest rank 1
            s = xsorted.iloc[t]
            nom_vector[t] = s**omega if keys_sorted[t] in gene_set else 0
    else:
        # muchfaster
        nom_vector = []
        for g in keys_sorted:
            if g in gene_set:
                s = xsorted.loc[g]
                increment = s**omega
            else:
                increment = 0
            nom_vector.append(increment)
    # Fi is the cumululative divided by the full sum
    Fi = np.cumsum(nom_vector)

    if Fi[-1] < 0:
        print(f"division by zero, somethings wrong here: {gene_set}: {Fi}")
    Fi = Fi / Fi[-1]

    # now for the term summarizing the genes NOT in the vector
    NO_vector = np.zeros(len(x))
    for t in range(len(x)):
        NO_vector[t] = 1 if not keys_sorted[t] in gene_set else 0
    Fi_not = np.cumsum(NO_vector) / (len(x) - len(gene_set))

    return max_deviation_from_zero(Fi - Fi_not)


def max_deviation_from_zero(x):
    """
    for vector x, get its biggest deviation from zero (inlcuding its sign)
    e.g:
    x = [0,1,-1, 2]  -> 2
    x = [0,1,-2, 1]  -> -2

    This is different than np.max(np.abs(x))!!
    """
    themax = np.max(x)
    themin = np.min(x)

    if np.abs(themax) > np.abs(themin):
        return themax
    else:
        return themin
