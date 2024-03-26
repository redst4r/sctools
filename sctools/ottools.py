"""
docstrign
"""

import ot
import tqdm
from sklearn.manifold import MDS
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform, pdist
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import numpy as np
import gc

def get_filtered_matrix(adata, matrix, g1, g2, group_name):
    """
    filters down the matrix according to the groups in adata to
    only contain cells from g1 in the rows and g2 in the colums.
    Basically creats the costmatrix for OT compring group1 to group2
    """
    ix1 = np.where(adata.obs[group_name]==g1)[0]
    ix2 = np.where(adata.obs[group_name]==g2)[0]

    filtered_D = matrix[np.ix_(ix1, ix2)] # pick only the relevant datapoints
    i1 = adata[ix1].obs.index.values
    i2 = adata[ix2].obs.index.values
    return filtered_D, i1, i2

def yield_pairs(adata, cost_matrix, group_name):
    """
    iterating over all pairs of groupings, returning the respective cost_matrix
    """
    groups = np.sort(np.unique(adata.obs[group_name].values))

    print(f'Comparing {len(groups)} groups')
    for i, g1 in enumerate(groups):
        for j, g2 in enumerate(groups):
            if i>=j: continue
            filtered_D, _, _ = get_filtered_matrix(adata, cost_matrix, g1, g2, group_name)
            yield g1, g2, filtered_D


def calc_patient_matrix(adata, cost_matrix, group_name, ot_method, ot_params=None):
    """
    main OT function: given an attribute [in .obs] calcualte OT beetween cells from those groups

    Parameters:
    :param adata: AnnData containing the scRNAseq data
    :param cost_matrix: The cost-matrix of transporting each cell to each other. Must be adata.shape[0] x adata.shape[0]
    :param group_name: A field in adata.obs[] used to define the OT groups to compare
    :param ot_method: ['exact', 'entropy']
    """

    if ot_params is None:
        ot_params = {}

    assert ot_method in ['exact', 'entropy']
    if ot_method == 'entropy':
        assert 'lambda' in ot_params

    dmat_emd_df = []

    def ot_function(mat):
        if ot_method =='exact':
            transport_plan = ot.emd([], [], mat, numItermax=500_000)
        else:
            # in ot <0.8 theres a bug where sinkhorn2 still returns the transport map
            # they fixed in in 0.8
            # to cirvumvent this issue, just use sinkhorn() and compute the cost ourselves
            #
            # we rescale mat by its maximum here for numeric stability; doesnt affect teh solutiuon
            transport_plan = ot.sinkhorn([], [], mat/mat.max(), reg=ot_params['lambda'], numItermax=5_000)

        return np.sum(mat*transport_plan), transport_plan

    debias_ot = {}  # store the OT(a,a) in here
    for g in tqdm.tqdm(adata.obs[group_name].unique(), desc='debiasing terms'):
        filtered_D, _, _ = get_filtered_matrix(adata, cost_matrix, g, g, group_name)
        debias_ot[g], _ = ot_function(filtered_D)
        dmat_emd_df.append({'group1': g, 'group2': g, 'distance': debias_ot[g], 'debiased_distance': 0, 'upper_tria': 'no'})

    transport_plans = {}
    for g1, g2, filtered_D in tqdm.tqdm(yield_pairs(adata, cost_matrix, group_name), desc='OT pairs'):
        wd2, transport_plans[(g1,g2)] = ot_function(filtered_D)
        wd2_debiased = wd2 - 0.5 * debias_ot[g1] - 0.5 * debias_ot[g2]
        dmat_emd_df.append({'group1': g1, 'group2': g2, 'distance': wd2, 'debiased_distance':wd2_debiased, 'upper_tria': 'yes'})
        dmat_emd_df.append({'group1': g2, 'group2': g1, 'distance': wd2, 'debiased_distance':wd2_debiased, 'upper_tria': 'no'})  # symmetric
    return pd.DataFrame(dmat_emd_df) , transport_plans


def yield_pairs_compute_cost(adata, cost_matrix_fn, group_name):
    """
    computing the cost matrix, instead of just subsetting the HUGE matrix
    """
    groups = np.sort(np.unique(adata.obs[group_name].values))

    print(f'Comparing {len(groups)} groups')
    for i, g1 in enumerate(groups):
        for j, g2 in enumerate(groups):
            if i>=j: continue

            ix1 = np.where(adata.obs[group_name]==g1)[0]
            ix2 = np.where(adata.obs[group_name]==g2)[0]

            filtered_D = cost_matrix_fn(adata[ix1], adata[ix2] )
            yield g1, g2, filtered_D


def calc_patient_matrix_cost_matrix_on_the_fly(adata, cost_matrix_fn, group_name, ot_method, ot_params=None, verbose=False):
    """
    this calculates the cost matrix separatley for each pair, instead of one HUGE matrix and indexing into it
    This function is usefull for LARGE datasets. Storing the transport plans (a series of n1 x n2 matrixes, n1, n2 are groupsizes) gets expensive, so we drop those

    :param adata: scRNAseq single cell data
    :param cost_matrix_fn: how to calculate the cost amtrix between two groups. must be of signature f(AnnData, AnnData), producing a cost matrix between the two data groups
    ;param group_name: field in adata.obs to indentify the OT groups
    :param ot_method: ['exact', 'entropy']


    :returns:
    - The OT distances in a flat dataframe
    """
    if ot_params is None:
        ot_params = {}

    assert ot_method in ['exact', 'entropy']
    if ot_method == 'entropy':
        assert 'lambda' in ot_params

    dmat_emd_df = []

    def ot_function(mat):
        if ot_method == 'exact':
            transport_plan = ot.emd([], [], mat, numItermax=500_000)
        else:
            # we rescale mat by its maximum here for numeric stability; doesnt affect teh solutiuon
            transport_plan = ot.sinkhorn([], [], mat/mat.max(), reg=ot_params['lambda'], numItermax=5_000)
        return np.sum(mat*transport_plan), transport_plan

    debias_ot = {}  # store the OT(a,a) in here
    for g in tqdm.tqdm(adata.obs[group_name].unique(), desc='debiasing terms'):
        a1 = adata[adata.obs[group_name] == g]
        if verbose:
            print(f'Comp. distance matrix: {g} vs {g}: {len(a1)}x{len(a1)}')
        filtered_D = cost_matrix_fn(a1, a1)
        if verbose:
            print(f'Comp. OT: {g} vs {g}: {len(a1)}x{len(a1)}')

        debias_ot[g], _ = ot_function(filtered_D)
        dmat_emd_df.append({'group1': g, 'group2': g, 'distance': debias_ot[g], 'debiased_distance': 0 , 'upper_tria': 'no'})
        gc.collect()

    # transport_plans = {}
    n_groups = len(adata.obs[group_name].unique())
    n_comparisions = int(n_groups * (n_groups -1 ) / 2 )
    for g1, g2, filtered_D in tqdm.tqdm(yield_pairs_compute_cost(adata, cost_matrix_fn, group_name), desc='OT pairs', total=n_comparisions):
        if verbose:
            print(f'{g1} vs {g2}: {filtered_D.shape}')
        # wd2, transport_plans[(g1,g2)] = ot_function(filtered_D)
        wd2, _ = ot_function(filtered_D)
        if verbose:
            print('done')
        wd2_debiased = wd2 - 0.5 * debias_ot[g1] - 0.5 * debias_ot[g2]
        dmat_emd_df.append({'group1': g1, 'group2': g2, 'distance': wd2, 'debiased_distance':wd2_debiased, 'upper_tria': 'yes'})
        dmat_emd_df.append({'group1': g2, 'group2': g1, 'distance': wd2, 'debiased_distance':wd2_debiased, 'upper_tria': 'no'})  # symmetric
        gc.collect()
    return pd.DataFrame(dmat_emd_df)


def calc_patient_matrix_old(adata, cost_matrix, group_name, ot_method, ot_params=None):

    assert ot_method in ['exact', 'entropy']
    if ot_method == 'entropy':
        assert 'lambda' in ot_params

    groups = np.sort(np.unique(adata.obs[group_name].values))
    dmat_emd = np.zeros([len(groups), len(groups)])

    print(f'Comparing {len(groups)} groups')
    for i, g1 in tqdm.tqdm(enumerate(groups)):
        for j, g2 in enumerate(groups):
            if i>=j: continue
            ix1 = np.where(adata.obs[group_name]==g1)[0]
            ix2 = np.where(adata.obs[group_name]==g2)[0]

            filtered_D = cost_matrix[np.ix_(ix1, ix2)] # pick only the relevant datapoints

            if ot_method =='exact':
                wd2 = ot.emd2([],[],filtered_D, numItermax=500000)
            elif ot_method == 'entropy':
                wd2 = ot.sinkhorn2([],[],filtered_D, reg=ot_params['lambda'] ,numItermax=5000)
            else:
                raise ValueError('unnknonwn OT method')
#             wd2 = ot.bregman.sinkhorn_stabilized([],[],filtered_D, reg=0.01)

            dmat_emd[i,j] = wd2
            dmat_emd[j,i] = wd2

    return dmat_emd, groups


def flatdf_to_distance(df, distance_field):
    """
    turing a dataframe representign a flat distance matrix into a square distance matrix
    """
    def _aggr(x):
        assert len(x) == 1
        return x

    distance_matrix = pd.crosstab(df.group1, df.group2, df[distance_field], aggfunc=_aggr)
    return distance_matrix


def display_distance_matrix(df, method='average', figsize=(5,5), distance_field='debiased_distance'):
    """
    clustermap would cluster on the distance-matrix, instead taking the distance-matrix as an already precomputed clustering
    hence we explicitly tell it to no cluster by itself
    """
    distance_matrix = flatdf_to_distance(df, distance_field)

    flat_dmat = squareform(distance_matrix) # weird, squareform is its own inverse
    row_linkage = hc.linkage(flat_dmat, method=method)
    df_ = distance_matrix
    g = sns.clustermap(df_, row_linkage=row_linkage, col_linkage=row_linkage, annot=True, figsize=figsize)


def display_distance_matrix_old(distance_matrix, sample_names, method='average'):
    """
    clustermap would cluster on the distance-matrix, instead taking the distance-matrix as an already precomputed clustering
    hence we explicitly tell it to no cluster by itself
    """
    assert distance_matrix.shape[0] == distance_matrix.shape[1]
    flat_dmat = squareform(distance_matrix) # weird, squareform is its own inverse
    row_linkage = hc.linkage(flat_dmat, method=method)
    df_ = pd.DataFrame(distance_matrix, columns=sample_names, index=sample_names)
    g = sns.clustermap(df_, row_linkage=row_linkage, col_linkage=row_linkage, annot=True, figsize=(5,5))
#     plt.figure()
#     g = sns.heatmap(df_, annot=True)

def mds_embedding(df_flat, distance_field='debiased_distance'):

    distance_matrix = flatdf_to_distance(df_flat, distance_field)
    mds = MDS(n_components=2, metric=True, dissimilarity='precomputed', n_init=200)
    mds_emb = mds.fit_transform(distance_matrix)

    # lets plot it a little nicer
    df = pd.DataFrame(mds_emb, columns=['mds1', 'mds2'])
    df['samplename'] = distance_matrix.index

    plt.figure(figsize=(10,5))
    plt.subplot(121)
    sns.scatterplot(x='mds1', y='mds2', hue='samplename', data=df, legend=False)
    plt.axis('equal')
    plt.subplot(122)
    sns.scatterplot(x='mds1', y='mds2', hue='samplename', data=df)
    plt.legend()
    plt.axis('equal')

    return mds_emb

def _get_edges_from_transport_plan(transport_plan, n_lines):
    """
    note: this gets only the most important (weighted) directed edges from SOURCE->TARGET
    to get the other way around, just run on the transpose
    """
    # lets collect all edges that we want to draw in a mask/adjacency matrix first
    adj_mat = np.zeros_like(transport_plan)
    for i in tqdm.trange(transport_plan.shape[0]):
        conditional_map = transport_plan[i]/transport_plan[i].sum()  # where does datapoint i get transformed to
        pdf = np.sort(conditional_map)[::-1]
#         cdf =  np.cumsum(pdf)
#         cutoff = np.min(pdf[cdf <= 0.95])
        cutoff = pdf[np.minimum(n_lines, len(pdf))-1]
        ix = conditional_map > cutoff
        adj_mat[i,:] = ix
    return adj_mat
        # ix2plot = np.where(ix)[0]

def visualize_transport(X1, X2, gamma, n_lines, size=5, background_X=None, cost_matrix=None, linewidth=1, linealpha=0.1, figsize=(10,10)):
    """
    cost_matrix: this would be a cell-by-cell distance matrix, used to color the edges according to the cost of moving source to target
    n_lines: for each datapoint in the source distriubtion, plot only the first n_lines targets (with strongest association)
    """
    import matplotlib.pyplot as plt
    from matplotlib import collections  as mc

    assert len(X1) == gamma.shape[0]
    assert len(X2) == gamma.shape[1]

    fig, ax = plt.subplots(figsize=figsize)
    if background_X is not None:
        plt.scatter(background_X[:, 0], background_X[:, 1], s=size, alpha=0.4, c='grey')

    plt.scatter(X1[:, 0], X1[:, 1], s=size)
    plt.scatter(X2[:, 0], X2[:, 1], s=size)

    # lets collect all edges that we want to draw in a mask/adjacency matrix first
    adj_mat1 = _get_edges_from_transport_plan(gamma, n_lines)
    adj_mat2 = _get_edges_from_transport_plan(gamma.T, n_lines).T
    print('n_entries1', np.sum(adj_mat1 > 0))
    print('n_entries2', np.sum(adj_mat2 > 0))

    adj_mat = adj_mat1 + adj_mat2
    print('n_entries', np.sum(adj_mat > 0))


    if cost_matrix is not None:
        emin = cost_matrix[adj_mat != 0].min()
        emax = cost_matrix[adj_mat != 0].max()
        edgecolor_matrix = cost_matrix * adj_mat - emin
        edgecolor_matrix /= emax

    colors = []
    lines = []

    for i in tqdm.trange(len(X1)):
        for j in range(len(X2)):
            if adj_mat[i, j] != 0:
                if cost_matrix is not None:
                    col_ix = edgecolor_matrix[i,j]
                    col = plt.cm.viridis(col)
                else:
                    col='black'
                # print(col)
                colors.append(col)
                lines.append(
                    [(X1[i,0], X1[i,1]),
                     (X2[j,0], X2[j,1])
                    ]
                )
    print('#lines:', len(lines))

    lc = mc.LineCollection(lines, colors=colors, linewidths=linewidth, alpha=linealpha)
    ax.add_collection(lc)
    plt.colorbar()
    plt.show()
