import tqdm
from sklearn.manifold import MDS
import seaborn as sns
import ot
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform, pdist
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import numpy as np

def get_filtered_matrix(adata, matrix, g1, g2, group_name):
    """
    filters down the matrix according to the groups in adata to
    only contain cells from g1 in the rows and g2 in the colums.
    Basically creats the costmatrix for OT compring group1 to group2
    """
    groups = np.sort(np.unique(adata.obs[group_name].values))
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
            filtered_D,_,_ = get_filtered_matrix(adata, cost_matrix, g1, g2, group_name)
            yield g1, g2, filtered_D


def calc_patient_matrix(adata, cost_matrix, group_name, ot_method, ot_params=None):
    """
    main OT function: given an attribute [in .obs] calcualte OT beetween cells from those groups
    """

    if ot_params is None:
        ot_params = {}

    assert ot_method in ['exact', 'entropy']
    if ot_method == 'entropy':
        assert 'lambda' in ot_params

    dmat_emd_df = []

    def ot_function(mat):
        if ot_method =='exact':
            return ot.emd2([],[],mat, numItermax=500_000)
        else:
            return ot.sinkhorn2([],[],mat, reg=ot_params['lambda'] ,numItermax=5_000)

    debias_ot = {}  # store the OT(a,a) in here
    for g in tqdm.tqdm(adata.obs[group_name].unique(), desc='debiasing terms'):
        filtered_D, _, _ = get_filtered_matrix(adata, cost_matrix, g, g, group_name)
        debias_ot[g] = ot_function(filtered_D)
        dmat_emd_df.append({'group1': g, 'group2': g, 'distance': debias_ot[g], 'debiased_distance': 0 , 'upper_tria': 'no'})


    for g1, g2, filtered_D in tqdm.tqdm(yield_pairs(adata, cost_matrix, group_name), desc='OT pairs'):
        wd2 = ot_function(filtered_D)
        wd2_debiased = wd2 - 0.5 * debias_ot[g1] - 0.5 * debias_ot[g2]
        dmat_emd_df.append({'group1': g1, 'group2': g2, 'distance': wd2, 'debiased_distance':wd2_debiased, 'upper_tria': 'yes'})
        dmat_emd_df.append({'group1': g2, 'group2': g1, 'distance': wd2, 'debiased_distance':wd2_debiased, 'upper_tria': 'no'})  # symmetric
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


def flatdf_to_distance(df):
    """
    turing a dataframe representign a flat distance matrix into a square distance matrix
    """
    def _aggr(x):
        assert len(x) == 1
        return x

    distance_matrix = pd.crosstab(df.group1, df.group2, df.debiased_distance, aggfunc=_aggr)
    return distance_matrix


def display_distance_matrix(df, method='average'):
    "clustermap would cluster on the distance-matrix, instead taking the distance-matrix as an already precomputed clustering"
    "hence we explicitly tell it to no cluster by itself"
 
    distance_matrix = flatdf_to_distance(df)

    flat_dmat = squareform(distance_matrix) # weird, squareform is its own inverse
    row_linkage = hc.linkage(flat_dmat, method=method)
    df_ = distance_matrix
    g = sns.clustermap(df_, row_linkage=row_linkage, col_linkage=row_linkage, annot=True, figsize=(5,5))
#     plt.figure()
#     g = sns.heatmap(df_, annot=True)

def display_distance_matrix_old(distance_matrix, sample_names, method='average'):
    "clustermap would cluster on the distance-matrix, instead taking the distance-matrix as an already precomputed clustering"
    "hence we explicitly tell it to no cluster by itself"
    assert distance_matrix.shape[0] == distance_matrix.shape[1]
    flat_dmat = squareform(distance_matrix) # weird, squareform is its own inverse
    row_linkage = hc.linkage(flat_dmat, method=method)
    df_ = pd.DataFrame(distance_matrix, columns=sample_names, index=sample_names)
    g = sns.clustermap(df_, row_linkage=row_linkage, col_linkage=row_linkage, annot=True, figsize=(5,5))
#     plt.figure()
#     g = sns.heatmap(df_, annot=True)

def mds_embedding(df_flat):
    
    distance_matrix = flatdf_to_distance(df_flat)
    
    mds = MDS(n_components=2, metric=True, dissimilarity='precomputed', n_init=200)

    mds_emb = mds.fit_transform(distance_matrix)
    
    # lets plot it a little nicer
    df = pd.DataFrame(mds_emb, columns=['mds1', 'mds2'])
    df['samplename'] = distance_matrix.index
#     df['samplename'] = [_.split('_')[0] for _ in sample_names]

    plt.figure(figsize=(10,5))
    plt.subplot(121)
    sns.scatterplot(x='mds1', y='mds2', hue='samplename', data=df, legend=False)
    
    plt.subplot(122)
    sns.scatterplot(x='mds1', y='mds2', hue='samplename', data=df)
    plt.legend()
