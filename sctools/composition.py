import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import gmean
from scipy.spatial.distance import pdist, squareform
from scipy.stats import dirichlet
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
import tqdm

def log_geometric_mean(x, axis):
    """
    logarithm of the geometric mean along the given axis
    """
    D = x.shape[axis]
    loggmean = 1/D * np.sum(np.log(x), axis=axis, keepdims=True)
    return loggmean

def clr_transform(x, axis):
    """
    center log ratio transform for x, the features to be considered (cell types...) are along axis
    """
    return np.log(x) - log_geometric_mean(x, axis)

def aichinson_distance(x, axis=1):
    """
    Aitinson distance matrix of the compositional data x (in raw counts)
    - does CLR transform and computes the distances in CLR-space
    """
    y = clr_transform(x, axis)

    d = pdist(y)
    return squareform(d)

# clr is scale invariant
# x = np.array([[1,2,0,4], [10,11,12,13]])
# x_norm = x/x.sum(1, keepdims=1)
# np.testing.assert_allclose(clr_transform(x, axis=1),  clr_transform(x_norm, axis=1))
#
def clr_transform_bayesian(x, axis, n_samples):
    """
    Center log ratio transform, taking into account the multinomial uncertainty in the counts
    It'll return samples from the CLR posterior in the 3rd dimension
    """
    # NOTE: this is WRONG, dont normalize before, we need the full counts for the dirichlet!!
    # WRONG:     x_norm = x/x.sum(1, keepdims=1)

    x_bayes = []
    for i in tqdm.trange(len(x)):
        x_bayes.append(dirichlet(0.5 + x[i]).rvs(n_samples).T)
#         x_bayes.append(x[i])

    ## samples x features x mcmc_samples
    x_bayes = np.stack(x_bayes)

    x_clr = clr_transform(x_bayes,axis=axis)
    return x_clr

def aichinson_distance_bayesian(x, axis=1):
    """
    Aitchinson distance with Bayesian uncertainty from the counting
    3rd dimension will be samples from the posterior of distance d[i,j]
    """
    # estimate the clr transform (samples from the posterior)
    x_clr_bayes = clr_transform_bayesian(x, axis)

    # obs x features x mc_samples
    # we need to calc an obs x obs distance matrix
    # to get d(x_i,x_j) calcualte the distance between the respective mc_samples

    n_obs, n_features, n_mc = x_clr_bayes.shape
    d = np.zeros((n_obs, n_obs, n_mc))
    for i in range(n_obs):
        for j in range(n_obs):
            xi = x_clr_bayes[i,:,:].T
            xj = x_clr_bayes[j,:,:].T
            _d_tmp = np.sqrt(np.sum((xi-xj)**2, axis=1))
            assert _d_tmp.shape == (n_mc,)
            d[i,j,:] = _d_tmp
    return d


def compositional_pca(adata_coda):
    """
    PCA on the compositional adata, based on the posterior mean of the data in CLR space
    """
    x_bayes_clr = clr_transform_bayesian(adata_coda.X, axis=1)
#     D_bayes = aichinson_distance_bayesian(adata_coda.X)
    clr_posterior = x_bayes_clr.mean(axis=2)

    scaler = StandardScaler()
    clr_posterior = scaler.fit_transform(clr_posterior)
    pca = PCA()
    x_pca_posterior = pca.fit_transform(clr_posterior)

    principal_components_df = pd.DataFrame(pca.components_, columns=adata_coda.var.index )
    return pca, scaler, x_pca_posterior, principal_components_df


def plot_pca_with_uncertainty(pca, scaler, x_clr_bayes, components, color_vector=None):
    # project the uncertainty
    pca_dim1, pca_dim2 = components

    for i in range(x_clr_bayes.shape[0]): # iteratign over samples
        _x = x_clr_bayes[i].T
        _x = scaler.transform(_x)
        _x_pca = pca.transform(_x)
        plt.scatter(_x_pca[:,pca_dim1], _x_pca[:,pca_dim2], s=1, alpha=0.5)

    # if color_vector is None:
    #     plt.scatter(x_pca_posterior[:,pca_dim1], x_pca_posterior[:,pca_dim2])
    # else:
    #     plt.scatter(x_pca_posterior[:,pca_dim1], x_pca_posterior[:,pca_dim2], c=color_vector)
#         [color_dict_diagnosis[_] for _ in data_scanpy_1.obs.diagnosis]

def plot_pca_loadings(principal_components_df, components=(0,1)):

    # each row vector pca.components[i,:] is a single PC
    c1, c2 = components
    pc1 = principal_components_df.iloc[c1]
    pc2 = principal_components_df.iloc[c2]

    plt.scatter(pc1, pc2)
    for i in range(len(pc1)):
        plt.arrow(0,0,pc1.values[i], pc2.values[i], color='k')
        plt.text(x=pc1.values[i], y=pc2.values[i], s=pc1.index[i])


# actual hiereachical clustering on the posterior-mean distances
def clustered_heatmap_from_sccoda_CLR(sccoda_adata, figsize=(15, 5), barcolormap=None):
#     X = sccoda_adata[:, cluster_order].X.copy()
#     X = X/X.sum(1, keepdims=True)  # divide by total number of cells per sample

    from sccoda.util import data_visualization as viz
    from crukiopy.colormaps import color_dict_diagnosis
    diag_colors = [color_dict_diagnosis[_] for _ in sccoda_adata.obs.diagnosis]
    procedure_colors = ['red' if _ =="biopsy" else "blue" for _ in sccoda_adata.obs.procedure]

    X = sccoda_adata.X.copy()
    D_bayes = aichinson_distance_bayesian(sccoda_adata.X)
    D_posterior_mean = D_bayes.mean(axis=2)

    df_cluster = pd.DataFrame(X/X.sum(axis=1, keepdims=True), columns=sccoda_adata.var.index, index=sccoda_adata.obs.index)

    Z = linkage(squareform(D_posterior_mean), method='ward', optimal_ordering=True)
    sns.clustermap(df_cluster.T, row_cluster=False, col_linkage=Z,
                   col_colors=[diag_colors, procedure_colors],
                   xticklabels=True, yticklabels=True, figsize=figsize,
                   method='ward', #dendrogram_ratio=0.2,
#                    cmap=sns.dark_palette("#69d", reverse=False, as_cmap=True),
                   cmap=sns.color_palette("Blues"),
)

    ix = leaves_list(Z)
    order = [sccoda_adata.obs.index[i] for i in ix]
    fig = viz.stacked_barplot(sccoda_adata, feature_name="samples", figsize=(10,5),
                              level_order=order, #cmap = godsnot_cmap
                              cmap=barcolormap
                              )
    plt.xticks(rotation=90);
    return df_cluster, Z


"""
Conpoistional algebra/geometry
"""
def C(x, axis=1):
    """
    rescale onto simplex
    """
    return x/x.sum(axis, keepdims=True)

def perturbation(x, p):
    assert np.all(x.shape == p.shape)
    return C(x*p)

def power(x, alpha):
    assert np.isscalar(alpha)
    return C(x**alpha)

def inner_product(x,y, axis=1):
    if False:
        gx = gmean(x, axis=axis).reshape(-1,1)
        gy = gmean(y, axis=axis).reshape(-1,1)
        I = np.sum(np.log(x/gx) * np.log(y/gy), axis=axis)
        return I

    else:
        log_gx = log_geometric_mean(x, axis=axis)
        log_gy = log_geometric_mean(y, axis=axis)

    #     np.testing.assert_allclose(log_gx,  np.log(gx))
        _t = (np.log(x) - log_gx) * (np.log(y)-log_gy)

    #     np.testing.assert_allclose(I, np.sum(_t, axis=axis))

        return np.sum(_t, axis=axis)

def get_straight_line(x0,p, t_space):
    """
    a "straight" line in simplex space
    """
    x_traj = []
    for t in t_space:
        xnew = perturbation(x0, power(p, t)) 
        x_traj.append(xnew)

    x_traj = np.vstack(x_traj)
    return x_traj


def aitchinson_distance_pairwise(x,y):
    np.testing.assert_allclose(x , C(x))
    np.testing.assert_allclose(y, C(y))

#     a = np.log(x / gmean(x, axis=1))
#     b = np.log(y / gmean(y, axis=1))
    
    a = np.log(x) - log_geometric_mean(x, axis=1)
    b = np.log(y) - log_geometric_mean(y, axis=1)
    
    d = np.sqrt(np.sum((a-b)**2))
    
    return d


"""
Isometric logratio transform
"""

def isometric_basis_vectors_RD(i, D):
    "basis for compositional data in real space"
    assert i>0, "index starts at 1"
    assert i<=D-1, "index max at D-1"
    pre = np.sqrt(i/(i+1))

    u = [1/i]*i + [-1] + [0]*(D-i-1)
    return pre*np.array(u)

def isometric_basis_vectors_SD(i, D):
    "basis for compositional data in simplex space"
    assert i>0, "index starts at 1"
    assert i<=D-1, "index max at D-1"

    A = np.sqrt(1/(i*(i+1)))
    e = [A]*i + [-np.sqrt(i/(i+1))] + [0]*(D-i-1)

    return C(np.exp(e), axis=0)

def irl_transform(x, mode='SD'):
    assert mode in ["SD", "RD"]
    N_obs, D = x.shape

    if mode=="RD":
        # moving to R^D projection on the RD basis vectors
        x_clr = clr_transform(x, axis=1)
        irl = []
        for i in range(1,D):
            u = isometric_basis_vectors_RD(i, D)
            irl.append(x_clr @ u)
        return np.vstack(irl).T

    elif mode=="SD":
        # staying in S, projecting onto the SD basis
        # be sure to use the aitchinson inner product for projection!
        irl = []
        for i in range(1,D):
            u = isometric_basis_vectors_SD(i, D).reshape(1,-1)
            _t = inner_product(x,u).T
            
            irl.append(_t)
        return np.vstack(irl).T
    else:
        raise ValueError(f'unknown mode {mode}')