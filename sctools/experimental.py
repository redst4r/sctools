from samalg import SAM
import pandas as pd
import scanpy as sc

"""
some experimental features/code not to be included in the main pipeline yet
"""

def SAM_get_gene_weights(adata, batch_field):
    """
    doing self-assembling manifold analysis (on each batch separately) to identify
    variable genes (called weights in their paper)

    returns a DataFrame, each sample is a row, each gene a column and the importance/weight of that gene as the value
    """
    df_weights = []
    for s in adata.obs[batch_field].unique():
        print(f'running SAM for {s}')
        b = adata[adata.obs[batch_field]==s]
        sadata = sc.AnnData(X=b.raw.X, obs=b.obs, var=b.raw.var)
        sam=SAM(counts=sadata)
        sam.preprocess_data() # log transforms and filters the data
        sam.run(projection=None) #skip umap projection
        weights = sam.adata.var['weights']
        weights.name = s
        df_weights.append(weights)

    df_weights =pd.DataFrame(df_weights)
    return df_weights

def SAM_recipe(adata, batch_field,  n_top_genes, log=True):
    """
    pretty much doing the same as `sc.pp.recipe_zheng17`, but using SAM to determine highly variable genes.
    SAM is run on each batch independently, and a genes weight is the maximum weight across all batches

    """
    assert adata.raw, "needs adata.raw"
    df_weights = SAM_get_gene_weights(adata, batch_field)
    aggr_weights = df_weights.max(0)  # rows are the samples!
    aggr_weights.name = 'SAM_weights'
    adata.var['SAM_weights'] = aggr_weights  # save the weights for later inspection
    HVG_list = aggr_weights.sort_values().tail(n_top_genes).index.values

    adata = adata[:, HVG_list]     # subset the genes
    sc.pp.normalize_per_cell(adata)                 # renormalize after filtering
    if log:
        sc.pp.log1p(adata)                      # log transform: adata.X = log(adata.X + 1)
    sc.pp.scale(adata)
    return adata

import numpy as np
import harmonypy as ha
def harmony_iterative(adata, vars_use, niter):
    """
    this applies the harmony correction iteratively and gives access to the
    corrected PCA projections (Z_corr_array) and the membership matrix (R)
    """
    data_mat = np.array(adata.obsm['X_pca'])

    # just initialize
    ho = ha.run_harmony(data_mat,  adata.obs, vars_use, max_iter_harmony=0)
    Z_corr_array = []
    R_array = []

    for i in range(niter):
        ho.harmonize(iter_harmony=1, verbose=True)
        Z_corr_array.append(ho.Z_corr.T)
        R_array.append(ho.R)

    return Z_corr_array, R_array
