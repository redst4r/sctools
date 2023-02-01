import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import numpy as np
import sys
# sys.path.append('/home/michi/ms_python_packages/anndata2ri')
sys.path.append('/home/mstrasse/anndata2ri')
import anndata2ri

import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from patsy import dmatrices
import pandas as pd
import scanpy as sc

import itertools
def DESeq2_pseudobulk_wrapper(adata_pseudo, formula: str, var_of_interest: str):
    """
    :param formula: something like "~treatment + batch", needs to start with a ~
    :param var_of_interest: which variable in the formula is actually of interest (the other are just nuissance, regressed out 
        Shrinkage is performed for all comparisons including this var_of_interest, and returned

    :return: a dict of dataframes, corresponding to the different comparisons involing var_of_interest
    """
    assert formula.startswith('~'), "Formula must start with ~"
    assert isinstance(var_of_interest, str), "var_of_interest must be a string"
    importr("DESeq2")
    importr("PCAtools")
    """
    if missing do 
    ------------------------------------------
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("PCAtools")
    ------------------------------------------
    """

    # move into R
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv['scanpy.data'] = adata_pseudo

    # creating the DESeq object
    ro.r('coldata = colData(scanpy.data)')
    ro.r('counts = assay(scanpy.data)')
    ro.r(f'dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design= {formula})')

    # striagten out the levels
    # note that RPy2 respects pd.Categorical levels, hence we can do taht in python already!
    # ro.r('dds$diagnosis <- factor(dds$diagnosis, levels = c("NE", "M", "D", "T", "NS"))')

    # doing DEseq computations
    print('Main DESeq conputation')
    ro.r('dds <- DESeq(dds)')

    # some visualization
    ro.r('vsd <- vst(dds, blind=FALSE)')
    ro.r('p <- PCAtools::pca(assay(vsd), metadata = colData(dds), removeVar = 0.1, scale=F)')

    pandas2ri.activate()
    result_names = list(ro.r('resultsNames(dds)'))
    results_of_interest = [_ for _ in result_names if _.startswith(var_of_interest)]
    result_dict = {}
    for r in results_of_interest:
        print(f'shrinkage for {r}')
        ro.r(f'resLFC <- lfcShrink(dds, coef="{r}", type="apeglm")')
        df = ro.r('as.data.frame(resLFC)')
        result_dict[r] = df

    # get the PCA 
    df_pca = ro.r('p')
    dict_pca = dict(df_pca.items())
    df_vsd = ro.r('assay(vsd)')

    dict_pca['metadata'] = ro.r('as.data.frame')(dict_pca['metadata'])

    pandas2ri.deactivate()

    return result_dict, dict_pca, df_vsd


"""
code to correlate Covariates to PCA components
"""


def linear_model_pc_vs_covariate(pca_score, covariate):
    """
    regressing a covariate against PC scores
    i.e. does the PC score correlate with this covariate

    if categorical, the model will encode it!
    """
    # pca_scores = DF['PC1']
    # cont_covariate = DF['n_molecules']

    _df = pd.DataFrame({'pc': pca_score, 'covariate':covariate})
    y, X = dmatrices('pc~covariate', data=_df, return_type='dataframe')

    mod = sm.OLS(y, X)
    res = mod.fit()

    q = res.pvalues.copy()
    q.index = q.index.map(lambda x: "pval_"+x)

    w = res.params.copy()
    w.index = w.index.map(lambda x: "coeff_"+x)
    _series = pd.concat([q, w])

#     plt.figure()
#     plt.scatter(covariate, pca_score)
#     plt.scatter(covariate, res.fittedvalues)

    return _series, res

def is_categorical(array_like):
    return array_like.dtype.name == 'category'

def covariate_pc_correlation(df_pca, df_meta):
    results = []
    for pcname, covname in itertools.product(df_pca.columns, df_meta.columns):
        score_vector = df_pca[pcname].values
        cov_vector = df_meta[covname].values

        _series, res = linear_model_pc_vs_covariate(score_vector, cov_vector)

        # if its a categorical, we'll get several pvals. lets just take the minumim
        pval = np.min(
            _series[[_ for _ in _series.index if _.startswith('pval_covariate')]]
        )

        results.append({
            'PC': pcname,
            'covariate': covname,
            'pval': pval,
            'Rsquared': res.rsquared,
            'Rsquared_adj': res.rsquared_adj,
            'categorical_covariate': is_categorical(cov_vector)
        })

    results = pd.DataFrame(results)

    # make sure the PCs are sorted PC1,PC2,...PC10, PC11
    results['PC'] = pd.Categorical(results['PC'], categories=df_pca.columns)
    return results

import itertools

def plot_pca_grid(DF, pc_list, meta_list):
    df = []
    for pci, pcj in itertools.combinations(pc_list, 2):

        dp = pd.DataFrame({'score_i': DF[pci], 'score_j': DF[pcj], 'i':pci, 'j':pcj})
        dm = DF[meta_list]


        df.append(pd.concat([dp, dm], axis=1))

    return pd.concat(df)

