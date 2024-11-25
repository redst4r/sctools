import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import numpy as np
import anndata2ri

import statsmodels.api as sm
from patsy import dmatrices
import pandas as pd
import scanpy as sc

import itertools
def DESeq2_pseudobulk_wrapper(adata_pseudo, formula: str, vars_of_interest):
    """
    :param formula: something like "~treatment + batch", needs to start with a ~
    :param var_of_interest: which variable in the formula is actually of interest (the other are just nuissance, regressed out
        Shrinkage is performed for all comparisons including this var_of_interest, and returned

    :return: a dict of dataframes, corresponding to the different comparisons involing var_of_interest
    """
    assert formula.startswith('~'), "Formula must start with ~"

    # turn into a list of only single string element is provided
    if isinstance(vars_of_interest, str):
        vars_of_interest = [vars_of_interest]

    assert isinstance(vars_of_interest, list), "var_of_interest must be a string"
    importr("DESeq2")

    # move into R
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv['scanpy.data'] = adata_pseudo

    # creating the DESeq object
    ro.r('coldata = colData(scanpy.data)')
    ro.r('counts = assay(scanpy.data)')
    ro.r(f'dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design= {formula})')

    # doing DEseq computations
    print('Main DESeq conputation')
    ro.r('dds <- DESeq(dds)')

    pandas2ri.activate()
    result_dict = {}
    result_names = list(ro.r('resultsNames(dds)'))
    for var in vars_of_interest:
        results_of_interest = [_ for _ in result_names if _.startswith(var)]
        for r in results_of_interest:

            # to get unshrunk fold changes
            ro.r(f'res <- results(dds, name="{r}")')
            df_unshrunk = ro.r('as.data.frame(res)')

            print(f'shrinkage for {r}')
            ro.r(f'resLFC <- lfcShrink(dds, coef="{r}", type="apeglm")')
            df = ro.r('as.data.frame(resLFC)')

            # merge the shrunk and unshrunk
            df_merged = df.merge(
                df_unshrunk[['log2FoldChange', 'lfcSE']].rename({'log2FoldChange': 'log2FoldChange_unshrunk', 'lfcSE': 'lfcSE_unshrunk'}, axis=1),
                left_index=True, right_index=True
            )
            result_dict[r] = df_merged

    ro.r('vsd <- vst(dds, blind=FALSE)')
    df_vsd = ro.r('assay(vsd)')
    adata_vsd = sc.AnnData(df_vsd.T, obs=adata_pseudo.obs, var=adata_pseudo.var)
    pandas2ri.deactivate()

    if False:
        # some visualization
        importr("PCAtools")
        ro.r('p <- PCAtools::pca(assay(vsd), metadata = colData(dds), removeVar = 0.1, scale=F)')
        # get the PCA
        # df_pca = ro.r('p')
        # dict_pca = dict(df_pca.items())
        pandas2ri.activate()

        dict_pca = {}
        dict_pca['metadata'] = ro.r('as.data.frame(p$metadata)')
        dict_pca['rotated'] = ro.r('as.data.frame(p$rotated)')
        dict_pca['loadings'] = ro.r('as.data.frame(p$loadings)')
        dict_pca['dispersions'] = ro.r('as.data.frame(mcols(dds))')
        dict_pca['size_factors'] = ro.r('sizeFactors(dds)')
        # dict_pca['variance'] = ro.r('as.data.frame(p$variance)')
        # dict_pca['sdev'] = ro.r('as.data.frame')(dict_pca['sdev'])
        # dict_pca['xvars'] = ro.r('as.data.frame')(dict_pca['xvars'])
        # dict_pca['yvars'] = ro.r('as.data.frame')(dict_pca['yvars'])
        # dict_pca['components'] = ro.r('as.data.frame')(dict_pca['components'])
        dict_pca['df_pca'] = dict_pca['rotated'].merge(dict_pca['metadata'], left_index=True, right_index=True)

        pandas2ri.deactivate()
        return result_dict, dict_pca, adata_vsd, ro.r('dds')
    else:
        return result_dict, adata_vsd, ro.r('dds')


def DESeq2_pseudobulk_wrapper_LRT(adata_pseudo, formula_full: str, formula_reduced: str):
    """
    Differential expression as determined by a LRT: wchich genes are significantly better explained by the full model than by the baseline
    :param formula_full: something like "~treatment + batch", needs to start with a ~
    :param formula_reduced: formula of the baseline model in the LRT
    :return: a dict of dataframes, corresponding to the different comparisons involing var_of_interest
    """
    assert formula_full.startswith('~'), "Formula must start with ~"
    assert formula_reduced.startswith('~'), "Formula must start with ~"
    importr("DESeq2")

    # move into R
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv['scanpy.data'] = adata_pseudo

    # creating the DESeq object
    ro.r('coldata = colData(scanpy.data)')
    ro.r('counts = assay(scanpy.data)')
    ro.r(f'dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design= {formula_full})')

    # striagten out the levels
    # note that RPy2 respects pd.Categorical levels, hence we can do taht in python already!
    # ro.r('dds$diagnosis <- factor(dds$diagnosis, levels = c("NE", "M", "D", "T", "NS"))')

    # doing DEseq computations
    print('Main DESeq conputation')
#     ro.r('dds <- DESeq(dds)')
    ro.r(f'dds <- DESeq(dds, test="LRT", reduced = {formula_reduced})')
    ro.r('vsd <- vst(dds, blind=FALSE)')

    pandas2ri.activate()
    df_DE = ro.r('as.data.frame(results(dds))')
    df_vsd = ro.r('assay(vsd)')
    pandas2ri.deactivate()

    adata_vsd = sc.AnnData(df_vsd.T, obs=adata_pseudo.obs, var=adata_pseudo.var)

    
    if False:
        """ currently, PCAtools has a weird error: 
        [matrixStats (>= 1.2.0)] useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE
        """
        importr("PCAtools")

        # some visualization
        ro.r('p <- PCAtools::pca(assay(vsd), metadata = colData(dds), removeVar = 0.1, scale=F)')

        pandas2ri.activate()
        dict_pca = {}
        dict_pca['metadata'] = ro.r('as.data.frame(p$metadata)')
        dict_pca['rotated'] = ro.r('as.data.frame(p$rotated)')
        dict_pca['loadings'] = ro.r('as.data.frame(p$loadings)')

        # get the PCA
        # df_pca = ro.r('p')
        # dict_pca = dict(df_pca.items())
        # dict_pca['metadata'] = ro.r('as.data.frame')(dict_pca['metadata'])

        pandas2ri.deactivate()

        dict_pca['df_pca'] = dict_pca['rotated'].merge(dict_pca['metadata'], left_index=True, right_index=True)
        return df_DE, dict_pca, adata_vsd, ro.r('dds')
    
    else:
        return df_DE, adata_vsd, ro.r('dds')


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


def plot_pca_grid(DF, pc_list, meta_list):
    df = []
    for pci, pcj in itertools.combinations(pc_list, 2):

        dp = pd.DataFrame({'score_i': DF[pci], 'score_j': DF[pcj], 'i':pci, 'j':pcj})
        dm = DF[meta_list]


        df.append(pd.concat([dp, dm], axis=1))

    return pd.concat(df)

