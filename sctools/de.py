import itertools
import scanpy as sc
import numpy as np
from scipy import stats
from scipy.sparse import spmatrix, csc_matrix
import pandas as pd
import tqdm
import plotnine as pn
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from patsy import dmatrices
"""
tools for differential expression in scanpy
"""


def get_de_genes(adata, cluster, as_list=False, key='rank_genes_groups'):
    df = scanpy_DE_to_dataframe_fast(adata, key)[cluster]

    genes = df['name']
    genes = genes.replace({np.nan: 'nan'})

    return genes if as_list else " ".join(genes)


def DEG_one_vs_all(Q):

    """
    t-test of each leiden cluster against all other clusters
    """

    # precompute the groups. as a matter of fact, compute all the stats for the ttest
    groupnames = Q.obs['leiden'].unique()

    print('precomputing')
    # prec_groups= {gname: Q[Q.obs.query(f'leiden==@gname').index,:] for gname in groupnames}

    prec_groups_means= {gname: Q[Q.obs.query(f'leiden==@gname').index,:].raw.X.A.mean(0) for gname in groupnames}
    prec_groups_vars = {gname: Q[Q.obs.query(f'leiden==@gname').index,:].raw.X.A.var(0) for gname in groupnames}
    prec_groups_nobs = {gname: len(Q[Q.obs.query(f'leiden==@gname').index,:]) for gname in groupnames}

    the_df = []

    # we deliberately iterate over the product instead of combinations (n^2 vs n*(n-1)/2)
    # such that each cluster is represented as group1
    # otherwise it gets very confusing to look for upregultated genes  of a specific cluster
    for i, j in itertools.product(groupnames, groupnames):
        if i == j: continue

        g1m = prec_groups_means[i]
        g2m = prec_groups_means[j]
        g1s = np.sqrt(prec_groups_vars[i])
        g2s = np.sqrt(prec_groups_vars[j])

        with np.errstate(invalid="ignore"):
            scores, pvals = stats.ttest_ind_from_stats(
                mean1=g1m, std1=g1s, nobs1=prec_groups_nobs[i],
                mean2=g2m, std2=g2s, nobs2=prec_groups_nobs[j],
                equal_var=False  # Welch's
            )
        ddd = [{'score': s,
                'pval': p,
                'name':n,
                'group1':i,
                'group2':j,
                'group1_mean': _g1m,
                'group2_mean': _g2m,
                'group1_std': _g1s,
                'group2_std': _g2s,
               }                                                      # WARNING: the RAW is very important here!!
               for s,p,n, _g1m, _g2m, _g1s, _g2s in zip(scores, pvals, Q.raw.var.index.values, g1m, g2m, g1s,g2s) if not np.isnan(s) and p<0.1 and s>0] # s>0 such that we only look at upregulated genes
        the_df.extend(ddd)

    return pd.DataFrame(the_df)


def DEG_one_vs_all_aggregate(de):
    """
    aggregate the genes: look for ones that are significant against ever other cluster
    i.e. three clusters A,B,C: A-B, A-C is significant,
    then we call the gene a marker of cluster A
    """
    df_all_siginificant = []
    groups = de.group1.unique()
    for cl in tqdm.tqdm(groups):
        tmp1 = de.query('group1==@cl')  # all tests involving cluster cl
        # see wichi genes had significant tests across all clusters
        # problem is: we already filtered some genes out in `DEG_one_vs_all`, hence not all groups will be present
        # so, we only look at genes that have a test vs all remaining groups
        all_significant = tmp1.groupby('name').apply(lambda g: np.all(g.pval < 0.001) and len(g) == len(groups)-1)
        all_significant = all_significant[all_significant==True].index.values

        max_pval = tmp1.groupby('name').apply(lambda g: np.max(g.pval))

        df_all_siginificant.extend([{'group': cl, 'name': _, 'pval': max_pval[_]} for _ in all_significant])
    df_all_siginificant = pd.DataFrame(df_all_siginificant)
    return df_all_siginificant


def scanpy_DE_to_dataframe_fast(adata, key='rank_genes_groups'):
    """
    to just get DE genes and their scores/pvals out of scanpy

    faster version to `scanpy_DE_to_dataframe`

    ------------------------
    A note on the differential expression in scanpy, esp fold changes:
    for a single gene X and groups A and B the fold change is calculated as:

    # log(1+x)
    yA = np.log1p(X[group==A])
    yB = np.log1p(X[group==B])

    # averaging
    mA = yA.mean()
    mB = yB.mean()

    # undo log1p
    nA = np.exp(nA)-1
    nB = np.exp(nB)-1

    # calc log2 foldchange
    np.log2(nA/nB)
    """

    rank_dict = adata.uns[key]
    df = []

    groupby = adata.uns[key]['params']['groupby']
    groupnames = adata.uns[key]['names'].dtype.names

    for i in range(len(rank_dict['scores'])):
        # the items always come in pairs: up in group 1, up in group2
        # each of the following is of size: 1x(groups)
        """
        the storage format is kind of annoying:
        Its n-tuples (n=number of groups) and a differen gene for each group
        ('MTATP6P1', 'EEF1A1', 'HSP90AA1', 'HIST1H1E'),
        ('MPO', 'TPT1', 'TUBA1B', 'HIST1H1C')
        """
        s = rank_dict['scores'][i]
        n = rank_dict['names'][i]
        p = rank_dict['pvals'][i]
        q = rank_dict['pvals_adj'][i]
        fc = rank_dict['logfoldchanges'][i]

        for j in range(len(s)):  # iterationg over all groups, adding DE of group vs rest
            gname = groupnames[j]
            genename = n[j]
            tmp = { 'score': s[j],
                    'name': genename,
                    'pval': p[j],
                    'qval': q[j],
                    'logfoldchange': fc[j],
                    # 'symbol': adata.var.loc[genename].name, # ['symbol'],
                    'symbol': genename,  # ['symbol'],
                    'group': j,
                    'groupname': gname
                    }
            df.append(tmp)

    df = pd.DataFrame(df)
    # the df contains a mix of De for all grous
    # create on DF per group
    group_dfs = {}
    for g in df['groupname'].unique():
        tmp = df.query('groupname==@g').copy()
        group_dfs[g] = tmp

    return group_dfs


def split_adata_raw(adata, groupby):
    """
    splits the given adata into groups, expression matrix will be from the .raw.X,
    also converts X tp csc-format (fast access to single columns)
    """

    import scipy.sparse
    # percompute the cluster datasets
    # also change the raw matrix into sparse column format for fast indxein
    prec_datasets = {}  # gname: adata[adata.obs.query(f'{groupby}==@gname').index,:] for gname in groupnames}
    groupnames = adata.obs[groupby].unique()
    for gname in groupnames:
        _tmp_adata = adata[adata.obs.query(f'{groupby}==@gname').index,:]

        rawX = scipy.sparse.csc_matrix(_tmp_adata.raw.X)
        prec_datasets[gname] = sc.AnnData(rawX,
                                          obs=_tmp_adata.obs,
                                          var=_tmp_adata.raw.var)
    return prec_datasets


def scanpy_DE_to_dataframe(adata, key='rank_genes_groups'):
    """
    turn the differential expression from scanpy (rank...)
    into a more readable dataframe
    """
    rank_dict = adata.uns[key]
    df = []
    groupby = adata.uns[key]['params']['groupby']
    groupnames = adata.uns[key]['names'].dtype.names

    prec_datasets = split_adata_raw(adata, groupby)

    gene_to_index_raw = {g:i for i,g in enumerate(adata.raw.var.index)}

    for i in range(len(rank_dict['scores'])):
        # the items always come in pairs: up in group 1, up in group2

        # each of the following is of size: 1x(groups)
        s = rank_dict['scores'][i]
        n= rank_dict['names'][i]
        p= rank_dict['pvals'][i]
        q= rank_dict['pvals_adj'][i]
        fc = rank_dict['logfoldchanges'][i]

        for j in range(len(s)): # iterationg over all groups, adding DE of group vs rest

            # additionally count how many cells express that marker at all
            # in the cluster of interest
            # (avoids genes that are upregulated strognly in a single cell but otherwise not expressed at all)
            gname = groupnames[j]
            genename = n[j]

            if not isinstance(genename, str) and np.isnan(genename):
                continue # the min_percent expressing cells sets gene names to nan if they done fulfill the criterion

            adata_cluster = prec_datasets[gname]
            # X = adata_cluster[:, genename].X
            X = adata_cluster.X[:, gene_to_index_raw[genename]]
            percent_expressing = np.sum(X>0)/X.shape[0]
            raw_avg_expression = np.mean(X)
            tmp  = {'score': s[j],
                    'name': genename,
                    'pval': p[j],
                    'qval': q[j],
                    'logfoldchange': fc[j],
                    # 'symbol': adata.var.loc[genename].name, # ['symbol'],
                    'symbol': genename, # ['symbol'],
                    'group': j,
                    'groupname': gname,
                    'percent_expressing': percent_expressing,
                    'raw_avg_expression': raw_avg_expression}

            df.append(tmp)

    df = pd.DataFrame(df)
    # the df contains a mix of De for all grous
    # create on DF per group
    group_dfs = {}
    for g in df['groupname'].unique():
        tmp = df.query('groupname==@g').copy()
        group_dfs[g] = tmp

    return group_dfs


def gene_expression_to_flat_df_NEW(adata, genes, additional_vars, scale, use_raw):
    """
    scale: put all genes onto the same scale (mean0, std1)
    """
    flat_df = []

    for the_gene in genes:
        if use_raw:
            X = adata.raw[:, the_gene].X
        else:
            X = adata[:, the_gene].X
        if isinstance(X, spmatrix):
            X = X.A
        X = X.flatten()
        # standard scaling
        if scale:
            # this really scales the genes across groups,
            # as we do the group selection on X only below!
            X = X - X.mean()
            X = X / X.std()

        _tmp_df = pd.DataFrame()
        _tmp_df['expression'] = X
        _tmp_df['gene'] = the_gene
        for v in additional_vars:
            _tmp_df[v] = adata.obs[v].values

        flat_df.append(_tmp_df)

    flat_df = pd.concat(flat_df)
    return flat_df


def gene_expression_to_flat_df(adata, genes, grouping_var, scale, use_raw):
    """
    scale: put all genes onto the same scale (mean0, std1)
    """
    groups = adata.obs[grouping_var].unique()
    flat_df = []

    for the_gene in genes:
        for the_group in groups:
            ix_group = adata.obs[grouping_var] == the_group  # which condition the datapoint belongs to
            if use_raw:
                X = adata.raw[:, the_gene].X
            else:
                X = adata[:, the_gene].X
            if isinstance(X, spmatrix):
                X = X.A
            X = X.flatten()
            # standard scaling
            if scale:
                # this really scales the genes across groups,
                # as we do the group selection on X only below!
                X = X - X.mean()
                X = X / X.std()

            X = X[ix_group]

            _tmp_df = pd.DataFrame()
            _tmp_df['expression'] = X
            _tmp_df['gene'] = the_gene
            _tmp_df[grouping_var] = the_group

            flat_df.append(_tmp_df)

    flat_df = pd.concat(flat_df)
    return flat_df


# TODO this is highly redundant with gene_expression_to_flat_df
def de_adata_to_flat_df(adata, scale:bool, ngenes, qval_cutoff, key='rank_genes_groups'):
    _tmp = scanpy_DE_to_dataframe_fast(adata)

    de_genes = {g: df.query(f'qval<{qval_cutoff}').name.head(ngenes).tolist() for g,df in _tmp.items()}

    group = adata.uns[key]['params']['groupby']

    # building a long dataframe with
    # cell_index, genename, expression, which_de_group
    flat_df = []

    for c, genes in de_genes.items():
        for g in genes:  # de-genes for that group
            X = adata.raw[:, g].X  # thats the expression of a single gene all groups

            # annoying, sometimes is sparse sometimes its not!
            if isinstance(X, spmatrix):
                X = X.A
            # standard scaling
            if scale:
                X = X - X.mean()
                X = X / X.std()

            _tmp_df = adata.obs[[group]].copy()  # which condition the datapoint belongs to
            _tmp_df['expression'] = X
            _tmp_df['gene'] = g
            _tmp_df['which_de_group'] = c  # which group is this gene DE

            flat_df.append(_tmp_df)

    flat_df = pd.concat(flat_df)
    return flat_df, group


def plot_de(adata, scale=True, mode='boxplot', ngenes=5, qval_cutoff=0.2):

    """
    some ggplot-based display of the differential expression results.

    for each group of differnetially expressed genes (multiple genes per condition)
    do boxplots of how the cells in the different conditions compare

    big bonus: cells of different conditions are next to each other, easoer to compare than
    `sc.pl.rank_genes_groups_violin`

    mode: either `violin` or 'boxplot'
    """
    gg_df, group = de_adata_to_flat_df(adata, scale, ngenes, qval_cutoff)

    geom_dict = {
        'boxplot':pn.geom_boxplot,
        'violin': pn.geom_violin
        }

    plot = \
    pn.ggplot(gg_df, pn.aes(x='factor(gene)', y='expression', fill=f'factor({group})'))\
      + geom_dict[mode]() \
      + pn.facet_wrap('~which_de_group', scales='free') \
      + pn.theme(axis_text_x=pn.element_text(rotation=90, hjust=1))
      # + pn.scale_y_log10() # pn.geom_violin() # + pn.scale_y_log10()

    return plot


def _lm_de_helper(df, formula):
    """
    :param formula: something like: gene ~ diagnosis + patient
    """
    y, X = dmatrices(formula, data=df, return_type='dataframe')

    mod = sm.OLS(y, X)
    res = mod.fit()

    q = res.pvalues.copy()
    q.index = q.index.map(lambda x: "pval_"+x)

    w = res.params.copy()
    w.index = w.index.map(lambda x: "coeff_"+x)
    _series = pd.concat([q, w])
    return _series

def differential_expression_lm(adata, formula):
    """
    formula: left hand side: gene, log_gene, rank_gene, norm_gene
    """
    X = csc_matrix(adata.raw.X)  # usually much faster in csc
    results = []
    for i, gene in tqdm.tqdm(enumerate(adata.raw.var.index), total=adata.raw.shape[1]):
        df_tmp = adata.obs.copy()
        df_tmp['gene'] = X[:, i].A
        df_tmp['log_gene'] = np.log10(1+df_tmp['gene'])
        df_tmp['rank_gene'] = df_tmp['gene'].rank()
        df_tmp['norm_gene'] = df_tmp['gene'] / adata.obs['n_molecules']

        if df_tmp['gene'].sum() == 0:  # gene not expressed
            continue

        _series = _lm_de_helper(df_tmp, formula)
        _series['gene'] = gene
        _series['avg_expression'] = df_tmp['gene'].mean()
        _series['avg_normed_expression'] = df_tmp['norm_gene'].mean()

        results.append(_series)
    results = pd.DataFrame(results)

    # mutiple testing correaction
    pval_fields = [_ for _ in results.columns if _.startswith('pval_')]
    for pf in pval_fields:
        _, q, _, _ = multipletests(results[pf], method='fdr_bh')

        # rename to qval_ ...
        qname = 'q'+pf[1:]
        results[qname] = q

    return results



def differential_expression_lm_parallel(adata, formula, cores=4):
    """
    chunk the expression matrix in multiple column blocks, do DE on each block in parallel
    """


import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import sys
# sys.path.append('/home/michi/ms_python_packages/anndata2ri')
sys.path.append('/home/mstrasse/anndata2ri')
import anndata2ri

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

    # move into R
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv['scanpy.data'] = adata_pseudo

    # creating the DESeq object
    ro.r('coldata = colData(scanpy.data)')
    ro.r('counts = assay(scanpy.data)')
    ro.r('dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design= ~ patient + diagnosis)')

    # striagten out the levels
    # note that RPy2 respects pd.Categorical levels, hence we can do taht in python already!
    # ro.r('dds$diagnosis <- factor(dds$diagnosis, levels = c("NE", "M", "D", "T", "NS"))')

    # doing DEseq computations
    print('Main DESeq conputation')
    ro.r('dds <- DESeq(dds)')

    pandas2ri.activate()
    result_names = list(ro.r('resultsNames(dds)'))
    results_of_interest = [_ for _ in result_names if _.startswith(var_of_interest)]
    result_dict = {}
    for r in results_of_interest:
        print(f'shrinkage for {r}')
        ro.r(f'resLFC <- lfcShrink(dds, coef="{r}", type="apeglm")')
        df = ro.r('as.data.frame(resLFC)')
        result_dict[r] = df
    pandas2ri.deactivate()

    return result_dict
