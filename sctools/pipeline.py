import scanpy as sc
from sctools import annotate_qc_metrics, annotate_coding_genes
from sctools.annotations import annotate_cellcycle
import numpy as np
import scrublet as scr
import pandas as pd
import gc
import harmonypy as ha


class Verbose(object):
    """
    context manager for verbose output around statements
    """

    def __init__(self, msg, verbose):
        self.msg = msg
        self.verbose = verbose

    def __enter__(self):
        if self.verbose:
            print(self.msg)
        return

    def __exit__(self, type, value, traceback):
        if self.verbose:
            print(f'Done {self.msg}')


def standard_processing(adata, detect_doublets=True, MITO_CUTOFF=0.4):
    """
    Wrapper around `michi_kallisto_recipe()` with standardized values for
    QC cutoffs, does Harmony correction and doublet detection!
    """
    UMI_CUTOFF = 1000
    VARIABLE_GENES = 4000

    DE_min = 0.25  # minimum fraction of cells that must express the gene in the cluster
    DE_max = 1  # max % cells outside of the cluster that express the gene
    adata.raw = adata.copy()
    adata = michi_kallisto_recipe(adata,
                                  umi_cutoff=UMI_CUTOFF,
                                  n_top_genes=VARIABLE_GENES,
                                  percent_mito_cutoff=MITO_CUTOFF,
                                  annotate_cellcycle_flag=True,
                                  harmony_correction='samplename',
                                  verbose=True)

    if detect_doublets:
        adata = annotate_doublets(adata, groupby='samplename')

    differential_expression_michi_kallisto_recipe(adata,
                                                  groupby='leiden',
                                                  n_genes=100,
                                                  method='wilcoxon',
                                                  min_in_group_fraction=DE_min,
                                                  max_out_group_fraction=DE_max)
    return adata


def michi_kallisto_recipe(adata, umi_cutoff=1000, n_top_genes=4000, percent_mito_cutoff=1,
                          annotate_cellcycle_flag=True,
                          harmony_correction=None,
                          harmony_clusters=None,
                          verbose=True):
    """
    filters for coding genes, adds QC, filters cells based on UMI, applies
    Zheng recipe calulates pca/nn/leiden/paga/umap

    adata: should be an object derived from kallisto
    """

    # make sure we get an unprocessed adata!
    assert adata.raw, ".raw does not exist!"
    assert np.all(adata.X.data == np.round(adata.X.data)), ".X must be ints/counts (maybe log1p'd?)"
    assert np.all(adata.raw.X.data == np.round(adata.raw.X.data)), ".raw.X must be ints/counts (maybe log1p'd?)"

    adata.uns['log_X'] = False  # keep track of whats logarithmized or not
    adata.uns['log_raw.X'] = False
    adata = preprocessing_michi_kallisto_recipe(adata, umi_cutoff, percent_mito_cutoff, verbose=True)

    """
    annotating the cellcycle on the reduced data (filtered for potential cells, not just CB)
    since preprocessing doenst change anything (except kicking out genes/cells)
    its alot more feasable after filtering!
    Note; this might still take ALOT of memory if more then 10k cells remain after filtering
    not sure what the cell-cycle scoring function in scanpy is doing internally!
    """
    # if verbose: print('annotating cell cycle')
    if annotate_cellcycle_flag:
        if verbose:
            print('Annotating Cell cycle')
        adata = annotate_cellcycle(adata)
        if verbose:
            print('Done: Annotating Cell cycle')

    if verbose:
        print('Zheng recipe')
    sc.pp.recipe_zheng17(adata, n_top_genes=n_top_genes, log=True, plot=False, copy=False)

    if verbose:
        print('Done: Zheng recipe')
    # zheng17 log1p the .X data, but doesnt touch .raw.X!
    adata.uns['log_X'] = True

    adata = postprocessing_michi_kallisto_recipe(adata, harmony_correction, harmony_clusters=harmony_clusters,verbose=verbose)

    assert np.all(adata.raw.X.data == np.round(adata.raw.X.data)), ".raw.X got log'd!"

    return adata


def preprocessing_michi_kallisto_recipe(adata, umi_cutoff, percent_mito_cutoff, verbose=True):
    """
    filtering (coding/n_genes/n_umis/mito), annotations, but no transformations of the data
    """
    if verbose:
        print('annotating QC')
    adata = annotate_qc_metrics(adata)

    if verbose:
        print('annotating and filtering for coding genes')
    adata = annotate_coding_genes(adata)
    adata = adata[:, adata.var.is_coding == True].copy()  # copying to avoid getting a view, which has some issues when copying later
    gc.collect()
    if verbose:
        print('filtering cells for UMI content')
    cells_before = adata.shape[0]
    adata = adata[adata.obs.query('n_molecules>@umi_cutoff').index].copy()  # copying to avoid getting a view, which has some issues when copying later
    gc.collect()
    cells_after = adata.shape[0]
    if verbose:
        print(f'Cells: {cells_before} -> {cells_after}')

    if verbose:
        print('filtering cells for mito content')
    cells_before = adata.shape[0]
    adata = adata[adata.obs.query('percent_mito<@percent_mito_cutoff').index].copy()  # copying to avoid getting a view, which has some issues when copying later
    gc.collect()

    cells_after = adata.shape[0]
    if verbose:
        print(f'Cells: {cells_before} -> {cells_after}')

    return adata


def postprocessing_michi_kallisto_recipe(adata, harmony_correction, harmony_clusters=None, verbose=True):
    """
    doing batch correction, PCA/UMAP etc
    """
    if verbose:
        print('PCA')
    sc.pp.pca(adata)

    # storing the original/uncorrected PCA
    adata.obsm['X_pca_original'] = adata.obsm['X_pca'].copy()
    gc.collect()
    """
    harmony batch correction if desired
    this works on the PCA proejction

    Harmony args:
    data_mat = V, ## PCA embedding matrix of cells
    meta_data = meta_data, ## dataframe with cell labels
    theta = 1, ## cluster diversity enforcement
    vars_use = 'dataset', ## variable to integrate out
    nclust = 5, ## number of clusters in Harmony model: corresponds to the parameter K in the manuscript.
    max.iter.harmony = 0, ## stop after initialization
    return_object = TRUE, ## return the full Harmony model object
    do_pca = FALSE ## don't recompute PCs

    PHI: design matrix

    """
    if harmony_correction:
        # assert harmony_clusters, "must set harmony_clusters!"
        n_batches = len(adata.obs[harmony_correction].unique())
        # if there's only a single sample, no batch correction needed
        # actually, harmony crashes with some error if you run it on a single batch (harmonypy 0.0.4)
        if n_batches > 1:
            if verbose:
                print('Harmony batch correction')
            vars_use = [harmony_correction]  # samplenames for harmony
            assert harmony_correction in adata.obs.columns
            # get out the PCA matrix
            data_mat = np.array(adata.obsm['X_pca'])
            # and harmonize
            ho = ha.run_harmony(data_mat,  adata.obs, vars_use, max_iter_harmony=25, nclust=harmony_clusters)
            adata.obsm['X_pca'] = np.transpose(ho.Z_corr)

    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=1)
    sc.tl.louvain(adata)
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata, color=['leiden'], show=False)
    paga_init_pos = sc.tl._utils.get_init_pos_from_paga(adata)
    sc.tl.umap(adata, init_pos=paga_init_pos)

    # sc.tl.draw_graph(adata, init_pos='paga') # somehow the init_pos works different here

    return adata


def differential_expression_michi_kallisto_recipe(adata, groupby, n_genes=100, method='wilcoxon', min_in_group_fraction=0, max_out_group_fraction=1, use_raw=True):
    """
    scanpy by default runs the differential expression on .raw.X, but also assumes
    logarithmized data. Otherwise, pvals come out wrong!

    To leave .raw unchanged, we log-transform, do DE-analysis and undo the log_transform

    TODO this better be an atomic operation, otherwise the log-flag gets screwed up
    or stick it into a decorator!

    The filtering for fraction of cells expressing is a bit hacky atm. It relies on sc.tl.filter_rank_genes... which just sets genenames to nan instead of removing them


    :param adata: AnnData object
    :param groupby: .obs column denoting the grouping of differential expression
    :param n_genes: number of DE genes to report at max
    :param method: what statistical test to do, see sc.tl.rank_genes_groups for options
    :param min_in_group_fraction: \in [0,1]. Additionally filter the DE genes: At least x% of cells must express this gene in the upregulated cluster.
    Prevents genes that are extremly high in a few cells from dominating the DE list. Useful for marker genes (every cell in a cluster should express that marker!)
    :param max_out_group_fraction: similar, just force genes to be exclusively expressed in the cluster: Filter genes that have more then x% cells expressing it outside the cluster
    """
    assert not adata.raw is None, "no data is present in the .raw storage. Differential expression will is only done on the .raw data!"
    assert not adata.uns['log_raw.X']
    adata.raw.X.data = np.log1p(adata.raw.X.data)
    sc.tl.rank_genes_groups(adata, groupby=groupby, n_genes=n_genes, method=method, use_raw=use_raw)
    # undoing the log
    adata.raw.X.data = np.round(np.exp(adata.raw.X.data) - 1)

    # filtering
    sc.tl.filter_rank_genes_groups(adata,
                                   log=False,  # since we undid the log!
                                   key_added='rank_genes_groups_filtered',
                                   min_in_group_fraction=min_in_group_fraction,
                                   min_fold_change=0,  # not filteirng for fold_change
                                   max_out_group_fraction=max_out_group_fraction) # not filetering for %expressing outside the cluster
    # now make the filtered genes the default DE genes
    adata.uns['rank_genes_groups_unfiltered'] = adata.uns['rank_genes_groups']
    adata.uns['rank_genes_groups'] = adata.uns['rank_genes_groups_filtered']


def export_for_cellxgene(adata, annotations):
    # for Export we have to pull all the genes back into the adata.X!
    _tmp = sc.AnnData(adata.raw.X,
                      obs=adata.obs[annotations],
                      var=adata.raw.var.copy(),
                      uns=adata.uns,
                      obsm=adata.obsm,
                      layers=adata.layers,
                      obsp=adata.obsp,
                      # these two wont work: adata.raw is N x 35000, but adata.var and varm are N x 4000
                      # varm=adata.varm,
                      # varp=adata.varp
                      )
    return _tmp


def annotate_doublets(adata, groupby='samplename', PLOTTING=False):
    doublet_vectors = []
    for s in adata.obs[groupby].unique():
        print(s)
        b = adata[adata.obs[groupby] == s]

        if len(b) > 50:
            scrub = scr.Scrublet(b.raw.X, expected_doublet_rate=0.06, sim_doublet_ratio=5)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                                      min_cells=3,
                                                                      min_gene_variability_pctl=85,
                                                                      n_prin_comps=50)
            if PLOTTING:
                scrub.plot_histogram();
        else:
            print(f'Warning: too few cells to determine doublets! {len(b)}')
            doublet_scores = np.full(len(b), np.nan)

        doublet_vectors.append(pd.Series(data=doublet_scores, index=b.obs.index))

    doublet_vectors = pd.concat(doublet_vectors)
    adata.obs['doublet_score'] = doublet_vectors[adata.obs.index]

    return adata
