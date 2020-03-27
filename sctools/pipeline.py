import scanpy as sc
from sctools import annotate_qc_metrics, annotate_cellcycle, annotate_coding_genes
import numpy as np

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

def standard_processing(adata):
    """
    Wrapper around `michi_kallisto_recipe()` with standardized values for
    QC cutoffs
    """
    UMI_CUTOFF = 1000
    VARIABLE_GENES = 4000
    MITO_CUTOFF = 0.4

    adata.raw = adata.copy()
    adata = michi_kallisto_recipe(adata,
                                  umi_cutoff=UMI_CUTOFF,
                                  n_top_genes=VARIABLE_GENES,
                                  percent_mito_cutoff=MITO_CUTOFF,
                                  annotate_cellcycle_flag=True,
                                  verbose=True )
    return adata

def michi_kallisto_recipe(adata, umi_cutoff=1000, n_top_genes=4000, percent_mito_cutoff=1, annotate_cellcycle_flag=False, verbose=True):

    """
    filters for coding genes, adds QC, filters cells based on UMI, applies
    Zheng recipe calulates pca/nn/leiden/paga/umap

    adata: should be an object derived from kallisto
    """

    # make sure we get an unprocessed adata!
    assert adata.raw, ".raw does not exist!"
    assert np.all(adata.X.data== np.round(adata.X.data)), ".X must be ints/counts (maybe log1p'd?)"
    assert np.all(adata.raw.X.data== np.round(adata.raw.X.data)), ".raw.X must be ints/counts (maybe log1p'd?)"

    adata.uns['log_X'] = False # keep track of whats logarithmized or not
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
        adata = annotate_cellcycle(adata)

    if verbose:
        print('Zheng recipe')

    # sc.pp.filter_genes(adata, min_counts=1)
    # sc.pp.normalize_total(adata, key_added='n_counts_all')

    sc.pp.recipe_zheng17(adata, n_top_genes=n_top_genes, log=True, plot=True, copy=False)

    # zheng17 log1p the .X data, but doesnt touch .raw.X!
    adata.uns['log_X'] = True

    adata = postprocessing_michi_kallisto_recipe(adata, verbose)

    assert np.all(adata.raw.X.data== np.round(adata.raw.X.data)), ".raw.X got log'd!"


    return adata


def preprocessing_michi_kallisto_recipe(adata, umi_cutoff, percent_mito_cutoff, verbose=True):

    """
    filtering (coding/n_genes/n_umis/mito), annotations, but no transformations of the data
    """
    adata = annotate_qc_metrics(adata)
    adata = annotate_coding_genes(adata)
    adata = adata[:, adata.var.is_coding==True]

    if verbose:
        print('filtering cells for UMI content')
    cells_before = adata.shape[0]
    adata = adata[adata.obs.query('n_molecules>@umi_cutoff').index]
    cells_after = adata.shape[0]
    if verbose:
        print(f'Cells: {cells_before} -> {cells_after}')

    if verbose:
        print('filtering cells for mito content')
    cells_before = adata.shape[0]
    adata = adata[adata.obs.query('percent_mito<@percent_mito_cutoff').index]
    cells_after = adata.shape[0]
    if verbose:
        print(f'Cells: {cells_before} -> {cells_after}')

    return adata

def postprocessing_michi_kallisto_recipe(adata, verbose=True):
    """
    doing PCA/UMAP etc
    """
    if verbose:
        print('PCA')
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=1)
    sc.tl.louvain(adata)
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata, color=['leiden'] ,show=False)

    paga_init_pos = sc.tl._utils.get_init_pos_from_paga(adata)

    sc.tl.umap(adata, init_pos=paga_init_pos)

    # sc.tl.draw_graph(adata, init_pos='paga') # somehow the init_pos works different here

    return adata


def differential_expression_michi_kallisto_recipe(adata, groupby, n_genes=100):
    """
    scanpy by default runs the differential expression on .raw.X, but also assumes
    logarithmized data. Otherwise, pvals come out wrong!

    To leave .raw unchanged, we log-transform, do DE-analysis and undo the log_transform

    TODO this better be an atomic operation, otherwise the log-flag gets screwed up
    or stick it into a decorator!
    """

    assert not adata.uns['log_raw.X']
    adata.raw.X.data = np.log1p(adata.raw.X.data)
    sc.tl.rank_genes_groups(adata, groupby=groupby, n_genes=n_genes)
    adata.raw.X.data = np.round(np.exp(adata.raw.X.data) - 1)



def export_for_cellxgene(adata, annotations):
    # for Export we have to pull all the genes back into the adata.X!
    _tmp = sc.AnnData(adata.raw.X, obs=adata.obs, var=adata.raw.var)
    _tmp.uns = adata.uns
    _tmp.obsm = adata.obsm
    _tmp.obs =  _tmp.obs[annotations]
    return _tmp
