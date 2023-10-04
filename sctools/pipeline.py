import scanpy as sc
import numpy as np
import scrublet as scr
import pandas as pd
import gc
import harmonypy as ha
import warnings
import logging
from sctools.annotations import annotate_cellcycle, annotate_qc_metrics, annotate_coding_genes
from sctools.transforms import split_adata, adata_merge
from sctools.misc import _downsample_total_counts


logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)


def _downsample(adata:sc.AnnData, split_field:str, umi_per_sample:dict):
    """
    Downsample each sample in the adata to the specified amount in umi_per_sample
    """
    adata_down = []
    for g, ad in split_adata(adata, split_field):
        target_nmol = umi_per_sample[g]
        if target_nmol != ad.raw.X.sum():
            newX = _downsample_total_counts(ad.raw.X, target_nmol, random_state=0, replace=False)
            adata_down.append(sc.AnnData(newX, obs=ad.obs, var=ad.raw.var))
        else:
            adata_down.append(sc.AnnData(ad.raw.X, obs=ad.obs, var=ad.raw.var))

    AD = adata_merge(adata_down)
    AD.raw = AD
    AD = annotate_qc_metrics(AD)
    return AD

def downsample(adata:sc.AnnData, split_field:str):
    """
    Downsample each sample in the adata to the same amount of counts (the minimum of all samples).
    recompute quality stats and return the downsampled adata
    """

    # determine counts per sample and find minimum
    nmols = {}
    for g, ad in split_adata(adata, split_field):
        nmols[g] = ad.raw.X.sum()
    min_nmol = min(nmols.values())

    print(f'downsampling to {min_nmol}')
    n_target = {g: min_nmol for g in nmols.keys()}  # set each target amount to min
    AD = _downsample(adata, split_field, n_target)
    return AD


def annotate_topN_barcodes(adata, k, split=None):
    """
    mark the k cells with the highest UMI counts per sample
    :param split: split adata according to adata.obs[split]
    """
    assert 'n_molecules' in adata.obs.columns, "n_molecules not defined"

    if split is None:  # no splitting
        ix = np.argsort(adata.obs['n_molecules'].values)[::-1][:k]
        cbs = adata.obs.index.values[ix]
    else:
        cbs = []
        for g in adata.obs[split].unique():
            df = adata.obs.query(f'{split}==@g')
            ix = np.argsort(df.n_molecules.values)[::-1][:k]
            cbs.extend(df.index.values[ix].tolist())
    adata.obs['is_cell'] = 'no'
    adata.obs.loc[cbs, 'is_cell'] ='yes'
    return adata


def standard_processing(adata, detect_doublets=True, MITO_CUTOFF=0.4, UMI_CUTOFF=1000):
    """
    Wrapper around `michi_kallisto_recipe()` with standardized values for
    QC cutoffs, does Harmony correction and doublet detection!
    """
    VARIABLE_GENES = 4000

    DE_min = 0.25  # minimum fraction of cells that must express the gene in the cluster
    DE_max = 1.0  # max % cells outside of the cluster that express the gene
    adata.raw = adata.copy()
    adata = michi_kallisto_recipe(adata,
                                  umi_cutoff=UMI_CUTOFF,
                                  n_top_genes=VARIABLE_GENES,
                                  percent_mito_cutoff=MITO_CUTOFF,
                                  annotate_cellcycle_flag=True,
                                  harmony_correction='samplename',
                                  verbose=True)

    if detect_doublets:
        logging.info('Annotating doublets')
        adata = annotate_doublets(adata, groupby='samplename')

    logging.info('Differential Experssion (batch corrected clusters)')
    # DE between the batch correted leiden clusters
    differential_expression_michi_kallisto_recipe(adata,
                                                  groupby='leiden',
                                                  n_genes=100,
                                                  method='wilcoxon',
                                                  min_in_group_fraction=DE_min,
                                                  max_out_group_fraction=DE_max,
                                                  key_added='rank_genes_groups')

    logging.info('Differential Experssion (raw (no batch correction) clusters)')
    # DE between the UNCORRECTED leiden clusters
    differential_expression_michi_kallisto_recipe(adata,
                                                  groupby='nobatch_leiden',
                                                  n_genes=100,
                                                  method='wilcoxon',
                                                  min_in_group_fraction=DE_min,
                                                  max_out_group_fraction=DE_max,
                                                  key_added='nobatch_rank_genes_groups')
    return adata


def michi_kallisto_recipe(adata, umi_cutoff=1000, n_top_genes=4000, percent_mito_cutoff=1.0,
                          annotate_cellcycle_flag=True,
                          harmony_correction=None,
                          harmony_clusters=None,
                          n_pcs=50,
                          verbose=True,
                          ignore_int_test=False):
    """
    filters for coding genes, adds QC, filters cells based on UMI, applies
    Zheng recipe calulates pca/nn/leiden/paga/umap

    adata: should be an object derived from kallisto
    """

    # make sure we get an unprocessed adata!
    assert adata.raw, ".raw does not exist!"
    if not ignore_int_test:
        assert np.all(adata.X.data == np.round(adata.X.data)), ".X must be ints/counts (maybe log1p'd?)"
        assert np.all(adata.raw.X.data == np.round(adata.raw.X.data)), ".raw.X must be ints/counts (maybe log1p'd?)"

    adata.uns['log_X'] = False  # keep track of whats logarithmized or not
    adata.uns['log_raw.X'] = False
    adata = preprocessing_michi_kallisto_recipe(adata, umi_cutoff, percent_mito_cutoff, verbose=True)

    """
    annotating the cellcycle on the reduced data (filtered for potential cells,
    not just CB) since preprocessing doenst change anything (except kicking out
    genes/cells). its alot more feasable after filtering!
    Note; this might still take ALOT of memory if more then 10k cells remain
    after filtering not sure what the cell-cycle scoring function in scanpy is
    doing internally!
    """
    if annotate_cellcycle_flag:
        logging.info('Annotating Cell cycle')
        adata = annotate_cellcycle(adata)
        logging.info('Done: Annotating Cell cycle')

    logging.info('Zheng recipe')
    sc.pp.recipe_zheng17(adata, n_top_genes=n_top_genes, log=True, plot=False, copy=False)

    logging.info('Done: Zheng recipe')
    # zheng17 log1p the .X data, but doesnt touch .raw.X!
    adata.uns['log_X'] = True

    adata = postprocessing_michi_kallisto_recipe(adata,
                                                 harmony_correction,
                                                 n_pcs=n_pcs,
                                                 harmony_clusters=harmony_clusters,
                                                 verbose=verbose)
    if not ignore_int_test:
        assert np.all(adata.raw.X.data == np.round(adata.raw.X.data)), ".raw.X got log'd!"

    return adata


def preprocessing_michi_kallisto_recipe(adata, umi_cutoff, percent_mito_cutoff, verbose=True):
    """
    filtering (coding/n_genes/n_umis/mito), annotations, but no transformations of the data
    """
    logging.info('annotating QC')
    adata = annotate_qc_metrics(adata)

    logging.info('annotating and filtering for coding genes')
    adata = annotate_coding_genes(adata)

    ix_coding = adata.var.is_coding == True
    # adata = adata[:, ix_coding]
    adata._inplace_subset_var(ix_coding) # inplace to not create a view
    gc.collect()

    logging.info('filtering cells for UMI content')
    cells_before = adata.shape[0]
    # adata = adata[adata.obs.query('n_molecules>@umi_cutoff').index]
    ix_umi = adata.obs.query('n_molecules>@umi_cutoff').index
    adata._inplace_subset_obs(ix_umi) # inplace to not create a view
    gc.collect()

    cells_after = adata.shape[0]
    logging.info(f'Cells: {cells_before} -> {cells_after}')

    logging.info('filtering cells for mito content')
    cells_before = adata.shape[0]
    ix_mito = adata.obs.query('percent_mito<@percent_mito_cutoff').index # copying to avoid getting a view, which has some issues when copying later
    adata._inplace_subset_obs(ix_mito)
    gc.collect()

    cells_after = adata.shape[0]
    logging.info(f'Cells: {cells_before} -> {cells_after}')

    return adata


def postprocessing_michi_kallisto_recipe(adata, harmony_correction, harmony_clusters=None, n_pcs=50, verbose=True):
    """
    doing batch correction, PCA/UMAP etc
    """
    logging.info('PCA')
    sc.pp.pca(adata, n_comps=n_pcs)

    gc.collect()

    if harmony_correction:
        logging.info('Harmony batch correction: uncorrected layout')
        # create an uncorrected layout first, for comparison!

        # storing the original/uncorrected PCA
        adata.obsm['X_pca_original'] = adata.obsm['X_pca'].copy()
        adata.varm['PCs_original'] = adata.varm['PCs'].copy()

        nobatch_key = 'nobatch'  # dont change this, its referenced in cellxgene export!
        logging.info('Neighbors')
        sc.pp.neighbors(adata, key_added=nobatch_key)
        logging.info('Leiden')
        sc.tl.leiden(adata, resolution=1, neighbors_key=nobatch_key, key_added=f'{nobatch_key}_leiden')
        logging.info('PAGA')
        sc.tl.paga(adata, groups=f'{nobatch_key}_leiden', neighbors_key=nobatch_key)
        sc.pl.paga(adata, color=[f'{nobatch_key}_leiden'], show=False)
        paga_init_pos = sc.tl._utils.get_init_pos_from_paga(adata, neighbors_key=nobatch_key)

        # unfort we cant tell umap where to store the result,
        # it'll go to .obsm['X_umap'] always
        logging.info('UMAP')
        sc.tl.umap(adata, init_pos=paga_init_pos, neighbors_key=nobatch_key)
        adata.obsm[f'X_umap_{nobatch_key}'] = adata.obsm['X_umap']

        # now the actual batch correction
        logging.info('Harmony batch correction: applying harmony')
        corrected_pca = _do_harmony(adata, harmony_correction, harmony_clusters)
        adata.obsm['X_pca'] = corrected_pca

    logging.info('sc.pp.neighbours')
    sc.pp.neighbors(adata)
    logging.info('sc.tl.leiden')
    sc.tl.leiden(adata, resolution=1)
    logging.info('sc.tl.paga')
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata, color=['leiden'], show=False)
    paga_init_pos = sc.tl._utils.get_init_pos_from_paga(adata)
    logging.info('sc.tl.umap')
    sc.tl.umap(adata, init_pos=paga_init_pos)
    logging.info('sc.tl.umap done')

    return adata


def _do_harmony(adata, harmony_correction: str, harmony_clusters):
    """
    apply harmony on the adata: uses the current .obsm[X_pca] to correct batch
    effects and creates a new PCA which is corrected for batches

    :params adata: AnnData object to correct
    :params harmony_correction: Batch variable
    :params harmony_clusters: Number of seeding clusters. If None, automatically selected
    :returns: np.array with corrected PC scores for each cell
    """
    # Harmony args:
    # data_mat = V, ## PCA embedding matrix of cells
    # meta_data = meta_data, ## dataframe with cell labels
    # theta = 1, ## cluster diversity enforcement
    # vars_use = 'dataset', ## variable to integrate out
    # nclust = 5, ## number of clusters in Harmony model: corresponds to the parameter K in the manuscript.
    # max.iter.harmony = 0, ## stop after initialization
    # return_object = TRUE, ## return the full Harmony model object
    # do_pca = FALSE ## don't recompute PCs
    # PHI: design matrix

    assert harmony_correction in adata.obs.columns
    n_batches = len(adata.obs[harmony_correction].unique())

    # if there's only a single sample, no batch correction needed
    # actually, harmony crashes with some error if you run it on a single batch (harmonypy 0.0.4)
    if n_batches > 1:

        vars_use = [harmony_correction]  # samplenames for harmony
        # get out the PCA matrix
        data_mat = np.array(adata.obsm['X_pca'])
        # and harmonize
        ho = ha.run_harmony(data_mat, adata.obs, vars_use,
                            max_iter_harmony=50,
                            nclust=harmony_clusters)
        corrected_X_pca = np.transpose(ho.Z_corr)
    else:
        warnings.warn('Only single batch in the data, skipping harmonh')
        corrected_X_pca = np.array(adata.obsm['X_pca'])

    return corrected_X_pca


def differential_expression_michi_kallisto_recipe(adata, groupby, n_genes=100, method='wilcoxon', min_in_group_fraction=0.0, max_out_group_fraction=1.0, use_raw=True, key_added='rank_genes_groups', csc_transform=True):
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
    :param csc_transform: for big matrices, things go alot faster if we have the .X matrix in csc format instead of csr. However this takes some memory! After we're done, we undo it
    """
    assert not adata.raw is None, "no data is present in the .raw storage. Differential expression will is only done on the .raw data!"
    assert not adata.uns['log_raw.X']
    adata.raw.X.data = np.log1p(adata.raw.X.data)

    if csc_transform:
        logging.info('doing csr->csc')
        # we have to trick scanpy a little, since we cant set .raw.X directly
        X = adata.raw.X.tocsc()
        _raw = sc.AnnData(X, obs=adata.obs, var=adata.raw.var)
        adata.raw = _raw
        logging.info('done csr->csc')


    sc.tl.rank_genes_groups(adata, groupby=groupby, n_genes=n_genes, method=method, use_raw=use_raw, key_added=key_added)
    # undoing the log
    adata.raw.X.data = np.round(np.exp(adata.raw.X.data) - 1)

    if csc_transform:  # undo it
        logging.info('doing csc->csr')
        X = adata.raw.X.tocsr()
        _raw = sc.AnnData(X, obs=adata.obs, var=adata.raw.var)
        adata.raw = _raw
        logging.info('done csc->csr')


    # filtering
    logging.info('Filtering DE')
    sc.tl.filter_rank_genes_groups(adata,
                                   key=key_added,
                                   key_added=f'{key_added}_filtered',
                                   min_in_group_fraction=min_in_group_fraction,
                                   min_fold_change=0,  # not filteirng for fold_change
                                   max_out_group_fraction=max_out_group_fraction) # not filetering for %expressing outside the cluster
    # now make the filtered genes the default DE genes
    adata.uns[f'{key_added}_unfiltered'] = adata.uns[key_added]
    adata.uns[key_added] = adata.uns[f'{key_added}_filtered']
    logging.info('Done Filtering DE')


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


def annotate_doublets(adata, groupby='samplename', PLOTTING=False, expected_doublet_rate=0.06):
    """
    run scrublet (Klein et al) on the samples within adata.
    Applies to each samplename separately (we can process merged datasets this way)

    """

    doublet_vectors = []
    for s in adata.obs[groupby].unique():
        print(s)
        b = adata[adata.obs[groupby] == s]

        if len(b) > 50:
            scrub = scr.Scrublet(b.raw.X, expected_doublet_rate=expected_doublet_rate, sim_doublet_ratio=5)
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


def annotate_celltype(adata):
    import scHCLpy.adata
    dfadata, _ = scHCLpy.adata.scHCL_adata(adata, verbose=True, n_cores=4)
    adata.obs = adata.obs.merge(dfadata, left_index=True, right_index=True)
    return adata
