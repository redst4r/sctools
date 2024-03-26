import scanpy as sc
import numpy as np
import pandas as pd
from sctools import kallisto 
# from sctools import adata_merge, annotate_cellcycle, enrich, annotate_qc_metrics

def load_from_kallisto(h5_filename: str):
    """
    load the output of kallisto-velocity (which puts out an anndata.AnnData)
    """
    Q = sc.read_h5ad(h5_filename)
    Q = kallisto.annotate_gene_symbols(Q)
    return Q

def annotate_velo(A):
    """
    get the output of kallisto-velocity ready for scanpy
    - renaming layers
    - annotate some metadata to .obs
    """
    # scVelo assumes the layers to be named
    # - 'unspliced'
    # - 'spliced'
    A.layers['spliced'] = A.layers['mature']
    A.layers['unspliced'] = A.layers['nascent']

    ## something weird: all the MT-genes have zero counts in unspliced or spliced layers,
    # they're all `ambigous`
    n_molecules_unspliced = A.layers['unspliced'].sum(1).A1
    n_molecules_spliced = A.layers['spliced'].sum(1).A1
    n_molecules_ambiguous = A.layers['ambiguous'].sum(1).A1

    A.obs['nuclear_fraction'] = n_molecules_unspliced / (n_molecules_unspliced + n_molecules_spliced)
    A.obs['n_molecules_unspliced'] = n_molecules_unspliced
    A.obs['n_molecules_spliced'] = n_molecules_spliced
    A.obs['n_molecules_ambiguous'] = n_molecules_ambiguous
    A.obs['n_molecules_total'] = n_molecules_ambiguous + n_molecules_unspliced + n_molecules_spliced

    # ix_mito = A.var.query('is_mito').index
    ix_mito = [_ for _ in A.var_names if _.startswith('MT-')]

    A.obs['n_mito'] = A[:, ix_mito].layers['ambiguous'].sum(1).A1 + A[:, ix_mito].layers['spliced'].sum(1).A1 + A[:, ix_mito].layers['unspliced'].sum(1).A1
    A.obs['percent_mito'] = A.obs['n_mito'] / A.obs['n_molecules_total']
    return A


def load_kallisto_and_velocity(fname_kallisto, fname_velocity):
    """
    jointly loads the kallisto qunaitifcation and the corresponding RNA-velocity.
    (note that while RNA-velocity is also an AnnData, it is slightly different from the kallisto qunaitifaction)
    """

    adata = kallisto.load_from_kallisto(fname_kallisto, kallisto_prefix='gene')

    # load the velo object too
    adata_velo = load_from_kallisto(fname_velocity)
    adata_velo = annotate_velo(adata_velo)

    # pull in the nuclear fraction
    adata.obs = adata.obs.merge(
        adata_velo.obs['nuclear_fraction'], left_index=True, right_index=True, how='left'
    )

    # carry over the velocity layers
    # note: the overlap of cells is not 100% (theres some cells in adata that dont have a "sister" in adata_velo)
    # but the overlap should be substantial!
    # TODO check overlap size to rule out bugs

    shared_cells = sorted(list(set(adata.obs_names ) & set(adata_velo.obs_names)))
    adata = adata[shared_cells]
    adata_velo = adata_velo[shared_cells]

    adata.layers['unspliced'] = adata_velo.layers['unspliced']
    adata.layers['spliced'] = adata_velo.layers['spliced']
    adata.layers['ambiguous'] = adata_velo.layers['ambiguous']

    return adata

