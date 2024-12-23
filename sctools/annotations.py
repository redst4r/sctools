import pandas as pd
import numpy as np
from scipy.sparse import spmatrix
from rnaseqtools.biomart_mapping import biomart_query_all
from sctools.score_genes import score_genes_cell_cycle, get_cellcycle_genes


def annotate_cellcycle(adata):
    """
    adds .obs['phase'] to the adata, indicating the estimated cell cycle phase
    also adds S_Score and G1M score to .obs
    """
    s_genes, g2m_genes = get_cellcycle_genes()
    s_genes = [_ for _ in s_genes if _ in adata.var.index]
    g2m_genes = [_ for _ in g2m_genes if _ in adata.var.index]
    # use my internal, memory saving method for now
    score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata


def annotate_coding_genes(adata):
    """
    flags all coding genes in .var
    also counts the number of reads from coding genes in .obs

    this should be run on raw counts data
    """

    if 'is_coding' in adata.var.columns:
        print('Coding genes already annotated. skipping this!')
        return adata

    df_biomart = biomart_query_all()

    def _is_coding(x):
        coding_transcript_biotype = [
            'protein_coding',
            'TR_D_gene',
            'TR_C_gene',
            'IG_J_gene',
            'IG_C_gene',
            'IG_D_gene',
            'TR_J_gene',
            'TR_V_gene',
            'IG_V_gene',
        ]
        # check if biotypes of the gene intersect with coding
        return len(set(x) & set(coding_transcript_biotype)) > 0

    df_coding = pd.DataFrame(df_biomart.groupby('hgnc_symbol')['transcript_biotype'].apply(_is_coding)).rename({'transcript_biotype': 'is_coding'}, axis=1)

    # _shared_genes = sorted(list(set(adata.var.index) & set(df_coding.index)))
    # adata = adata[:, _shared_genes] # often some genes are not in df_biomart
    adata.var = adata.var.merge(df_coding, left_index=True, right_index=True, how='left')
    # adata.raw.var= adata.raw.var.merge(df_coding, left_index=True, right_index=True)

    adata.obs['n_coding_molecules'] = np.array(adata.raw[:, adata.var.is_coding==True].X.sum(1)).flatten()
    adata.obs['percent_coding_molecules'] = adata.obs['n_coding_molecules'] / np.array(adata.X.sum(1)).flatten()
    # weirdly, sometimes barcodes have n_molecules==0, to avoid Nan, expliciatly set to 0
    adata.obs.loc[adata.obs['n_molecules']==0, 'percent_coding_molecules'] = 0

    return adata


def annotate_qc_metrics(adata):
    """
    annotate the number of molecules/UMI, mitochondrial genes, and genes
    expressed to each cell

    this should be run on raw count data
    """
    is_sparse = isinstance(adata.X, spmatrix)
    # some numbers on the cells
    adata.obs['n_molecules'] = adata.X.sum(1).A.flatten() if is_sparse else adata.X.sum(1).flatten()
    adata.obs['n_genes'] = (adata.X >0).sum(1).A.flatten() if is_sparse else (adata.X >0).sum(1).flatten()

    # annotating mitochondiral content
    """
    for some reason `adata[:,MT_genes].X.sum(1)` creates a memory leak, (also some pandas warnings about is_categorical())
    however, using the .raw version works fine?!
    """

    adata.var['is_mito'] = adata.var.index.map(lambda x: x.startswith('MT-'))

    if len(adata.var.query('is_mito==True')) > 0:
        # MT_genes = [_ for _ in adata.var.index if _.startswith('MT-')]
        # # memory leak:
        # adata.obs['n_mito'] = adata[:,MT_genes].X.sum(1).A.flatten() if is_sparse else adata[:,MT_genes].X.sum(1).flatten()

        # this works fine, due to .raw
        ix = adata.var.query('is_mito==True').index
        adata.obs['n_mito'] = adata.raw[:,ix].X.sum(1).A.flatten() if is_sparse else adata.raw[:,ix].X.sum(1).flatten()

    else:
        # in case the dat adoesnt contain a single mitochondrial gene.
        # Unless we specifically filtered for this, that shoul never happen!!
        adata.obs['n_mito'] = 0

    adata.obs['percent_mito'] = adata.obs['n_mito'] / adata.obs['n_molecules']

    # annotating robosomal content
    adata.var['is_ribo'] = adata.var.index.map(lambda x: x.startswith('RPS') or x.startswith('RPL'))
    if len(adata.var.query('is_ribo==True')) > 0:
        adata.obs['n_ribo'] = adata.raw[:,adata.var.query('is_ribo==True').index].X.sum(1).A.flatten()
    else:
        adata.obs['n_ribo'] = 0
    adata.obs['percent_ribo'] = adata.obs['n_ribo'] / adata.obs['n_molecules']

    nmg = [_ for _ in nuclear_mito_genes if _ in adata.var_names]
    if len(nmg) > 0:
        adata.obs['n_nuclear_mito'] = adata[:,nmg].X.sum(1).A.flatten() if is_sparse else adata[:,nmg].X.sum(1).flatten()
    else:
        adata.obs['n_nuclear_mito'] = 0 
    adata.obs['percent_nuclear_mito'] = adata.obs['n_nuclear_mito'] / adata.obs['n_molecules']

    return adata


nuclear_mito_genes = [
    # complex 1
    'NDUFS7', 'NDUFS8', 'NDUFV2', 'NDUFS3', 'NDUFS2', 'NDUFV1', 'NDUFS1', 'NDUFS6', 'NDUFA12', 'NDUFS4', 'NDUFA9', 'NDUFAB1', 'NDUFA2', 'NDUFA1', 'NDUFB3', 'NDUFA5', 'NDUFA6', 'NDUFA11', 'NDUFB11', 'NDUFS5', 'NDUFB4', 'NDUFA13', 'NDUFB7', 'NDUFA8', 'NDUFB9', 'NDUFB10', 'NDUFB8', 'NDUFC2', 'NDUFB2', 'NDUFA7', 'NDUFA3', 'NDUFA4', 'NDUFB5', 'NDUFB1', 'NDUFC1', 'NDUFA10', 'NDUFA4L2', 'NDUFV3', 'NDUFB6', 'NDUFAF1', 'NDUFAF2', 'NDUFAF3', 'NDUFAF4',
    # complex2
    'SDHA', 'SDHB', 'SDHC', 'SDHD',

    # complex 3
    'ETFDH',

    #complex 4
    'CYC1', 'UQCRFS1', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRB', 'UQCRH', 'UQCRQ', 'UQCR10', 'UQCR11',

    #complex 5
    'COX4I1', 'COX4I2', 'COX5A', 'COX5B', 'COX6A1', 'COX6A2', 'COX6B1', 'COX6B2', 'COX6C', 'COX7A1', 'COX7A2', #'COX7A3', 'COX7B', 'COX7C', 'COX7A2L', 'COX8A', 'COX8C', 'COA1', 'COA3', 'COA4', 'COA5', 'COA6', 'COA7', 'COX11', 'COX14', 'COX15', 'COX16', 'COX17', 'COX18', 'COX19', 'COX20'
]
