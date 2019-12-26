import pandas as pd
import numpy as np
import scanpy as sc
import tqdm
from scipy.sparse import spmatrix
from rnaseqtools.biomart_mapping import biomart_query_all

from de import *

def annotate_gene_symbols(Q):
    df_biomart = biomart_query_all()
    var = Q.var.merge(df_biomart[['hgnc_symbol', 'ensembl_gene_id_version']].drop_duplicates('ensembl_gene_id_version'),
                      left_index=True, right_on='ensembl_gene_id_version', how='left')
    new_symbols = []
    for _, row in var.iterrows():
        sym = row['hgnc_symbol']
        if not isinstance(sym, str):
            sym=row['ensembl_gene_id_version']
        new_symbols.append(sym)
    var.hgnc_symbol = new_symbols

    var.set_index('hgnc_symbol', inplace=True)

    Q.var = var
    Q.var_names_make_unique()
    return Q


def groupby_columns(adata, groupby_field:str, aggr_fun):
    """
    aggregate the genes in the adata
    """

    assert groupby_field in adata.var or adata.var.index.name == groupby_field
    s = adata.var.groupby(groupby_field)
    x_new = []
    var_new = []
    for symbol, indices in s.indices.items():
        xx = adata.X[:, indices]
        # the_values = np.sum(xx, axis=1)
        the_values = aggr_fun(xx)

        assert the_values.shape == (adata.shape[0], 1), f"{the_values.shape} vs {(adata.shape[0], 1)}"
        x_new.append(the_values)
        var_new.append(symbol)

    adata_aggr = sc.AnnData(np.concatenate(x_new, axis=1),
                            var=pd.DataFrame(var_new, columns=[groupby_field]),
                            obs=adata.obs.copy())
    return adata_aggr


def groupby_rows(adata, groupby_field:str, aggr_fun):

    """
    aggregate the observations, i.e. to turn the singlke cell into bulk experiments
    """
    assert groupby_field in adata.obs
    s = adata.obs.groupby(groupby_field)
    x_new = []
    obs_new = []
    for observation, indices in tqdm.tqdm(s.indices.items()):
        xx = adata.X[indices, :]
        # the_values = np.sum(xx, axis=0)
        the_values = aggr_fun(xx)
        assert the_values.shape == (1, adata.shape[1]) # one entry for each obserbation

        x_new.append(the_values)
        obs_new.append(observation)

    # the indexs of obs should be string, be fefault its int however
    obs_df = pd.DataFrame(obs_new, columns=[groupby_field])
    obs_df.index = obs_df.index.astype('str')

    adata_aggr = sc.AnnData(np.concatenate(x_new),
                           obs=obs_df,
                           var=adata.var.copy())
    return adata_aggr


def adata_merge(adatas, security_check=True):
    "merging the datasets, assuming the .var are identical"

    for a in adatas:
        a.var.drop(['_id', '_score', 'notfound'], inplace=True, axis=1, errors='ignore') # ignore if droppig non-existing fields
        if 'symbol' in a.var:
            a.var['symbol'].replace({np.nan: 'Nan'}, inplace=True)

    # make sure the columns of all datasets are the same
    # sort genes alphabetically, so that each adat, has the same order
    adatas = [a[:, a.var.sort_index().index] for a in adatas]

    reference_var = adatas[0].var
    for a in adatas:
        # either equal or both nan
        assert all(reference_var.index == a.var.index) , "indices differ!"

        if security_check:  # also check that the entire var-dataframe matches exactly (sometimes doesnt due to different ensebml IDs versions)
            if not np.all(a.var == reference_var.var):
                # chekc if the mismatches are Nan
                for i, j in zip(*np.where(a.var != reference_var)):
                    assert np.isnan(a.var.values[i, j])
                    assert np.isnan(reference_var.values[i, j])

    # cannot concat if the .var dont match between adata objects

    q = adatas[0].concatenate(adatas[1:])

    q.var = reference_var.copy()
    return q


def get_cellcycle_genes():
    cell_cycle_genes =[
        'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2',
        'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2',
        'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7',
        'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1',
        'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B',
        'BRIP1', 'E2F8', 'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2',
        'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF',
        'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1',
        'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1',
        'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2',
        'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR',
        'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA'
        ]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]

    return s_genes, g2m_genes


def annotate_cellcycle(adata):
    """
    adds .obs['phase'] to the adata, indicating the estimated cell cycle phase
    also adds S_Score and G1M score to .obs
    """
    s_genes, g2m_genes = get_cellcycle_genes()
    s_genes = [_ for _ in s_genes if _ in adata.var.index]
    g2m_genes = [_ for _ in g2m_genes if _ in adata.var.index]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata

def get_diffusion_pseudotime(adata, n_dcs, min_group_size):
    dpt = sc.tools._dpt.DPT(adata,
                            n_dcs=n_dcs,
                            min_group_size=min_group_size,
                            n_branchings=0,
                            allow_kendall_tau_shift=True)

    D= np.stack([dpt.distances_dpt[_] for _ in range(len(adata))])
    return D



def annotate_coding_genes(adata):
    """
    flags all coding genes in .var
    also counts the number of reads from coding genes in .obs

    this should be run on raw counts data
    """
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

    _shared_genes = sorted(list(set(adata.var.index) & set(df_coding.index)))
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
    annotate the number of molecules/UMI, mitochondrial genes, and genes expressed
    to each cell

    this should be run on raw count data
    """
    is_sparse = isinstance(adata.X, spmatrix)
    ## some numbers on the cells
    MT_genes = [_ for _ in adata.var.index if _.startswith('MT-')]
    adata.obs['n_molecules'] = adata.X.sum(1).A.flatten() if is_sparse else adata.X.sum(1).flatten()
    adata.obs['n_mito'] = adata[:,MT_genes].X.sum(1).A.flatten() if is_sparse else adata[:,MT_genes].X.sum(1).flatten()
    adata.obs['percent_mito'] = adata.obs['n_mito'] / adata.obs['n_molecules']
    adata.obs['n_genes'] = (adata.X >0).sum(1).A.flatten()  if is_sparse else (adata.X >0).sum(1).flatten()

    adata.var['is_mito'] = adata.var.index.map(lambda x: x.startswith('MT-'))
    adata.var['is_ribo'] = adata.var.index.map(lambda x: x.startswith('RPS') or x.startswith('RPL'))

    adata.obs['n_ribo'] = adata.raw[:,adata.var.query('is_ribo==True').index].X.sum(1)
    adata.obs['percent_ribo'] = adata.obs['n_ribo'] / adata.obs['n_molecules']

    return adata
