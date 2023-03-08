from rnaseqtools import biomart_mapping
import scanpy as sc
import pandas as pd
import json

def annotate_gene_symbols(Q):

    df_biomart = biomart_mapping.biomart_query_all()

    """
    sometimes the ensembl_gene_id_version doesnt match up with what kallisto puts down
    (maybe a newr/older version of the transcriptome than used in biomarkt)
    hence lets use the UNVERSIONED id
    """

    # some weirdness usingthe kallisto-velocity transcriptome-index
    Q.var['ensembl_gene_id'] = Q.var.index.map(lambda x: x.split('.')[0] if not x.endswith('_PAR_Y') else x.split('.')[0]+'_PAR_Y')

    assert len(Q.var) == len(Q.var['ensembl_gene_id'].unique())

    with_version = False
    ensebml_id_field = 'ensembl_gene_id_version' if with_version else 'ensembl_gene_id'
    df_tmp = df_biomart[['hgnc_symbol', ensebml_id_field]].drop_duplicates(ensebml_id_field)

    if with_version:
        var = Q.var.merge(df_tmp,
                          left_index=True,
                          right_on=ensebml_id_field,
                          how='left')
    else:
        var = Q.var.merge(df_tmp,
                          left_on='ensembl_gene_id',
                          right_on=ensebml_id_field,
                          how='left')

    new_symbols = []
    for _, row in var.iterrows():
        sym = row['hgnc_symbol']
        if not isinstance(sym, str):
            sym = row[ensebml_id_field]
        new_symbols.append(sym)
    var.hgnc_symbol = new_symbols

    var.set_index('hgnc_symbol', inplace=True)

    Q.var = var
    Q.var_names_make_unique()
    return Q


def load_from_kallisto(folder: str, kallisto_prefix='genecount'):
    """
    kallisto_prefix: usually the kallisto files are named
    - genecount.mtx
    - genecount.barcodes.txt
    - genecount.genes.txt

    but we can change that prefix
    """
    mtx_file = f'{folder}/{kallisto_prefix}.mtx'
    obs_file = f'{folder}/{kallisto_prefix}.barcodes.txt'
    var_file = f'{folder}/{kallisto_prefix}.genes.txt'

    Q = _load_kallisto_to_adata(mtx_file, var_file, obs_file, metricsfile=None)
    Q = annotate_gene_symbols(Q)
    return Q


def _load_kallisto_to_adata(matrixfile, genefile, bcfile, metricsfile=None):
    """
    basic form of loading kallisto files into anndata.AnnData objects
    - note that this does NOT annotate the genes
    """
    _tmp = sc.read_mtx(matrixfile)
    obs = pd.read_csv(bcfile,   header=None, index_col=0)
    var = pd.read_csv(genefile, header=None, index_col=0)
    obs.index.name = 'CB'  # needed for anndata>=0.7 which doesnt allow int as name
    var.index.name = None  # needed for anndata>=0.7 which doesnt allow int as name

    _tmp.obs = obs
    _tmp.var = var
    if metricsfile:
        metric_dict = pd.read_json(metricsfile, orient='columns', typ='series').to_dict()
        _tmp.uns['kallisto_metrics'] = metric_dict
    return _tmp


def parse_kallisto_log(filename):
    with open(filename,'r') as fh:
        l = ''.join(fh.readlines())
    j = json.loads(l)
    return j