from rnaseqtools import biomart_mapping
import scanpy as sc
import pandas as pd


def annotate_gene_symbols(Q):

    df_biomart=biomart_mapping.biomart_query_all()

    """
    sometimes the ensembl_gene_id_version doesnt match up with what kallisto puts down
    (maybe a newr/older version of the transcriptome than used in biomarkt)
    hence lets use the UNVERSIONED id
    """
    Q.var['ensembl_gene_id'] = Q.var.index.map(lambda x: x.split('.')[0])
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

def load_from_kallisto(folder:str, kallisto_prefix='genecount'):
    """
    kallisto_prefix: usually the kallisto files are named
    - genecount.mtx
    - genecount.barcodes.txt
    - genecount.genes.txt

    but we can change that prefix
    """
    mtx_file =  f'{folder}/{kallisto_prefix}.mtx'
    obs_file =  f'{folder}/{kallisto_prefix}.barcodes.txt'
    var_file =  f'{folder}/{kallisto_prefix}.genes.txt'
    Q = sc.read_mtx(mtx_file)
    obs = pd.read_csv(obs_file, header=None, index_col=0)
    var = pd.read_csv(var_file, header=None, index_col=0)

    Q.obs = obs
    Q.var = var

    Q = annotate_gene_symbols(Q)
    return Q
