import pandas as pd

"""
some code to parse output of Kraken/Bracken metagenomics analysis
"""

def parse_bracken(fname):
    df = pd.read_csv(fname, sep='\t')

    assert set(['name', 'taxonomy_id', 'taxonomy_lvl', 'kraken_assigned_reads', 'added_reads', 'new_est_reads', 'fraction_total_reads']) == set(df.columns)

    df['fraction_total_reads'] = df.new_est_reads / df.new_est_reads.sum()
    df['name'] = df['name'].apply(lambda x: x.strip())
    return df

def parse_kraken_report(fname):
    df = pd.read_csv(fname, sep='\t', header=None)
    df.columns = ['percent', 'reads', 'unique_reads','hierarchy', 'tax_id', 'name']

    df['name'] = df['name'].apply(lambda x: x.strip())
    return df
