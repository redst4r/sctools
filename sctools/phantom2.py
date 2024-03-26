"""
some more code around PhantomPurger, this one's based on the rust implementation
as used by pybustools.pybustools
"""

from pybustools import pybustools as rustbus
import time
import pandas as pd

def cb_overlap_in_frequencies(busfolders):
    """
    overlap of CellBarcodes across busfiles
    for each barcode, coutn the number of occurances in the busfiles (occurance is either #umi or sum of reads)

    :return:
    - df_umi_rust : DataFrame reporting the UMI abundances across samples
    - df_read_rust: DataFrame reporting the read abundances across samples
    """
    now = time.time()
    samplenames, umi_histogram, read_histogram = rustbus.phantom_fingerprint_cb(busfolders)
    elapse = time.time() - now
    print(f'took {elapse} sec')

    df_umi_rust = pd.Series(umi_histogram).to_frame().reset_index().rename({f'level_{i}': f'{sname}' for i,sname in enumerate(samplenames)}, axis=1).rename({0:'frequency'}, axis=1)
    df_read_rust = pd.Series(read_histogram).to_frame().reset_index().rename({f'level_{i}': f'{sname}' for i,sname in enumerate(samplenames)}, axis=1).rename({0:'frequency'}, axis=1)

    return df_umi_rust, df_read_rust

def cbumi_overlap_in_frequencies(busfolders):
    """
    overlap of CB/UMI across busfiles
    for each barcode/umi combination, coutn the number of occurances in the busfiles (occurance is either #umi or sum of reads)

    :return:
    - df_umi_rust : DataFrame reporting the UMI abundances across samples
    - df_read_rust: DataFrame reporting the read abundances across samples
    """
    now = time.time()
    samplenames, umi_histogram, read_histogram = rustbus.phantom_fingerprint_cbumi(busfolders)
    elapse = time.time() - now
    print(f'took {elapse} sec')
    # df_umi_rust = pd.Series(umi_histogram).to_frame().reset_index().rename({f'level_{i}': f'umi_{sname}' for i,sname in enumerate(samplenames)}, axis=1).rename({0:'frequency'}, axis=1)
    # df_read_rust = pd.Series(read_histogram).to_frame().reset_index().rename({f'level_{i}': f'read_{sname}' for i,sname in enumerate(samplenames)}, axis=1).rename({0:'frequency'}, axis=1)

    df_umi_rust = pd.Series(umi_histogram).to_frame().reset_index().rename({f'level_{i}': f'{sname}' for i,sname in enumerate(samplenames)}, axis=1).rename({0:'frequency'}, axis=1)
    df_read_rust = pd.Series(read_histogram).to_frame().reset_index().rename({f'level_{i}': f'{sname}' for i,sname in enumerate(samplenames)}, axis=1).rename({0:'frequency'}, axis=1)        
    return df_umi_rust, df_read_rust



def marginalize_fingerprint(df_fp, colnames):
    """
    project the fingerprint obtained by `cb_overlap_in_frequencies`,`cbumi_overlap_in_frequencies`
    onto the specified samples, adding up counts.
    """
    return df_fp.groupby(colnames).frequency.sum().reset_index()

import itertools
def do_pairwise(df_fp, samplenames):
    """
    pairwise comparison of fingerprints/frequencies:w
    """
    df_total = []
    for sname1,sname2 in itertools.combinations(samplenames,2):

        # df_margin = marginalize_fingerprint(df_fp, [f'umi_{sname1}', f'umi_{sname2}' ])
        # df_margin=  df_margin.rename({f'umi_{sname1}': 'count1', f'umi_{sname2}': 'count2'}, axis=1)
        df_margin = marginalize_fingerprint(df_fp, [sname1, sname2 ])
        df_margin=  df_margin.rename({sname1: 'count1', sname2: 'count2'}, axis=1)
        df_margin['sample1'] = sname1
        df_margin['sample2'] = sname2
        df_total.append(df_margin)
    df_total = pd.concat(df_total)
    return df_total
