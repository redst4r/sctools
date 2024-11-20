import pandas as pd
import numpy as np
import anndata2ri
from rpy2.robjects import r
import plotnine as pn


def emptyDroplets_Rwrapper(adata):
    """
    running the emptyDrops R package to determine which barcodes correspodn to cells
    and which ones are just empty

    :params adata: unprocessed AnnData from a single 10x experiment. This should be UNFILTERED for barcodes (whihc means it will be usually severel millions of "cells")

    :returns:

    """
    anndata2ri.activate()
    r('library(DropletUtils)')
    x = anndata2ri.py2rpy(adata)
    counts = r('assay')(x)  # the sparse matrix of counts
    eout = r('emptyDrops')(counts)
    anndata2ri.deactivate()

    eout['is_cell'] = eout.FDR < 0.001
    eout['Total'] = eout['Total'].astype(float)
    df_tmp = eout[~np.isnan(eout.LogProb)]
    p = pn.ggplot(df_tmp, pn.aes(x='Total', y='-LogProb', color='is_cell')) + \
        pn.geom_point(size=1, alpha=0.5) + \
        pn.scale_x_log10() + pn.scale_y_log10()

    return eout, p


def estimate_ambient_profile(adata):
    r('library(edgeR)')
    anndata2ri.activate()
    x = anndata2ri.py2rpy(adata)
    goodTuringProportions = r('goodTuringProportions')
    counts = r('assay')(x)  # the sparse matrix of counts
    rsum = r('rowSums')(counts)
    amb = goodTuringProportions(rsum)
    ambient_df = pd.DataFrame(amb.flatten(), index=adata.var.index, columns=['ambient_fraction'])

    anndata2ri.deactivate()
    return ambient_df



