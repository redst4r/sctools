import h5py
from scipy import sparse
import pandas as pd
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from anndata import AnnData
import rpy2.robjects as ro
import scanpy as sc
import numpy as np

"""
some code around getting kallisto-bulk quantificcations into adatas using R's tximport
"""

def read_kallisto_h5(h5file, samplename):
    """
    along the lines of https://github.com/thelovelab/tximport/blob/7b1fbe4bffc98c21e891a37c5a0e1751a85045ba/R/helper.R#L79
    """
    with h5py.File(h5file) as f:
        transcript_ids = f['/aux/ids'][:]
        eff_lengths = f['/aux/eff_lengths'][:]
        lengths = f['/aux/lengths'][:]
        counts = sparse.csr_matrix(f['est_counts'])
    
        assert len(transcript_ids) == counts.shape[1]
        assert len(transcript_ids) == len(eff_lengths)
        
        n_bootstraps = f['/aux/num_bootstrap'][0]
    
        bootstraps = []
        for i in range(n_bootstraps):
            # read the coutns into a sparse mat
            counts_bs = sparse.csr_matrix(
                f[f'/bootstrap/bs{i}'][:]
            )
            bootstraps.append(counts_bs)
        bootstraps = sparse.vstack(bootstraps)
    
        normfac = (1e6)/np.sum(counts.A.flatten()/eff_lengths)
        tpm = normfac * (counts.A.flatten() / eff_lengths)
        # assert 1==0
        var = pd.DataFrame({
            'target_id': transcript_ids,
            'eff_length': eff_lengths,
            'lengths': lengths,
            'est_counts': counts.A.flatten(),
            'tpm': tpm
        })
        var.target_id =  var.target_id.apply(lambda x: x.decode())
        X = counts
        X_bs = bootstraps

        """
        uses (rounded) est_counts as entires in the count matrix
        """
        adata = sc.AnnData(
            np.round(var['est_counts'].values).astype(int).reshape(1,-1),
            var = var[['target_id','eff_length','lengths']].set_index('target_id'),
            obs=pd.DataFrame([{'samplename': samplename}]).set_index('samplename'),
            layers={'tpm': var['tpm'].values.reshape(1, -1)}
        )
    
    return X_bs, adata


def read_kallisto_h5_via_tximport(design_file, formula, design_column_filepath='path', txout=True, t2g=None):
    """
    reads an entire design using tximport in R
    returns an AnnData

    if `txout = False`, we aggregate to gene level and need to supply a transcript2gene file
    (e.g. `/home/michi/mounts/TB4drive/kallisto_resources_v50/t2g.txt`

    :param design_column_filepath: which column in design_file contains the filepath info
    """
    assert formula.startswith('~')
    if not txout:
        assert t2g is not None
    ## prep
    ro.r(
        f"""
        library(tximport)
        library(dplyr)
        s2c <- read.table('{design_file}', header=TRUE, sep='\t')
        s2c <- dplyr::mutate(s2c, deseq_path = file.path({design_column_filepath}, 'abundance.h5'))  # for deseq, we need the complete path to the h5 files
        """
    )

    ## tximport
    if txout:
        ro.r(
        """       
        xxx = tximport(s2c$deseq_path, type="kallisto", txOut=TRUE, ignoreAfterBar=TRUE)
        # xx contains:
        # - xxx$abundance
        # - xxx$counts
        # - xxx$infReps
        # - xxx$length
        # - xxx$countsFromAbundance 
        # each of which is a genes x samples table
        """
        )
    else:
        # mapping to genes instead
        ro.r(
        f"""
        # load the t2g
        t2gtable = read.table('{t2g}', header=FALSE, sep='\t')
        
        xxx = tximport(s2c$deseq_path, type="kallisto", txOut=FALSE, ignoreAfterBar=TRUE, tx2gene=t2gtable)
        # xx contains:
        # - xxx$abundance
        # - xxx$counts
        # - xxx$infReps
        # - xxx$length
        # - xxx$countsFromAbundance 
        # each of which is a genes x samples table
        """
        )
        
    ## make the DEseq2 object
    ro.r(
    f"""
    library(DESeq2)
    dds <- DESeqDataSetFromTximport(xxx,
                                    colData = s2c,
                                    design = {formula})
    """
    )   
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        df = ro.r("as.data.frame(assay(dds))")
        obs = ro.r("as.data.frame(colData(dds))")
    
    adata = AnnData(
        X=sparse.csr_matrix(df.T.values),
        obs=obs,
        var=pd.DataFrame(index=df.index)
    )    
    return adata