import scanpy as sc
import numpy as np
import numba
# import scipy as sp
from scipy.sparse import issparse, isspmatrix_csr, csr_matrix
# from sklearn.utils import sparsefuncs, check_array
# from pandas.api.types import is_categorical_dtype
# from anndata import AnnData
import h5py
from anndata._io.h5ad import read_dataframe


def get_diffusion_pseudotime(adata, n_dcs, min_group_size):
    dpt = sc.tools._dpt.DPT(adata,
                            n_dcs=n_dcs,
                            min_group_size=min_group_size,
                            n_branchings=0,
                            allow_kendall_tau_shift=True)

    D = np.stack([dpt.distances_dpt[_] for _ in range(len(adata))])
    return D




"""
adata converts .X to float32 ALWAYS.
this causes numerical issues when doing the .cumsum() below

i.e. this sum can come out smaller then the actual sum
"""


def _downsample_total_counts(X, total_counts, random_state, replace):
    total_counts = int(total_counts)
    total = X.sum()
    if total < total_counts:
        return X
    assert issparse(X)
    original_type = type(X)
    if not isspmatrix_csr(X):
        X = csr_matrix(X)

    X = X.astype(int)
    _downsample_array(
        X.data,
        total_counts,
        random_state=random_state,
        replace=replace,
        inplace=True,
    )
    X.eliminate_zeros()
    if original_type is not csr_matrix:
        X = original_type(X)

    return X


@numba.njit(cache=True)
def _downsample_array(
    col: np.ndarray,
    target: int,
    random_state = 0,
    replace: bool = True,
    inplace: bool = False,
):
    """\
    Evenly reduce counts in cell to target amount.
    This is an internal function and has some restrictions:
    * total counts in cell must be less than target
    """
    np.random.seed(random_state)
    cumcounts = col.cumsum()
    if inplace:
        col[:] = 0
    else:
        col = np.zeros_like(col)
    total = np.int_(cumcounts[-1])
    sample = np.random.choice(total, target, replace=replace)
    sample.sort()
    geneptr = 0
    for count in sample:
        while count >= cumcounts[geneptr]:
            geneptr += 1
        col[geneptr] += 1
    return col


def load_obs(h5ad_filename):
    """
    loads the obs-dataframe only
    """
    with h5py.File(h5ad_filename, 'r') as f:
        obs = read_dataframe(f['/obs'])
    return obs

def load_dataframe(h5ad_filename, path):
    """
    loads any dataframe in the h5ad. path specifics the location within the
    h5ad (e.g. /raw/var, /obs, /var)
    """
    with h5py.File(h5ad_filename, 'r') as f:
        df = read_dataframe(f[path])
    return df
