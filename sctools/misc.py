import scanpy as sc
import numpy as np

class Verbose(object):
    """
    context manager for verbose output around statements
    """
    def __init__(self, msg, verbose):
        self.msg = msg
        self.verbose = verbose
    def __enter__(self):
        if self.verbose:
            print(self.msg)
        return
    def __exit__(self, type, value, traceback):
        if self.verbose:
            print(f'Done {self.msg}')



def get_diffusion_pseudotime(adata, n_dcs, min_group_size):
    dpt = sc.tools._dpt.DPT(adata,
                            n_dcs=n_dcs,
                            min_group_size=min_group_size,
                            n_branchings=0,
                            allow_kendall_tau_shift=True)

    D = np.stack([dpt.distances_dpt[_] for _ in range(len(adata))])
    return D
