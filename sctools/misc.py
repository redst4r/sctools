import scanpy as sc
import numpy as np

def get_diffusion_pseudotime(adata, n_dcs, min_group_size):
    dpt = sc.tools._dpt.DPT(adata,
                            n_dcs=n_dcs,
                            min_group_size=min_group_size,
                            n_branchings=0,
                            allow_kendall_tau_shift=True)

    D = np.stack([dpt.distances_dpt[_] for _ in range(len(adata))])
    return D
