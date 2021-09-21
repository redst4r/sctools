
from anndata import AnnData
from scipy import sparse
from scipy.stats import chi2
import numpy as np


def calc_kBET(
        data: AnnData,
        attr: str,
        nn_key: 'distances'
):
    """
    TODO: the chi2 approximation is only valid for large number of neighbours.
          we instead do it on the standard 15-nn from scanpy, which might be an issue!

    :param data: anndata.AnnData of single cells, containing a neighbourhood graph
    :param attr: batch variable whose local distribution we'll check
    :param nn_key: the NN-matrix und data.obsp[nn_key] will be used

    :returns: stat, p_value; The chi2 test-statistic and pvalue of uneven
        distribution of the attribute in a cells neighborhood. Small pval indicates
        a possible batch effect (the attribute is not evenly distributed in the neighborhood)
    """
    assert attr in data.obs
    assert nn_key in data.obsp

    # build a 1-hot encoding matrix for each category/sample
    # todo: replace with sklearn.preprocessing.OneHotEncoder?? faster?
    c, mat = [], []
    for cat in data.obs[attr].unique():
        c.append(cat)
        mat.append((data.obs[attr] == cat).astype(int))

    mat = np.stack(mat,1)
    mat = sparse.csr_matrix(mat)


    ideal_dist = mat.sum(0).A.flatten()
    ideal_dist = ideal_dist/ideal_dist.sum()


    A = (data.obsp[nn_key] > 0).astype(int)  # binarized nearest neighbour matrix
    # add the diagnoal; the datapoint itself is represented in the distribution!
    # todo: does the datapoint itself really contribute??
    A = A + sparse.diags(np.ones(A.shape[0]))

    # simple matrix multiplication onto the 1-hot matrix counts the categories
    # in the neighbourhood of each datapoint
    observed_counts = A @ mat

    if False:  # this code only works if K is constant, somehow its not!!
        K = A.sum(1).A.flatten()  # number of nearest neighbors per cell; actually, this includes the cell itself (added the diagonal before)
        assert np.all(K == K[0])
        K = K[0]
        expected_counts = ideal_dist * K
    else:
        K = A.sum(1).A.flatten()
        expected_counts = (K.reshape(-1, 1) * ideal_dist)

    # Caculate the chi2 statistics
    stat = np.sum(np.square(np.subtract(observed_counts.A, expected_counts)) / expected_counts, 1)
    dof = expected_counts.shape[1] -1
    p_value = 1 - chi2.cdf(stat, dof)

    return stat, p_value
