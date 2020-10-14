import pandas as pd
import numpy as np
import scanpy as sc
import tqdm


def groupby_columns(adata, groupby_field:str, aggr_fun):
    """
    aggregate the genes in the adata
    """

    assert groupby_field in adata.var or adata.var.index.name == groupby_field
    s = adata.var.groupby(groupby_field)
    x_new = []
    var_new = []
    for symbol, indices in s.indices.items():
        xx = adata.X[:, indices]
        # the_values = np.sum(xx, axis=1)
        the_values = aggr_fun(xx)

        assert the_values.shape == (adata.shape[0], 1), f"{the_values.shape} vs {(adata.shape[0], 1)}"
        x_new.append(the_values)
        var_new.append(symbol)

    adata_aggr = sc.AnnData(np.concatenate(x_new, axis=1),
                            var=pd.DataFrame(var_new, columns=[groupby_field]),
                            obs=adata.obs.copy())
    return adata_aggr


def groupby_rows(adata, groupby_field:str, aggr_fun):

    """
    aggregate the observations, i.e. to turn the singlke cell into bulk experiments
    """
    assert groupby_field in adata.obs
    s = adata.obs.groupby(groupby_field)
    x_new = []
    obs_new = []
    for observation, indices in tqdm.tqdm(s.indices.items()):
        xx = adata.X[indices, :]
        # the_values = np.sum(xx, axis=0)
        the_values = aggr_fun(xx)
        assert the_values.shape == (1, adata.shape[1]) # one entry for each obserbation

        x_new.append(the_values)
        obs_new.append(observation)

    # the indexs of obs should be string, be fefault its int however
    obs_df = pd.DataFrame(obs_new, columns=[groupby_field])
    obs_df.index = obs_df.index.astype('str')

    adata_aggr = sc.AnnData(np.concatenate(x_new),
                            obs=obs_df,
                            var=adata.var.copy())
    return adata_aggr


def adata_merge(adatas, security_check=True, verbose=True):
    "merging the datasets, assuming the .var are identical"

    assert isinstance(adatas, list), "adatas hsa to be supplied as list!"

    # TODO get rid of this custom code
    for a in adatas:
        a.var.drop(['_id', '_score', 'notfound'], inplace=True, axis=1, errors='ignore') # ignore if droppig non-existing fields
        if 'symbol' in a.var:
            a.var['symbol'].replace({np.nan: 'Nan'}, inplace=True)

    # make sure the columns of all datasets are the same
    # sort genes alphabetically, so that each adat, has the same order
    adatas = [a[:, a.var.sort_index().index] for a in adatas]

    reference_var = adatas[0].var
    for a in adatas:
        # either equal or both nan
        assert all(reference_var.index == a.var.index) , "indices differ!"

        if security_check:  # also check that the entire var-dataframe matches exactly (sometimes doesnt due to different ensebml IDs versions)
            if not np.all(a.var == reference_var):
                # chekc if the mismatches are Nan
                for i, j in zip(*np.where(a.var != reference_var)):
                    assert np.isnan(a.var.values[i, j]), f"{i},{j}, {a.var.loc[i,j]}"
                    assert np.isnan(reference_var.values[i, j])

    # cannot concat if the .var dont match between adata objects

    q = adatas[0].concatenate(adatas[1:])

    q.var = reference_var.copy()
    return q
