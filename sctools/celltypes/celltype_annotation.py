import numpy as np
from sklearn.metrics import roc_auc_score
import pandas as pd
import tqdm
import toolz
import sctools.celltypes.panglaodb as panglaodb


def build_marker_dict(adata):
    """
    my custom list of markers, partially based on panglaodb, partially on personal communication
    """

    celltype_marker_dict = panglaodb.get_celltype_markers(
        sens=0.8, spec=0.8,
        list_of_organs=['Connective tissue', 'Immune system',
                        'GI tract', 'Blood', 'Vasculature',
                        'Smooth muscle', 'Epithelium', 'Skin']
        )

    # now fix a few by hand and expert knowledge
    celltype_marker_dict['Gastric chief cells'] = ['PGA4', 'PGA3', 'PGA5', 'PGC']
    celltype_marker_dict['Enteroendocrine cells'].extend(['CHGA', 'CHGB'])
    celltype_marker_dict['Parietal cells'] = ['ATP4A', 'ATP4B']
    celltype_marker_dict['Endothelial cells'].append('PECAM1')
    celltype_marker_dict['Goblet cells'] = ['MUC2','ITLN1', 'SPINK4']
    celltype_marker_dict['Intestinal epithelial cells'] = ['CDX2']
    celltype_marker_dict['Foveolar cells 1'] = ['MUC5AC', 'TFF1' ]
    celltype_marker_dict['Foveolar cells 2'] = ['TFF1' ,'GKN1']
    celltype_marker_dict['Mucus neck cells'] = ['TFF1' ,'MUC6']
    celltype_marker_dict['MSC'] = ['EPHB2', 'SOX9','OLFM4','CD44','AXIN2', 'SOX4']
    celltype_marker_dict['Paneth cells'] = ['DEFA5', 'DEFA6']

    # delete a few that we dont see anyways
    del celltype_marker_dict['Platelets']
    del celltype_marker_dict['Neutrophils']
    del celltype_marker_dict['Pericytes'] # <-?? maybe should be included
    del celltype_marker_dict['Adipocytes']
    del celltype_marker_dict['Eosinophils']

    #filter out all genes that we dont measure in adata
    celltype_marker_dict = toolz.valmap(lambda v: [_ for _ in v if _ in adata.raw.var_names], celltype_marker_dict)

    empty_celltypes = [k for k,v in celltype_marker_dict.items() if len(v) == 0]
    if empty_celltypes:
        print('No markers for:')
        print(empty_celltypes)

    celltype_marker_dict = toolz.valfilter(lambda v: len(v) > 0, celltype_marker_dict)

    return celltype_marker_dict

def dict_argmax_val(d):
    "returns the key with the biggest value"
    max_val = -np.inf
    max_key = None
    for k,v in d.items():
#         print(v, max_val)
        if v > max_val:
            max_key = k
            max_val = v
    return max_key, max_val

# def create_meta_markers_pairwise(adata, marker_dict, cluster_name='leiden', auc_cutoff=0.7, SIMPLE=True):
#     assert adata.raw
#
#     # first, erase all previous marker annotations
#     to_drop = [_ for _ in adata.obs.columns if _.startswith('broad_score_')] + ['marker_based_celltype']
#     adata.obs.drop(to_drop, axis=1, inplace=True, errors='ignore')
#
#     mol_per_cell = np.array(adata.raw.X.sum(1)).flatten()
#     for celltype, markers in marker_dict.items():
#         if len(markers) == 1:
#             metagene = np.array(adata.raw[:,markers].X).flatten()
#         else:
# #             metagene = np.array(adata.raw[:,markers].X.sum(1)).flatten()
#             metagene = np.array(adata.raw[:,markers].X.mean(1)).flatten()
#
#         s = np.log(1e-5 + 10**4 * metagene/mol_per_cell)
#         adata.obs[f'broad_score_{celltype}'] = s
#
#     # use the AUC to figure out which celltype each cluster belongs to
#     # for each cluster look which cell-score has the best AUC
#
#     auc_dict = {}
#     df = []
#     import itertools
#     cluster_pairs = itertools.combinations(adata.obs[cluster_name].unique(), 2)
#     for cl1, cl2 in tqdm.tqdm(cluster_pairs):
#         celltypes = marker_dict.keys()
#         ix = adata.obs.query(f'{cluster_name} == @cl1 or {cluster_name}==@cl2').index
#         adata_sub = adata[ix]
#         labels = adata_sub.obs[cluster_name]==cl2  # a positive prediction signifies cluster 2
#
#         for celltype in celltypes:
#
#             if SIMPLE:
#                 scores = adata_sub.obs[f'broad_score_{celltype}']
#             else:
#                 from sklearn.ensemble import RandomForestClassifier
#                 cls = RandomForestClassifier(n_estimators=50, oob_score=True)
#                 x = adata_sub.raw[:,marker_dict[celltype]].X
#
#                 if not isinstance(x, np.ndarray):
#                     x = x.A
#
#                 if x.ndim==1:
#                     x = x.reshape(-1,1)
#
#                 cls.fit(x, labels)
#                 scores = cls.oob_decision_function_[:,1]
#
#             a = roc_auc_score(labels, scores)
#             auc_dict[(cl1, cl2, celltype)] = a
#             df.append({'c1': cl1, 'c2': cl2, 'celltype': celltype, 'auc': a})
#             df.append({'c2': cl1, 'c1': cl2, 'celltype': celltype, 'auc': 1-a})  # for symmetry
#
#     df = pd.DataFrame(df)
#     return auc_dict, df



def create_meta_markers(adata, marker_dict, cluster_name='leiden', auc_cutoff=0.7, SIMPLE=True):
    """
    Using a scoring system from Broad (cant find the paper ATM) to give each cell and celltype a score

        marker_dict: str->list[str]: celltype with a list of markers

        auc_cutoff: to annotate a cluster as a specific celltype, it has to have AUC>auc_cutoff
            otherwise we consider it a to ambigous choice
    """

    assert adata.raw, "works only on raw.X data"

    # first, erase all previous marker annotations
    to_drop = [_ for _ in adata.obs.columns if _.startswith('broad_score_')] + ['marker_based_celltype']
    adata.obs.drop(to_drop, axis=1, inplace=True, errors='ignore')

    mol_per_cell = np.array(adata.raw.X.sum(1)).flatten()
    for celltype, markers in marker_dict.items():
        if len(markers) == 1:
            metagene = adata.raw.obs_vector(markers[0])
        else:
#             metagene = np.array(adata.raw[:,markers].X.sum(1)).flatten()
            metagene = np.array(adata.raw[:,markers].X.mean(1)).flatten()

        s = np.log(1e-5 + 10**4 * metagene/mol_per_cell)
        adata.obs[f'broad_score_{celltype}'] = s

    # use the AUC to figure out which celltype each cluster belongs to
    # for each cluster look which cell-score has the best AUC
    df_cluster_celltype = []
    cluster_celltype_mapping = {}
    for cl in tqdm.tqdm(adata.obs[cluster_name].unique()):
        celltypes = marker_dict.keys()

        aucs= {}
        for celltype in celltypes:
            labels = adata.obs[cluster_name]==cl

            if SIMPLE:
                scores = adata.obs[f'broad_score_{celltype}']
            else:
                from sklearn.ensemble import RandomForestClassifier
                cls = RandomForestClassifier(n_estimators=50, oob_score=True)
                x = adata.raw[:,marker_dict[celltype]].X

                if not isinstance(x, np.ndarray):
                    x = x.A

                if x.ndim==1:
                    x = x.reshape(-1,1)

                cls.fit(x, labels)
                scores = cls.oob_decision_function_[:,1]

            a = roc_auc_score(labels, scores)
            aucs[celltype] = a

        argmax_celltype, argmax_auc = dict_argmax_val(aucs)
        if argmax_auc > auc_cutoff:
            # print(argmax_celltype)
            cluster_celltype_mapping[cl] = argmax_celltype
        else:
            print(f'couldnt assign cluster {cl}; best AUC only {argmax_auc}')
            cluster_celltype_mapping[cl] = 'unknown'

        aucs['cluster'] =cl # for the dataframe, also annotate the cluster we're looking at
        df_cluster_celltype.append(aucs)
    df_cluster_celltype = pd.DataFrame(df_cluster_celltype).set_index('cluster').sort_index()

    adata.obs['marker_based_celltype'] = [cluster_celltype_mapping[_] for _ in adata.obs[cluster_name]]

    return df_cluster_celltype, cluster_celltype_mapping




def build_table_literature_vs_differentially_expressed(cluster_celltype_mapping, data_marker_df, literature_marker_dict):
    the_table = []
    for cl in cluster_celltype_mapping.keys():
        putative_celltype = cluster_celltype_mapping[cl]
        literature_markers = ['-'] if putative_celltype == 'unknown' else literature_marker_dict[putative_celltype]
        data_markers = data_marker_df.query('group==@cl').sort_values('pval').head(50).name.values

        data_literature_overlap = len(set(literature_markers) & set(data_markers))

        the_table.append({
            'cluster': cl,
            'putative_celltype': putative_celltype,
            'literature markers': ", ".join(literature_markers),
            'data markers': ', '.join(np.sort(data_markers)),
            'data_literature_overlap': data_literature_overlap
        })

    return pd.DataFrame(the_table).sort_values('putative_celltype')
