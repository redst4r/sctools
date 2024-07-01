import numpy as np
from sklearn.metrics import roc_auc_score
import pandas as pd
import tqdm
import toolz
import sctools.celltypes.panglaodb as panglaodb


def build_marker_dict(adata, tissue):
    """
    my custom list of markers, partially based on panglaodb, partially on personal communication
    """

    assert tissue in ['Gastric', 'Esophagus', 'Lung'], f"Tissue {tissue} unknown. Should be Gastric, Esophagus or Lung"

    """
    Gastric and esophagus have overlapping celltypes, so treat them as one.
    Especially due to the metaplasia, where Esophagus starts to look like Gastric/Intestinal

    Lung celtypes are different however!
    """
    organ_dict = {
        'Gastric': ['Connective tissue', 'Immune system', 'GI tract', 'Blood', 'Vasculature', 'Smooth muscle', 'Epithelium'],
        'Esophagus': ['Connective tissue', 'Immune system', 'GI tract', 'Blood', 'Vasculature', 'Smooth muscle', 'Epithelium'],
        'Lung': ['Connective tissue', 'Immune system', 'Lungs', 'Blood', 'Vasculature', 'Smooth muscle', 'Epithelium']
    }

    celltype_marker_dict = panglaodb.get_celltype_markers(
        sens=0.8, spec=0.8,
        list_of_organs=organ_dict[tissue]
        )

    # now fix a few by hand and expert knowledge
    if tissue in ['Gastric', 'Esophagus']:
        celltype_marker_dict['Gastric chief cells'] = ['PGA4', 'PGA3', 'PGA5', 'PGC']
        celltype_marker_dict['Enteroendocrine cells'].extend(['CHGA', 'CHGB'])
        celltype_marker_dict['Parietal cells'] = ['ATP4A', 'ATP4B']
        celltype_marker_dict['Endothelial cells'].append('PECAM1')
        celltype_marker_dict['Goblet cells'] = ['MUC2', 'ITLN1', 'SPINK4']
        celltype_marker_dict['Intestinal epithelial cells'] = ['CDX2']
        celltype_marker_dict['Foveolar cells 1'] = ['MUC5AC', 'TFF1' ]
        celltype_marker_dict['Foveolar cells 2'] = ['TFF1', 'GKN1']
        celltype_marker_dict['Mucus neck cells'] = ['TFF1', 'MUC6']
        celltype_marker_dict['MSC'] = ['EPHB2', 'SOX9', 'OLFM4', 'CD44', 'AXIN2', 'SOX4']
        celltype_marker_dict['Paneth cells'] = ['DEFA5', 'DEFA6']
    elif tissue == 'Lung':
        # these are from the Teichman/HCA paper
        celltype_marker_dict['Basal epithelial'] = ['KRT5', 'TP63']
        celltype_marker_dict['Ciliated'] = ['FOXJ1', 'KCTD12', 'PIFO']
        celltype_marker_dict['Secretory/goblet_1'] = ['MUC5AC', 'KRT4', 'CD36']
        celltype_marker_dict['Secretory/goblet_2'] = ['MUC5AC', 'CXCL10', 'IDO1', 'NOS2', 'IL19']
        celltype_marker_dict['Aviolar-T2'] = ['SFTPC']
        celltype_marker_dict['Aviolar-T2'] = ['AGER']
        celltype_marker_dict['Club/Clara'] = ['MUC5AC', 'MUC5B']
    else:
        raise ValueError(f'unknown tissue: {tissue}')

    # delete a few that we dont see anyways
    del celltype_marker_dict['Platelets']
    del celltype_marker_dict['Neutrophils']
#     del celltype_marker_dict['Pericytes'] # <-?? maybe should be included
    del celltype_marker_dict['Adipocytes']
    del celltype_marker_dict['Eosinophils']

    # filter out all genes that we dont measure in adata
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
    for k, v in d.items():
        if v > max_val:
            max_key = k
            max_val = v
    return max_key, max_val


def annotate_celltypes_broad(adata, marker_dict, cluster_name='leiden', auc_cutoff=0.7, SIMPLE=True, CELLTYPE_COLNAME = 'marker_based_celltype'):
    """
    Using a scoring system from Broad (cant find the paper ATM) to give each cell and celltype a score.
    To annotate a cluster of cells with a cell type, the AUC is used:
    How well can a certain set of markers distinguish the cluster from all other cells.

    Parameters:
        marker_dict: str->list[str]: celltype with a list of markers

        auc_cutoff: to annotate a cluster as a specific celltype, it has to have AUC>auc_cutoff
            otherwise we consider it a to ambigous choice
    Return:
        - modifies the adata, adding the scores and final celltype as columns into adata.obs
        - df_cluster_celltype: each cluster/celltye AUC as a pd.DataFrame
        - cluster_celltype_mapping (dict): mapping cluster --> cell type
    """

    assert adata.raw, "works only on raw.X data"

    PREFIX = 'broad_score_'

    # first, erase all previous marker annotations
    to_drop = [_ for _ in adata.obs.columns if _.startswith(PREFIX)] + [CELLTYPE_COLNAME]
    adata.obs.drop(to_drop, axis=1, inplace=True, errors='ignore')

    # each list of markers (for one celltype) is turned into a single score,
    # which is just the average of all markers in that cell
    # it gets normalized by cellsize and log-transformed
    mol_per_cell = np.array(adata.raw.X.sum(1)).flatten()
    for celltype, markers in marker_dict.items():
        if len(markers) == 1:
            metagene = adata.raw.obs_vector(markers[0])
        else:
            metagene = np.array(adata.raw[:, markers].X.mean(1)).flatten()

        s = np.log(1e-5 + 10**4 * metagene/mol_per_cell)
        adata.obs[f'{PREFIX}{celltype}'] = s

    # use the AUC to figure out which celltype each **cluster** belongs to
    # for each cluster look which cell-score has the best AUC
    df_cluster_celltype = []
    cluster_celltype_mapping = {}
    for cl in tqdm.tqdm(adata.obs[cluster_name].unique()):
        celltypes = marker_dict.keys()

        aucs = {}
        for celltype in celltypes:
            labels = adata.obs[cluster_name] == cl

            if SIMPLE:
                scores = adata.obs[f'{PREFIX}{celltype}']
            else:
                from sklearn.ensemble import RandomForestClassifier
                cls = RandomForestClassifier(n_estimators=50, oob_score=True)
                x = adata.raw[:, marker_dict[celltype]].X

                if not isinstance(x, np.ndarray):
                    x = x.A

                if x.ndim == 1:
                    x = x.reshape(-1, 1)

                cls.fit(x, labels)
                scores = cls.oob_decision_function_[:, 1]

            a = roc_auc_score(labels, scores)
            aucs[celltype] = a

        argmax_celltype, argmax_auc = dict_argmax_val(aucs)
        if argmax_auc > auc_cutoff:
            # print(argmax_celltype)
            cluster_celltype_mapping[cl] = argmax_celltype
        else:
            print(f'couldnt assign cluster {cl}; best AUC only {argmax_auc}')
            cluster_celltype_mapping[cl] = 'unknown'

        aucs['cluster'] = cl  # for the dataframe, also annotate the cluster we're looking at
        df_cluster_celltype.append(aucs)
    df_cluster_celltype = pd.DataFrame(df_cluster_celltype).set_index('cluster').sort_index()

    adata.obs[CELLTYPE_COLNAME] = [cluster_celltype_mapping[_] for _ in adata.obs[cluster_name]]

    return df_cluster_celltype, cluster_celltype_mapping

# some legacy naming:
create_meta_markers = annotate_celltypes_broad


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


import celltypist
from celltypist import models
from sctools import pipeline
import scanpy as sc
def do_celltypist(adata, modelname='Immune_All_Low.pkl'):
    """
    needs to have raw data in adata.raw
    """

    adata_for_celltypist = pipeline.export_for_cellxgene(adata.copy(), adata.obs.columns)
    sc.pp.normalize_total(adata_for_celltypist, target_sum=10000)
    sc.pp.log1p(adata_for_celltypist)
    predictions = celltypist.annotate(
        adata_for_celltypist, 
        model = modelname, 
        # model = 'Immune_All_High.pkl', 
        majority_voting = True,
        # for multilabel prediction (also "Unassigned")
        # mode = 'prob match', p_thres = 0.5
    )
    # A = predictions.to_adata()
    # just pick out the relevant attributes
    df_celltypist = predictions.to_adata().obs[['predicted_labels', 'over_clustering', 'majority_voting','conf_score']]
    return df_celltypist
