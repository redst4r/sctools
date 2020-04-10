import pandas as pd
import os
"""
ROUGH wrapper around PanglaoDB for cell type markers
"""
PANGLAO_GZ = os.path.dirname(__file__) +'/PanglaoDB_markers_21_Oct_2019.tsv.gz'
df_panglao = pd.read_csv(PANGLAO_GZ, sep='\t')

# only markers working in Human (or mouse and human)
df_panglao = df_panglao.query('species == "Mm Hs" or species=="Hs"')
df_panglao = df_panglao[df_panglao['canonical marker'] == 1]

# if sensitivity for human is not present take mouse instead, same for spec
ix = df_panglao['sensitivity_human'] == 0
df_panglao.loc[ix, 'sensitivity'] =  df_panglao.loc[ix, 'sensitivity_mouse']
df_panglao.loc[~ix, 'sensitivity'] =  df_panglao.loc[~ix, 'sensitivity_human']

ix = df_panglao['specificity_human'] == 0
df_panglao.loc[ix, 'specificity'] =  df_panglao.loc[ix, 'specificity_mouse']
df_panglao.loc[~ix, 'specificity'] =  df_panglao.loc[~ix, 'specificity_human']

 # for some weird reason they define spec as the fraction of other cells that express the marker (i.e. 1==all other cells epxress teh marker)
df_panglao['specificity_human'] = 1 - df_panglao['specificity_human']
df_panglao['specificity_mouse'] = 1 - df_panglao['specificity_mouse']
df_panglao['specificity'] = 1 - df_panglao['specificity']


def get_celltype_markers(sens, spec, list_of_organs):
    """
    extract all celltype markers from panglaodb for the given organ and sensitivities, specificities

    sens (float): Sensitivity of the markers
    spec (float): Specificity of the markers
    list_of_organs: all celltypes from these organs are returned
    """
    q = df_panglao.query('sensitivity > @sens and specificity>@spec')
    q = q[q['organ'].isin(list_of_organs)]

    return q.groupby('cell type')['official gene symbol'].apply(lambda x: list(x)).to_dict()
