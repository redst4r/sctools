from sctools.transforms import split_adata
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import sys
# sys.path.append('/home/michi/ms_python_packages/anndata2ri')
sys.path.append('/home/mstrasse/anndata2ri')
import anndata2ri
import scanpy as sc
import numpy as np
import pandas as pd
import plotnine as pn

def sctransform(adata, n_genes=2000, return_hvg_only=True):
    """
    apply sctransform to the given adata (actually to adata.X)

    transfers the adata to R, and runs the R-version of sctransform
    returns a new adata, with the pearson residuals in .X and some info about the gene-models in .uns
    """
    importr("Seurat")
    importr("sctransform")

    hvg_bool = 'TRUE' if return_hvg_only else 'FALSE'
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv['sce'] = adata
    ro.r('seurat_obj = as.Seurat(sce, counts="X", data = NULL)')
    ro.r(f'res <- SCTransform(object=seurat_obj, assay="originalexp", vst.flavor = "v2",  method = "glmGamPoi", return.only.var.genes = {hvg_bool}, do.correct.umi = TRUE, n_genes={n_genes})')  # , vars.to.regress = "percent_mito",

    with localconverter(ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter):
        # res = ro.r('as.SingleCellExperiment(res)')
        model_fit = ro.r('as.data.frame(res$SCT@SCTModel.list$model1@feature.attributes)')
        pearson_resid = ro.r('as.data.frame(res@assays$SCT@scale.data)').T  # not sparse!
        # corrected_UMIs = ro.r('as.data.frame(res@assays$SCT@counts)').T     # sparse, but returned as DF
        """
        some clarification:
        - `corrected UMIs` will always be full dimension (all genes)
        - `pearson_resid` will vary: if return_hvg_only, this will have less genes
        Hence we cant always store those two in the same adata.

        TODO: pull corrected_umis not as a dataframe, but as a sparse matrix!!
        """

        at = sc.AnnData(pearson_resid.values,
                        obs=adata.obs.loc[pearson_resid.index],
                        var=model_fit.loc[pearson_resid.columns.values])
        at.uns['sctransform_model'] = model_fit

        """
        some explanation on the model_fit DataFrame:
        We fit a NB-model where the linear part is:
        log(mu) = beta0 + beta1 * log(n_umi)

        - log_umi: this is the COEFFICIENT beta1 (confusing naming!!). This is actually fixed to log(10)
        - (Intercept): beta0
        """

    return at


def plot_variance(adata):
    """
    plot the mean variance relation and the residual variance (after sctransform)
    """
    df_sctransform = adata.uns['sctransform_model'].copy()
    hvg_genes = set(adata.var.index)
    df_sctransform['is_HVG'] = df_sctransform.index.map(lambda x: 'yes' if x in hvg_genes else 'no')

    points =  pn.geom_point(alpha=0.2, size=0.2)
    p1 = pn.ggplot(df_sctransform, pn.aes('gmean', 'variance', color='is_HVG')) + points + pn.geom_abline(slope=1, color='red') + pn.scale_x_log10()+ pn.scale_y_log10()
    p2 = pn.ggplot(df_sctransform, pn.aes('gmean', 'residual_variance', )) + points + pn.scale_x_log10()+ pn.scale_y_log10()  + pn.geom_density_2d(color='red')

    return p1, p2


import plydata
def plot_model(adata):
    """
    plot the fitted model parameters
    """
    df_sctransform = adata.uns['sctransform_model'].copy()
    # mark variable genes
    hvg_genes = set(adata.var.index)
    df_sctransform['is_HVG'] = df_sctransform.index.map(lambda x: 'yes' if x in hvg_genes else 'no')
    # define OD
    df_sctransform = df_sctransform >> plydata.define(step1_OD='1+ gmean/step1_theta', OD='1+ gmean/theta')

    theme = pn.theme(figure_size=(3,3))
    # Mean vs Intercept
    p1 = pn.ggplot(df_sctransform, pn.aes('gmean', 'step1_(Intercept)', color='is_HVG')) + pn.geom_point(size=0.5) + pn.scale_x_log10() + pn.geom_point(pn.aes('gmean', y='(Intercept)'), size=0.1, color='black') + theme + pn.labs(title='Mean vs Intercept')

    # Mean vs Theta
    p2 = pn.ggplot(df_sctransform, pn.aes('gmean', 'step1_theta', color='is_HVG')) + pn.geom_point(size=0.5) + pn.scale_x_log10() + pn.scale_y_log10() + pn.geom_point(pn.aes('gmean', y='theta'), size=0.1, color='black') + theme+ pn.labs(title='Mean vs Theta')

    # Mean vs Overdispersion (1+mean/theta)
    p3 = pn.ggplot(df_sctransform, pn.aes('gmean', 'step1_OD', color='is_HVG')) + pn.geom_point(size=0.5, alpha=0.5) + pn.scale_x_log10() + pn.scale_y_log10() + pn.geom_point(pn.aes('gmean', y='OD'), size=0.1, color='black') + theme+ pn.labs(title='Mean vs Overdispersion')

    return p1, p2, p3


def predict_model(sct_df, gene, n_molecules):
    """
    predict the gene expression value under the SCT model for a given gene and a vector of sequencing depth per cell (n_molecules)

    :param sct_df: The model output of sctransform. Usually stored in adata.uns['sctransform_model']
    """

    # the model fit/prediction of the expression values (null model)
    gene_parameters = sct_df.loc[gene]
    # the response of the linear model: log_mu = intercept + b1 * log10(n_molecules)
    log_mu = gene_parameters['(Intercept)'] + gene_parameters['log_umi'] * np.log10(n_molecules)
    mu = np.exp(log_mu)
    theta = sct_df.loc[gene]['theta']
    var = mu + (mu**2 / theta)
    std = np.sqrt(var)

    model_data = pd.DataFrame({
        'n_molecules': n_molecules,
        'mu': mu,
        'theta': theta,
        'std': std
    }).sort_values('n_molecules')
    return model_data

def plot_model_fit(adata_sc, adata_raw, gene, color=None):
    """
    plot the n_molecules vs expression relation for the actual data and the model prediction
    """
    sct_df = adata_sc.uns['sctransform_model']

    # the dataframe for the raw/unnormed count data
    count_data = adata_raw.obs.copy()
    count_data['gene_expression'] = adata_raw[:, gene].X.A.flatten()
#     count_data['pearson_residual'] = adata_sc[:, gene].X.flatten()

    model_data = predict_model(sct_df, gene, adata_raw.obs.n_molecules.values)

    # the model: expected epxresion +/- std
    pn_model = pn.ggplot(model_data, pn.aes('n_molecules', '1+mu')) + pn.geom_line(size=2) + pn.geom_line(pn.aes('n_molecules', '1+mu+std'), color='grey', size=2) + pn.geom_line(pn.aes('n_molecules', '1+mu-std'),size=2, color='grey')

    # raw data plot
#     color = 'pearson_residual' if color is None else color
    if color is None:
        pn_rawdata =  pn.geom_point(pn.aes('n_molecules','1+gene_expression'), data=count_data,size=0.1, alpha=0.25)
        pn_rawdata_density = pn.geom_density_2d(pn.aes('n_molecules','1+gene_expression'), data=count_data, color='blue', levels=10)
        p = pn_model + pn_rawdata +pn_rawdata_density

    else:
        pn_rawdata =  pn.geom_point(pn.aes('n_molecules','1+gene_expression', color=color), data=count_data,size=0.1)
        p = pn_model + pn_rawdata

    return p + pn.scale_x_log10() + pn.scale_y_log10()
