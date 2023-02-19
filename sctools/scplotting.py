import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import tqdm
import plotnine as pn

godsnot_64 = [
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#ff000d", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C"]

COLORS = godsnot_64

expression_cmap = LinearSegmentedColormap.from_list('mycmap', [(0.7,0.7,0.7),plt.cm.Reds(0.3), plt.cm.Reds(0.9)])
godsnot_cmap = ListedColormap(godsnot_64, 'godsnot_64')

def rand_cmap(nlabels, type='bright', first_color_black=True, last_color_black=False, verbose=False):
    """
    Creates a random colormap to be used together with matplotlib. Useful for segmentation tasks
    :param nlabels: Number of labels (size of colormap)
    :param type: 'bright' for strong colors, 'soft' for pastel colors
    :param first_color_black: Option to use first color as black, True or False
    :param last_color_black: Option to use last color as black, True or False
    :param verbose: Prints the number of labels and shows the colormap. True or False
    :return: colormap for matplotlib
    """
    from matplotlib.colors import LinearSegmentedColormap
    import colorsys
    import numpy as np


    if type not in ('bright', 'soft'):
        print ('Please choose "bright" or "soft" for type')
        return

    if verbose:
        print('Number of labels: ' + str(nlabels))

    # Generate color map for bright colors, based on hsv
    if type == 'bright':
        randHSVcolors = [(np.random.uniform(low=0.0, high=1),
                          np.random.uniform(low=0.2, high=1),
                          np.random.uniform(low=0.9, high=1)) for i in range(nlabels)]

        # Convert HSV list to RGB
        randRGBcolors = []
        for HSVcolor in randHSVcolors:
            randRGBcolors.append(colorsys.hsv_to_rgb(HSVcolor[0], HSVcolor[1], HSVcolor[2]))

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]

        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Generate soft pastel colors, by limiting the RGB spectrum
    if type == 'soft':
        low = 0.6
        high = 0.95
        randRGBcolors = [(np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high),
                          np.random.uniform(low=low, high=high)) for i in range(nlabels)]

        if first_color_black:
            randRGBcolors[0] = [0, 0, 0]

        if last_color_black:
            randRGBcolors[-1] = [0, 0, 0]
        random_colormap = LinearSegmentedColormap.from_list('new_map', randRGBcolors, N=nlabels)

    # Display colorbar
    if verbose:
        from matplotlib import colors, colorbar
        from matplotlib import pyplot as plt
        fig, ax = plt.subplots(1, 1, figsize=(15, 0.5))

        bounds = np.linspace(0, nlabels, nlabels + 1)
        norm = colors.BoundaryNorm(bounds, nlabels)

        cb = colorbar.ColorbarBase(ax, cmap=random_colormap, norm=norm, spacing='proportional', ticks=None,
                                   boundaries=bounds, format='%1i', orientation=u'horizontal')

    return random_colormap

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb
def recolor(adata, field, colors):
    """
    recolor the column `adata.obs[field]` using the specified color_vecotr

    handy as scnapy sometimes screws up the coloring for columns with MANY categories
    """
    colorfield = f'{field}_colors'
    if colorfield in adata.uns:
        del adata.uns[colorfield]

    if isinstance(colors, dict):
        cats = adata.obs[field].cat.categories.values
        cvector = [colors[i] for i in cats]
    else:
        cats = adata.obs[field].unique().astype('str')
        if len(cats) > len(colors):
            colors = colors + colors
        cvector = [colors[i] for i in range(len(cats))]

    assert all(c.startswith('#') and (len(c) == 7 or len(c) == 9 )for c in cvector), "colors invalid"
    adata.uns[colorfield] = cvector


def kneeplot(adata, expected_num_cells):
    "also return the UMI threshold for true cells based on expt number"
    knee = np.sort((np.array(adata.X.sum(axis=1))).flatten())[::-1]
    ax =plt.gca()
    # fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(range(len(knee)), knee, label="kallisto", linewidth=5, color="k")
    ax.axvline(x=expected_num_cells, linewidth=3, color="g")
    ax.axhline(y=knee[expected_num_cells], linewidth=3, color="g")

    ax.set_xlabel("Set of Barcodes")
    ax.set_ylabel("UMI Counts")

    plt.grid(True, which="both")
    ax.legend()
    # plt.show()

    return knee[expected_num_cells]

def kneeplot_split_plotnine(adata, splitfield='samplename'):

    if isinstance(splitfield, str):
        splitfield = [splitfield]

    df_result = []
    for group, obs_group in adata.obs.groupby(splitfield):
        O = obs_group.sort_values('n_molecules', ascending=False)[['n_molecules']]
        O['x'] = np.arange(O.shape[0])
        O['norm_n_molecules'] = O['n_molecules'] / O['n_molecules'].sum()

        for s in splitfield:
            O[s] = obs_group[s].values
        df_result.append(O)

    return pd.concat(df_result)


def kneeplot_split(adata, splitfield='samplename'):
    splits = sorted(adata.obs[splitfield].unique())
    for s in splits:
        a_tmp = adata[adata.obs[splitfield]==s]
        "also return the UMI threshold for true cells based on expt number"
        knee = np.sort((np.array(a_tmp.X.sum(axis=1))).flatten())[::-1]
        ax =plt.gca()
        ax.loglog(range(len(knee)), knee, label=s, linewidth=5, )

        ax.set_xlabel("Set of Barcodes")
        ax.set_ylabel("UMI Counts")

        plt.grid(True, which="both")
        ax.legend()

def kneeplot_split_cumulative(adata, splitfield='samplename'):
    splits = sorted(adata.obs[splitfield].unique())
    for s in splits:
        a_tmp = adata[adata.obs[splitfield]==s]
        cumulative_knee = np.cumsum(np.sort((np.array(a_tmp.X.sum(axis=1))).flatten())[::-1])
        ax =plt.gca()
#         ax.loglog(range(len(cumulative_knee)), cumulative_knee, label=s, linewidth=5, )
        ax.plot(range(len(cumulative_knee)), cumulative_knee, label=s, linewidth=5, marker='x' )

        ax.set_xlabel("Set of Barcodes")
        ax.set_ylabel("cumulative UMI Counts")

        plt.grid(True, which="both")
        ax.legend()


from sctools.transforms import split_adata
from scipy.stats import poisson
import pandas as pd

def sparse_var(X, axis):
    """
    calc variance of sparse matrix
    """
    m = X.mean(axis).A1
    _t = X.copy()
    _t.data = _t.data**2
    return _t.mean(axis).A1 - m**2

def mean_var_relation(adata):
    mean_exp = adata.X.mean(0).A1
    var_exp = sparse_var(adata.X, axis=0)
    detection_rate = (adata.X>0).mean(0).A1

    mean_var_df = pd.DataFrame({
        'mean_exp': mean_exp,
        'var_exp': var_exp,
        'detection_rate': detection_rate,
        'gene': adata.var.index.values
    })
    mean_var_df['poisson_detection_rate'] = 1 - poisson.pmf(0, mu=mean_var_df.mean_exp)
    return mean_var_df


def mean_var_relation_split(adata, split):
    mean_var_df  = []
    for sname, a in tqdm.tqdm(split_adata(adata, split_field=split)):
        _d = mean_var_relation(a)
        _d[split] = sname
        mean_var_df.append(_d)

    mean_var_df = pd.concat(mean_var_df)
    return mean_var_df

def plot_mean_var(adata, split=None):
    if split is None:
        df = mean_var_relation(adata)
    else:
        df = mean_var_relation_split(adata, split)

    p1 = pn.ggplot(df, pn.aes('mean_exp', "var_exp")) + pn.geom_point(size=0.1, alpha=0.1) + pn.geom_abline(slope=1, color='red') +pn.scale_x_log10() + pn.scale_y_log10()
    p2 = pn.ggplot(df, pn.aes('mean_exp', "detection_rate")) + pn.geom_point(size=0.1, alpha=0.1) +pn.scale_x_log10() + pn.geom_line(pn.aes('mean_exp', "poisson_detection_rate"), color='red')
    if split is not None:
        p1 = p1 + pn.facet_wrap('samplename')
        p2 = p2 + pn.facet_wrap('samplename')
    return p1, p2


def extract_cmap_from_adata(adata, field: str):
    """
    turn the colors stored in adata.uns into a dict: category->color
    :param adata:
    :param field: which field in .obs to get the color coding for
    """
    assert adata.obs[field].dtype == "category"
    cm = {}
    for i,cat_name in enumerate(adata.obs[field].cat.categories):
        cm[cat_name] = adata.uns[f'{field}_colors'][i]
    return cm

def extract_colorvector_from_adata(adata, field: str):
    """
    for a categorical cell attribute in .obs, extract a vector that encodes
    each cell by the respecitve attribute of the cell (colors stored in adata.uns
    :param adata:
    :param field: which field to get the color vector for
    """
    cmap = extract_cmap_from_adata(adata, field)

    cvector = [cmap[_] for _ in adata.obs[field]]
    return cvector
