import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np

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
