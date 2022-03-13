import numpy as np
import matplotlib.pyplot as plt


def biplot(adata, components):
    """
    PCA biplot, annotating the directions of the genes in PCA space

    :param components: tuple. which PCs to plot. Indexing starts at 1, i.e. to get the first and third component use (1,3)
    """
    assert len(components) == 2 and components[0]>0 and components[1]>0
    components = [_-1 for _ in components]
    pcs = adata.varm['PCs'][:,components]
    # scale them
    pcs = pcs*3
    component_names = adata.var.index.values

    # take the ones with the biggest direction per dimesnion
    ix1 = np.argsort(pcs[:,components[0]])[::-1]
    ix2 = np.argsort(pcs[:,components[1]])[::-1]
    ix1_neg = np.argsort(pcs[:,components[0]])
    ix2_neg = np.argsort(pcs[:,components[1]])

    # take the ones with the biggest magnitude
    m = np.sqrt(np.sum(pcs**2, axis=1))
    ix3 = np.argsort(m)[::-1]


    ix = np.concatenate([ix3[:10], ix2[:5], ix1[:5], ix1_neg[:5], ix2_neg[:5]])

    pcs = pcs[ix]
    component_names = component_names[ix]

    xs = adata.obsm['X_pca'][:, components[0]]
    ys = adata.obsm['X_pca'][:, components[1]]

#     plt.scatter(xs, ys, s=1, alpha=0.5)

    for i in range(pcs.shape[0]):
    # arrows project features (ie columns from csv) as vectors onto PC axes
        plt.arrow(0, 0, pcs[i,0]*max(xs), pcs[i,1]*max(ys),
                  color='k', width=0.0005, head_width=0.0025)
        plt.text(pcs[i,0]*max(xs)*1.2, pcs[i,1]*max(ys)*1.2,
                 component_names[i], color='k', fontdict={'fontsize': 8})
    print(component_names)
