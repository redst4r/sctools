# README

sctools is a set of convenience functions for dealing with scRNAseq data within
scanpy.


## Pipeline
Functions to do fairly standard processing of scRNAseq data contained in scanpys AnnData objects.
These functions essentially turn raw count data into normalized, clustered and dimension reduced data.
```python
from sctools.pipeline import standard_processing
adata = standard_processing(adata)

sc.pl.umap(adata, color=['leiden', 'doublet_score'])

# differnetially expressed genes between clusters
sc.pl.rank_genes_groups(adata)
```

## Differntial expression
`scanpy_DE_to_dataframe_fast()` turns scanpys awkward way of storing differential expression
into a bunch of `pandas.DataFrame`s  
