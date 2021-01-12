import gzip

"""
a few functions handy for dealing with the more raw 10x data
"""


def _load_whitelist(fname):
    "loads a whitelist of cellbarcdes as a set"
    with open(fname, 'r') as fh:
        whitelist = set(_.strip() for _ in fh.readlines())
    return whitelist


def load_v3_CB_feature_translation(fname):
    """
    v3 uses different CB sequences for the expression library and the feature-library:
    even thugh its the same cell/bead, it has two different CB sequences!

    the mapping is stored in a file from cellranger: `cellranger-x.y.z/cellranger-cs/x.y.z/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz`
        - the first column is the expression CB
        - the second column is the feature CB

    see also https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-

   return:
       dict, mapping featureCB -> expressionCB


    load_v3_CB_feature_translation('/home/michi/mounts/TB4drive/kallisto_resources/translation_3M-february-2018.txt.gz')
   """
    the_dict = {}
    with gzip.open(fname, 'rt') as fh:
        for line in fh:
            exprCB, featureCB = line.strip().split()
            the_dict[featureCB] = exprCB
    return the_dict
