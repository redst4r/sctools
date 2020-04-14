"""
a chaotic collection of genesets, mostly for cell type markers
"""

# from teichman HCA lung paper
lung_markers = {
    'basal': ['KRT5', 'TP63'],
    'ciliated': ['FOXJ1', 'KCTD12', 'PIFO'],
    'secretory/goblet': ['CEACAM5', 'S100A4', 'MUC5AC'],
    'secretory/goblet_1': ['MUC5AC', 'KRT4', 'CD36'],
    'secretory/goblet_2': ['MUC5AC', 'CXCL10', 'IDO1', 'NOS2', 'IL19'],
    'aviolar-T2': ['SFTPC'],
    'aviolar-T1': ['AGER'],
    'club/clara': ['MUC5AC', 'MUC5B'],

}


other_signatures = {
    # from https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html
    'dissociation': ["ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA", "NR4A1", "PPP1R15A", "SOCS3", "ZFP36"]
}

# KRT19: surface mucus
stemness_genes = ['POU5F1','SOX2','MYC','KLF4','NANOG','ZFP42','PROM1','CTNNB1','ABCG2','CD34','CD44','ZSCAN4','LGR5','BMI1','BMI1','BMI1','EPAS1','EZH2','NES','TWIST1','HIF1A','NOTCH1','CTNNB1','KDM5B']

nuclear_mito_genes = [
    # complex 1
    'NDUFS7', 'NDUFS8', 'NDUFV2', 'NDUFS3', 'NDUFS2', 'NDUFV1', 'NDUFS1', 'NDUFS6', 'NDUFA12', 'NDUFS4', 'NDUFA9', 'NDUFAB1', 'NDUFA2', 'NDUFA1', 'NDUFB3', 'NDUFA5', 'NDUFA6', 'NDUFA11', 'NDUFB11', 'NDUFS5', 'NDUFB4', 'NDUFA13', 'NDUFB7', 'NDUFA8', 'NDUFB9', 'NDUFB10', 'NDUFB8', 'NDUFC2', 'NDUFB2', 'NDUFA7', 'NDUFA3', 'NDUFA4', 'NDUFB5', 'NDUFB1', 'NDUFC1', 'NDUFA10', 'NDUFA4L2', 'NDUFV3', 'NDUFB6', 'NDUFAF1', 'NDUFAF2', 'NDUFAF3', 'NDUFAF4',
    # complex2
    'SDHA', 'SDHB', 'SDHC', 'SDHD',

    # complex 3
    'ETFDH',

    #complex 4
    'CYC1', 'UQCRFS1', 'UQCRC1', 'UQCRC2', 'CYC1', 'UQCRB', 'UQCRH', 'UQCRQ', 'UQCR10', 'UQCR11',

    #complex 5
    'COX4I1', 'COX4I2', 'COX5A', 'COX5B', 'COX6A1', 'COX6A2', 'COX6B1', 'COX6B2', 'COX6C', 'COX7A1', 'COX7A2', #'COX7A3', 'COX7B', 'COX7C', 'COX7A2L', 'COX8A', 'COX8C', 'COA1', 'COA3', 'COA4', 'COA5', 'COA6', 'COA7', 'COX11', 'COX14', 'COX15', 'COX16', 'COX17', 'COX18', 'COX19', 'COX20'
]


keratins = ['CASP14', 'CSNK1A1', 'EPPK1', 'FAM83H', 'FBF1', 'GPER1', 'KRT1', 'KRT13', 'KRT14', 'KRT18', 'KRT2', 'KRT3', 'KRT4', 'KRT5', 'KRT6A', 'KRT6B', 'KRT6C', 'KRT7', 'KRT71', 'KRT72', 'KRT73', 'KRT74', 'KRT75', 'KRT76', 'KRT77', 'KRT78', 'KRT79', 'KRT8', 'KRT80', 'KRT81', 'KRT82', 'KRT83', 'KRT84', 'KRT85', 'KRT86', 'KRT87P', 'KRTAP1-1', 'KRTAP1-3', 'KRTAP1-4', 'KRTAP1-5', 'KRTAP10-1', 'KRTAP10-10', 'KRTAP10-11', 'KRTAP10-12', 'KRTAP10-2', 'KRTAP10-3', 'KRTAP10-4', 'KRTAP10-5', 'KRTAP10-6', 'KRTAP10-7', 'KRTAP10-8', 'KRTAP10-9', 'KRTAP11-1', 'KRTAP12-1', 'KRTAP12-2', 'KRTAP12-3', 'KRTAP12-4', 'KRTAP16-1', 'KRTAP2-1', 'KRTAP2-3', 'KRTAP2-4', 'KRTAP24-1', 'KRTAP29-1', 'KRTAP3-1', 'KRTAP3-2', 'KRTAP3-3', 'KRTAP4-1', 'KRTAP4-11', 'KRTAP4-12', 'KRTAP4-2', 'KRTAP4-3', 'KRTAP4-4', 'KRTAP4-5', 'KRTAP4-6', 'KRTAP4-8', 'KRTAP4-9', 'KRTAP5-1', 'KRTAP5-10', 'KRTAP5-11', 'KRTAP5-2', 'KRTAP5-3', 'KRTAP5-4', 'KRTAP5-5', 'KRTAP5-6', 'KRTAP5-7', 'KRTAP5-8', 'KRTAP5-9', 'KRTAP9-1', 'KRTAP9-2', 'KRTAP9-3', 'KRTAP9-4', 'KRTAP9-6', 'KRTAP9-7', 'KRTAP9-8', 'KRTAP9-9', 'TCHP']

Stuart_gastric_markers = {
    'goblet': ['MUC2'],
    'Foveolar cells': ['MUC5AC'],  # PMC
    'Parietal cells': ['ATP4A', 'ATP4B'],
    'chief': ['PGA4', 'PGA3', 'PGA5', 'PGC'],
    'Paneth cells': ['DEFA5',],
    'Mucous secreting cells': ['MUC6'],
    'Intestinal epithelial cells': ['CDX2']
}

dd = {
    'T-cell': ['CD2','CD3D','CD3E','CD3G'],
    'B-cell': ['CD19','CD79A','CD79B','BLK'],
    'Macrophage': ['CD163','CD14','CSF1R'],
    'Endothelial': ['PECAM1',   'VWF',   'CDH5'],
    'CAFs': ['FAP', 'THY1', 'DCN', 'COL1A1', 'COL1A2', 'COL6A1', 'COL6A2', 'COL6A3'],

# Table S12. Genes preferentially expressed by Tregs compared to CD4+ and CD8+ T-cells
    'Treg': ['IL2RA',    'FOXP3',    'S100A4',    'CCR8',    'TNFRSF1B',    'GBP5']
}

my_markers = {
    'enteroendocrine': ['GAST'],
    'chief': ['BHLHA15'],
    'proliferative': ['TOP2A'],        #  proliferative
    'goblet': ['SPINK4'],        # goblet
    'enterocytes': ['ALPL'],        # enterocytes
    'intestinal Stem cell': ['SOX4', 'AXIN2', 'CD44'],  # intestinal Stem cell markers
}

# from Single-cell transcriptomics of human t cells revelas tissue and activation signatarues in health and disease
tcell_markers = {
    'CD4_act': ['IL2', 'TNF', 'IL4R'],
    'CD4_Treg': ['FOXP3', 'IL2RA', 'CTLA4'],
    'CD4_TCM': ['CCR7', 'SELL', 'TCF7'],  # circulating memory T "resting"
    'CD4_TRM': ['CXCR6', 'ITGA1'],  # tissue resident memory T "resting"
    'CD8_TEM_TRM': ['CCL5','GZMB','GZMK','CXCR6', 'ITGA1'],
    'CD8_TEM_TRM_act': ['IFNG','CCL4', 'CCL3'],
    'CD8_TEMRA': ['PRF1', 'NKG7']  #terminally differentiated effector cells

    }

Zhang2019_gastric_markers = {
    # epithelial
    'epithelial': ['MUC1', 'KRT18', 'EPCAM'],
    # non-epithelial
    'non-epithelial': ['VIM', 'PTPRC'],
    ##############################
    ## epithelial subgroups:
    ##############################
        'GMC' :['MUC6', 'TFF2'],        # GMC
        'PMC/Foveolar':['MUC5AC', 'TFF1'],        # PMC
        'chief': ['PGA4', 'PGA3'],        # Chief cells
        'enteroendocrine': ['CHGA', 'CHGB'],          # enteroendocrine
        'proliferative': ['MKI67', 'BRIC5'],        #  proliferative
        'cancer cells': ['CEACAM5', 'CEACAM6', 'BAX', 'CCDND2',],    # cancer cells
        'goblet': ['ITLN1', 'MUC2'],        # goblet
        'enterocytes': ['FABP1', 'APOA1', 'APOA4'],        # enterocytes
        'intestinal Stem cell': ['OLFM4', 'EPHB2', 'SOX9'],  # intestinal Stem cell markers
        'Parietal cells': ['ATP4A', 'ATP4B'],
        'Paneth cells': ['DEFA5', 'DEFA6'],
        'Intestinal epithelial cells': ['CDX2'],
    ##############################
    ## NON EPHITHELIAL SUBTYPES
    ##############################
        'T-cell': ['CD2', 'CD3D'],        # t-cell
        'B-cell': ['CD79A', 'MS4A1'],        # b cell
        'Macrophage': ['CSF1R', 'CD68'],        # macrophahe
        'Fibro': ['DCN', 'PDPN'],        # fibroblasts
        'SMC': ['ACTA2'],        # smooth muscle
        'Endothelial': ['VWF', 'ENG'],        # endothelial
        'Mast cell': ['TPSAB1', 'TPSB2'],        #mast cells
        'RBC': ['HBA1', 'HBA2', 'HBB']

    #######
    # H.pylori infection
    #'REG3A', 'LCN2', 'COX7B', 'UQCRB'
}

# from the scanpy website
blood_markers = {
                     'B-cell': ['CD79A', 'MS4A1'],
                     'T-cell': ['CD3D'],
                     'T-cell CD8+': ['CD8A', 'CD8B'],
                     'NK': ['GNLY', 'NKG7', 'FCGR3A', 'FCGR3B','NCAM1'],  # FCGR3A CD16A
                     'Myeloid': ['CST3', 'LYZ'],
                     'Monocytes': ['FCGR3A'],
                     'Dendritic': ['FCER1A']}


# CD markers
CD_markers = {
    'T-cell': ['CD3D','CD4','CD8B', 'CD8A'],
    'B-cell': ['CD19', 'MS4A1'],  # MS4A1 = CD20
    'Dendritic': ['ITGAX', 'IL3RA'], # CD11c = ITGAX , IL3RA=CD123
    'NK': ['NCAM1'],  # NCAM1 = CD56
    'HSC': ['CD34'],
    'Macrophage': ['CD14','CD33'],
    #'Granulocyte': ['CEACAM8'],  # CD66b = CEACAM8
    'Platelet': ['ITGA2B', 'ITGB3','SELP'],  # CD41  =ITGA2B;  CD61=ITGB3, SELP=Cd62
    'Erythrocyte': ['GYPA'],  #'GYPA' = CD235a
    'Blood-Endothelial': ['MCAM'],  # CD146 = MCAM
    'Blood-Epithelial' : ['GYPC']  # CD236= GYPC
}


def filter_dict_for_measured(s, Q):
    new_d = {}
    for k,v in s.items():
        assert isinstance(v, list)
        new_d[k]= [_ for _ in v if _ in Q.var.index]
    return new_d

def merge_markers(dict1, dict2):
    new_d = dict1.copy()

    for k, markerlist in dict2.items():
        if k in new_d:
            # see if we can extend d1s list
            new_d[k] = list(set(new_d[k] + markerlist))
        else:
            new_d[k] = markerlist
    return new_d




"""
MArkers from the Lung,spleen and oesophagus tissue remains stable ...
"""

HCA_markers = {
    'Epi_basal': ['KRT5', 'SCPEP1', 'ITGA6'],
    'Epi_stratified': ['KRT4'],
    'Epi_upper': ['ECM1'],
    'Glands_duct': ['KRT7', 'KRT23'],
    'Glands_mucus': ['TFF3', 'MUC5B'],
    'Stroma': ['DCN', 'COL1A2', 'FBLN1'],
}
