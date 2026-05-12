# Shared, project-wide constants (lists/dicts) used across notebooks and scripts.


PROSEG_PATH = "/home/s_nakamura/Proseg/Xenium_IPF/"  # Xenium/ProSeg data directory
WORK_PATH = "/mnt/d/data001/xenium_IPF/work space/Proseg/"  # path to your workspace
PUBLIC_PATH = "/mnt/d/data001/public/"  # path to scRNA-seq public data

SAMPLES_SYMBOLS = ["AF_s2", "FR_s1", "FR_s3", "FR_s4", "HC_s1", "HC_s3"]

MAJOR_CLUSTER_COLORS={
    "Epithelial": "#E41A1C",  
    "Fibroblast": "#377EB8",  
    "SMC/Pericyte": "#4DAF4A", 
    "Endothelial": "#984EA3",  
    "Monocyte/MΦ": "#FF7F00",  
    "NK/T cell": "#FFFF33",  
    "B/Plasma cell": "#E78AC3",  
    "Dendritic cell": "#00FFFF",  
    "Granulocyte": "#66C2A5",  
    "FRC/FDC": "#A65628",
}

CELL_TYPE_3M = {
    "Epithelial": "Epithelial",
    "Fibroblast": "Stromal",
    "SMC/Pericyte": "Stromal",
    "Endothelial": "Stromal",
    "Monocyte/MΦ": "Immune",
    "NK/T cell": "Immune",
    "B/Plasma cell": "Immune",
    "Dendritic cell": "Immune",
    "Granulocyte": "Immune",
    "FRC/FDC": "Stromal", 
}

IMP_GENE_NAMES = {
    "Epithelial": ["SFTPC", "CAV1", "KRT5", "KRT17", "KRT8", "SCGB3A2", "SCGB1A1", "CAPS"],
    "Fibroblast": ["ACTA2", "COL1A1", "ITGA8"],
    "Macrophage": ["SPP1", "PLA2G7", "CHI3L1", "FABP4"]
}

SUBTYPE_HATCHES = {
    "AT2": '//',
    "AT1": '\\\\',
    "KRT5-/KRT17+": "*",
    "Basal": "O",
    "RAS": ".",
    "Goblet": "-",
    "Ciliated": "|",
        
    "Myofibroblast": "*",
    "COMP+ fibroblast": "O",
    "GREM1+ fibroblast": "+",
    "Adventitial fibroblast": "//",
    "Lipofibroblast": ".",
    
    "Smooth muscle cell": ".",
    "Pericyte": "+",
    "Alveolar Pericyte": "O",

    "Lymphatic": ".",
    "Arteriole": '//',
    "Aerocyte": '*',
    "General capillary": 'XX',
    "Venule": '\\\\',
    
    "Monocyte": ".",
    "Interstitial macrophage": "O",
    "SPP1+ macrophage": "*",
    "CXCL9+ macrophage": "//",
    "Alveolar macrophage": "+",
    
    "CD4+ T": '//',
    "CD8+ T": '\\\\',
    "Treg": ".",
    "NK": "*",
    
    "B cell": ".",
    "Plasma cell": "//",
    
    "pDC": "XX",
    "cDC1": "//",
    "cDC2": "\\\\",
    "mregDC": "*",
    
    "Neutrophil": "//",
    "Mast cell": "+",
    
    "FRC": "O",
    "FDC": "+",
}

ALL_SUBTYPE_LIST = [
    "AT2", "AT1", "KRT5-/KRT17+", "Basal", "RAS", "Goblet", "Ciliated",
    "Myofibroblast", "COMP+ fibroblast", "GREM1+ fibroblast", "Adventitial fibroblast", "Lipofibroblast",
    "Smooth muscle cell", "Pericyte", "Alveolar Pericyte",
    "Lymphatic", "Arteriole", "Aerocyte", "General capillary", "Venule",
    "Monocyte", "Interstitial macrophage", "SPP1+ macrophage", "CXCL9+ macrophage", "Alveolar macrophage",
    "CD4+ T", "CD8+ T", "Treg", "NK",
    "B cell", "Plasma cell",
    "pDC", "cDC1", "cDC2", "mregDC", 
    "Neutrophil", "Mast cell",
    "FRC", "FDC"
]

CMAP_COLORLIST=["#1D224C", "#692755", "#B33632", "#FF6622", "#FFFFBB"]

POINTS_DICT = {
    "FR_s1": {
        "HE_points":[[1712.9, 2743.9],[11621.3, 10690.4]],
        "GDF_points":[[854.3, 1334.5],[4583.1, 4317.5]]
    },
    "HC_s1": {
        "HE_points":[[4064.3, 3030.9],[15685.4, 11607.3]],
        "GDF_points":[[1798.9, 1297.4],[6174.6, 4518.0]]
    },
    "AF_s2": {
        "HE_points":[[1937.4, 5302.2],[15992.3, 12154.3]],
        "GDF_points":[[784.9, 2139.2],[6074.8, 4713.4]]
    },
    "FR_s3": {
        "HE_points":[[820.2, 1909.4],[18339.6, 8708.2]],
        "GDF_points":[[511.3, 1021.4],[7106.6, 3563.8]]
    },
    "HC_s3": {
        "HE_points":[[2838.7, 5106.5],[13934.3, 14168.3]],
        "GDF_points":[[1437.0, 1929.9],[5616.6, 5328.5]]
    },
    "FR_s4": {
        "HE_points":[[1111.6, 1157.7],[16567.1, 11930.7]],
        "GDF_points":[[736.8, 389.0],[6556.9, 4426.2]]
    },

}

TITLE_ABBREVIATION = {
    "KRT5-/KRT17+": "KRT5-/17+",
    "Myofibroblast": "Myofib.",
    "COMP+ fibroblast": "COMP+ fib.",
    "GREM1+ fibroblast": "GREM1+ fib.",
    "Adventitial fibroblast": "Adventitial fib.",
    "Lipofibroblast": "Lipofib.",
    "Smooth muscle cell": "SMC",
    "Interstitial macrophage": "Interstitial MΦ",
    "SPP1+ macrophage": "SPP1+ MΦ",
    "CXCL9+ macrophage": "CXCL9+ MΦ",
    "Alveolar macrophage": "Alveolar MΦ",
}

EPITHELIAL_SUBTYPE = ["AT2", "AT1", "KRT5-/KRT17+", "Basal", "RAS", "Goblet", "Ciliated"]
