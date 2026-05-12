Data structure:

1. Raw\_Xenium/
Anonymized Xenium output files.
2. ProSeg\_processed\_SpatialData/
SpatialData objects generated from ProSeg-processed Xenium data.
Each object contains ProSeg-derived cell boundaries and the cell-by-gene count matrix used for downstream analyses.
Morphology images and transcript-level point data are not included in these SpatialData objects to avoid redundancy and reduce file size.
3. Analysis\_ready/
Final analysis-ready AnnData object.

