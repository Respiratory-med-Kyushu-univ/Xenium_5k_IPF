Data structure:
1. Raw_Xenium/
   Anonymized Xenium output files.

2. ProSeg_processed_SpatialData/
   SpatialData objects generated from ProSeg-processed Xenium data.
   Each object contains ProSeg-derived cell boundaries and the cell-by-gene count matrix used for downstream analyses.
   Morphology images and transcript-level point data are not included in these SpatialData objects to avoid redundancy and reduce file size.

3. Analysis_ready/
   Final analysis-ready AnnData object and source tables used to generate the manuscript figures.