\# Xenium\_5k\_IPF



Code for spatial transcriptomics analysis of Xenium 5k-derived idiopathic pulmonary fibrosis datasets.



The code was used for the manuscript: Nakamura S, et al. Trajectory inference of epithelial-centered neighborhood profiles reconstructs a pseudo-temporal continuum in idiopathic pulmonary fibrosis.



The repository contains scripts for preprocessing, cell-type annotation, epithelial-centered neighborhood profiling, pseudo-temporal continuum inference, pathway activity analysis, spatial cell-cell communication analysis, and figure generation.



Requirements: Xenium/ProSeg-derived cell-level data, including cell-by-gene count matrices, spatial coordinates, cell metadata, cell-type annotations, and sample metadata. The full dataset is deposited in GEO under accession number GSE334099.



The analysis notebooks and scripts were tested in a Linux-based environment using Python and JupyterLab.



Tested environment:



\- OS: Ubuntu 22.04 under WSL2

\- Python: 3.10.x and/or 3.12.x

\- JupyterLab: x.x.x

\- CUDA-enabled GPU environment for steps using CuPy



Major Python packages tested include:



\- scanpy: 1.12.1

\- anndata: 0.12.11

\- pandas: 2.3.3

\- numpy: 2.4.4

\- scipy: 1.17.1

\- scikit-learn: 1.8.0

\- matplotlib: 3.10.9

\- seaborn: 0.13.2

\- squidpy: 1.8.1

\- decoupler: 2.1.6

\- scFates: 1.0.9

\- gseapy: 1.1.6

\- tangram-sc: 1.0.4

\- stAI: source code obtained from the public GitHub repository; no formal version number was available

\- geopandas: 1.1.3

\- shapely: 2.1.2

\- pycirclize: 1.10.1

\- cupy-cuda12x: 13.5.1



No specialized hardware is required for most scripts. Selected GPU-accelerated steps use CuPy and require a CUDA-compatible GPU environment with a compatible CuPy package, such as cupy-cuda12x. Analyses on the complete Xenium/ProSeg dataset may also require sufficient RAM and storage.

