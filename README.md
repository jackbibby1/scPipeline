WIP - don't use

# scPipeline
Pipeline for basic end-to-end scRNA-seq processing

### Simplified functions for end to end processing of scRNA-seq data

#### Functions:

- `pre_process_scrna()` covers: 
  - Reading in data, either cellranger, h5, or tab files
  - Filtering based on mito percentages
  - Adding metadata
  - Merging samples if >1
- `process_scrna()` covers:
  - Normalisation
  - Batch correction (harmony, RPCA, or CCA)
  - Dimensionality reduction (PCA, UMAP, and tSNE)
  - Clustering

