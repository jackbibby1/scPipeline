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
  - Batch correction and integration (harmony, RPCA, or CCA)
  - Dimensionality reduction (PCA, UMAP, and tSNE)
  - Clustering

#### Example:

```ruby
# set wdir and generate metadata ------------------------------------------
setwd("path-to/wdir")
folders <- list.dirs("data", recursive = F)
metadata <- data.frame(donor = str_extract(string = folders, pattern = "d[0-9]"),
                       disease_status = str_extract(string = folders, pattern = "healthy|disease"))

# pre-process data --------------------------------------------------------
df <- pre_process_scrna(filepath = "data",
                        file_type = "cellranger",
                        filename_pattern = "filtered",
                        mito_pattern = "^MT-",
                        plot_mito = TRUE,
                        add_metadata = TRUE,
                        metadata = metadata,
                        merge_data = TRUE)

# downstream processing ---------------------------------------------------
df <- process_scrna(seurat_object = df,
                    normalisation_method = "SCT",
                    num_sct_features = 3000,
                    generate_tsne = TRUE,
                    batch_correction = TRUE,
                    correction_method = "harmony",
                    batch_correction_group = "donor")
```

