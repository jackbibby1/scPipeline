# scPipeline

- Pipeline for basic end-to-end scRNA-seq processing

#### Installation:

```
devtools::install_github("jackbibby1/scPipeline")
```

#### Main functions:

- `read_scrna()` covers: 
  - Reading in data, either cellranger, h5, or tab files (multicore option)
  - Filtering based on mito percentages
  - Prompt user for input on mito cutoff before proceeding
  - Adding metadata
  - Merging samples if >1
- `norm_integration()` covers:
  - Normalisation (Log or SCT)
  - Getting variable features
  - Scaling
  - Dimensionality reduction (PCA)
  - Batch correction and integration (any method available in `Seurat::IntegrateLayers()`)
- `seurat_clustering()`
  - Defining number of PCs
  - Dimensionality reduction (UMAP and/or tSNE)
  - Clustering

#### Example:

```ruby

# set wdir and generate metadata ------------------------------------------
setwd("path-to/wdir")
folders <- list.dirs("data", recursive = F)
meta <- data.frame(time = rep(c(0, 12, 24), times = 4),
                   cell = rep(c("cd4", "cd8"), each = 6),
                   lineage = rep(c("naive", "memory"), each = 3, times = 2))

# reading in sc data ---------------------------------------------------
# generates single Seurat object containing data merged from all filepath folders, with annotated metadata
df <- read_scrna(filepath = "~/Desktop/test/data",
                 filename_pattern = "GSM",
                 cores = 4, 
                 metadata = meta,
                 mito_pattern = "^MT-",
                 merge_data = TRUE)
                 
# normalisation and batch correction ---------------------------------------------------           
# performs normalisation, var features finding, scaling, PCA and joins layers
df <- norm_integration(df,
                       normalisation_method = "LogNormalize",
                       batch_correction = T,
                       correction_method = "HarmonyIntegration")
                 

# clustering data ---------------------------------------------------
# generates single Seurat object with umap and clustering performed

df <- seurat_clustering(df, reduction = "integrated")
                    
```

