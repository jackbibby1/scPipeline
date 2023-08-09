library(Seurat)
library(tidyverse)
devtools::load_all()

files <- list.files("~/My Drive/analysis_for_people/jubayer/colitis_analysis/E-MTAB-9492_arthritis/raw_data")

metadata <- data.frame(sample = files,
                       tissue = c("sm", "sm", "blood", "sf", "blood", "sf", "blood", "sf"))
metadata <- metadata %>% filter(tissue == "sm")

test <- pre_process_scrna(filepath = "~/My Drive/analysis_for_people/jubayer/colitis_analysis/E-MTAB-9492_arthritis/raw_data",
                          filename_pattern = "4040|BX254",
                          metadata = metadata)

batch_data <- process_scrna(seurat_object = test,
                            batch_correction = FALSE)

cca_data <- process_scrna(seurat_object = test,
                          batch_correction = TRUE,
                          correction_method = "scvi",
                          batch_correction_group = "sample")

rpca_data <- process_scrna(seurat_object = test,
                           batch_correction = TRUE,
                           correction_method = "rpca",
                           batch_correction_group = "sample")

harmony_data <- process_scrna(seurat_object = test,
                              batch_correction = TRUE,
                              correction_method = "harmony",
                              batch_correction_group = "sample")


reticulate::use_condaenv(condaenv = "r-reticulate", conda = "/opt/anaconda3/bin/conda")
reticulate::py_install(packages = c("scanpy", "scvi"))




