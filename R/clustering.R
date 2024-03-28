#' Clustering and dimensional reduction
#'
#' This function takes a Seurat object as an input
#' and performs downstream clustering and dimensionality reduction
#'
#' @param seurat_object Seurat object with normalised data
#' @param export_elbow Should the elbow plot be exported? Mainly for use in non-GUI sessions
#' @param generate_tsne Should tSNE be calculated?
#'
#' @examples \dontrun{
#'   data <- seurat_clustering(seurat_object = seurat_object)
#' }
#'
#' @return A Seurat object that contains batch-corrected data
#'
#' @export
#'
#'

seurat_clustering <- function(seurat_object = NULL,
                              reduction = "integrated",
                              export_elbow = FALSE,
                              generate_tsne = FALSE) {

  ##---------- choosing dimensions

  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("\n---------- Step 1: Choosing dimensions \n")

  cat("--- Generating elbow plot of PCs \n")
  print(Seurat::ElbowPlot(seurat_object, ndims = 50))

  if (export_elbow == TRUE) {

    cat("--- Exporting elbow plot to check PC number")
    ggplot2::ggsave("elbow_plot.png", width = 5, height = 5, dpi = 600)

  }

  elbow_value <- readline(prompt = "Choose number of PCs based on elbow plot: ")
  elbow_value <<- as.numeric(elbow_value)


  if (generate_tsne == TRUE) {

    message("--- Calculating tSNE using dims 1:", elbow_value)
    seurat_object <- Seurat::RunTSNE(seurat_object, dims = 1:elbow_value)

  }

  message("Performing UMAP using dims 1:", elbow_value)

  reduction_name <- grep(reduction, names(int_data@reductions), value = T)

  seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:elbow_value, reduction = reduction_name, verbose = FALSE) %>%
    Seurat::FindNeighbors(dims = 1:elbow_value, reduction = reduction_name, verbose = F) %>%
    Seurat::FindClusters(resolution = 0.5, verbose = F)

  cat("\n...Done \n")

  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")

}








