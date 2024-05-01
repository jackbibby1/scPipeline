#' Perform normalisation, batch correction, clustering, and dimred
#'
#' This function takes a Seurat object as an input and performs normalisation
#' (log or sctransform), batch correction (harmony, RPCA, or CCA), clustering, and
#' dimred (tSNE or UMAP).
#'
#' @param seurat_object Single Seurat object with annotated metadata
#' @param normalisation_method Method for normalisation (LogNormalize or SCTransform from Seurat)
#' @param generate_tsne Should a tSNE be generated?
#' @param batch_correction Should batch correction be done?
#' @param correction_method Which batch correction method to use. Options are: CCAIntegration, RPCAIntegration (default),
#'    HarmonyIntegration, FastMNNIntegration, and scVIIntegration
#'
#' @examples \dontrun{
#'   data <- norm_integration(seurat_object = seurat_data,
#'                         generate_tsne = TRUE,
#'                         batch_correction = TRUE,
#'                         correction_method = "harmony",
#'                         batch_correction_group = "donor")
#' }
#'
#' @return A Seurat object that contains normalised and batch corrected data, as well
#'     as clusterting and dimred information.
#'
#' @export
#'
#'

norm_integration <- function(seurat_object = NULL,
                             normalisation_method = "LogNormalize",
                             generate_tsne = FALSE,
                             batch_correction = TRUE,
                             correction_method = "RPCAIntegration") {


  ##---------- normalisation
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("\n---------- Step 1: Normalisation and scaling \n")

  if (batch_correction == FALSE) {

    cat("--- Not performing batch correction \n")
    cat("--- Joining layers \n")

    seurat_object <- SeuratObject::JoinLayers(seurat_object)

    seurat_object <- normalise_data(seurat_object = seurat_object,
                                    normalisation_method = normalisation_method)

    cat("\n...Done \n")
    return(seurat_object)

  } else {

    seurat_object <- normalise_data(seurat_object = seurat_object,
                                    normalisation_method = normalisation_method)

  }


  ##---------- batch correction

  cat("\n---------- Step 2: Batch correction \n")

  if (batch_correction == TRUE & normalisation_method != "SCT") {

    cat("--- Performing batch correction using", correction_method, "on RNA data", "\n")

    seurat_object <- Seurat::IntegrateLayers(object = seurat_object,
                                             method = correction_method,
                                             orig.reduction = "pca",
                                             new.reduction = paste("integrated", correction_method, sep = "_"),
                                             verbose = FALSE)

    cat("--- Joining layers \n")
    seurat_object <- SeuratObject::JoinLayers(seurat_object)

  } else if (batch_correction == TRUE & normalisation_method == "SCT") {

    cat("--- Performing batch correction using", correction_method, "on SCT data", "\n")

    seurat_object <- Seurat::IntegrateLayers(object = seurat_object,
                                             method = correction_method,
                                             normalization.method = "SCT",
                                             verbose = FALSE)

    cat("\n...Done \n")
    cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")

  }

  return(seurat_object)

}


#' Perform normalisation via log or SCT
#'
#' This function takes a Seurat object as an input and performs normalisation
#' (log or sctransform), batch correction (harmony, RPCA, or CCA), clustering, and
#' dimred (tSNE or UMAP).
#'
#' @param seurat_object Seurat object with non-normalised values
#' @param normalisation_method Method for normalisation (LogNormalize or SCTransform from Seurat)
#' @param num_var_features Number of variable features if using LogNormalize
#'
#' @examples \dontrun{
#'   data <- normalise_data(seurat_object = seu_obj)
#' }
#'
#' @return A Seurat object that contains normalised data
#'
#' @export
#'
#'

normalise_data <- function(seurat_object = NULL,
                           normalisation_method = NULL,
                           num_var_features = NULL) {

  ## normalise via log
  if (normalisation_method == "LogNormalize") {

    cat("--- Normalizing data using LogNormalize \n")

    seurat_object <- Seurat::NormalizeData(seurat_object) %>%
      Seurat::FindVariableFeatures() %>%
      Seurat::ScaleData() %>%
      Seurat::RunPCA(verbose = FALSE)

  } else if (normalisation_method == "SCT") {

    ## normalise by sct
    cat("--- Normalizing data using SCT \n")
    seurat_object <- Seurat::SCTransform(seurat_object,
                                         method = "glmGamPoi",
                                         vst.flavor = "v2",
                                         verbose = FALSE) %>%
      Seurat::RunPCA(verbose = FALSE)

  } else {

    stop("Choose either LogNormalize or SCT for normalisation_method")

  }

}









