#' Perform normalisation, batch correction, clustering, and dimred
#'
#' This function takes a Seurat object as an input and performs normalisation
#' (log or sctransform), batch correction (harmony, RPCA, or CCA), clustering, and
#' dimred (tSNE or UMAP).
#'
#' @param seurat_object Single Seurat object with annotated metadata
#' @param normalisation_method Method for normalisation (LogNormalize or SCTransform from Seurat)
#' @param num_var_features Number of variable features if using LogNormalize
#' @param num_sct_features Number of features to use if using SCTransform
#' @param generate_tsne Should a tSNE be generated?
#' @param batch_correction Should batch correction be done?
#' @param correction_method Which batch correction method to use. Either "harmony", "rpca", "cca", or "scvi"
#' @param batch_correction_group Group used to correct for batch effects. Grouping data should
#'    reference a column in the Seurat object e.g. "donor" or "donor_and_stim"
#' @param nintegration_features Number of features used for integration
#' @param integration_strength Strength of integration in RPCA. Refers to k
#'
#' @examples \dontrun{
#'   data <- process_scrna(seurat_object = seurat_data,
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

process_scrna <- function(seurat_object = NULL,
                          normalisation_method = "SCT",
                          num_var_features = 2000,
                          num_sct_features = 3000,
                          generate_tsne = FALSE,
                          batch_correction = TRUE,
                          export_elbow = TRUE,
                          correction_method = "harmony",
                          batch_correction_group = NULL,
                          nintegration_features = 3000,
                          integration_strength = 5) {


  `%notin%` <- Negate(`%in%`)
  cat("---------- Running pipeline for initial normalisation and scaling \n")

  if (batch_correction == TRUE & correction_method %in% c("cca", "scvi")) {
    cat("--- Data will be normalised later during integration \n")
  }

  if (batch_correction == FALSE |
      batch_correction == TRUE &
      correction_method == "harmony" |
      correction_method == "rpca") {

    if (batch_correction == TRUE & correction_method %in% c("harmony", "rpca")) {
      cat("--- Performing SCTransform and PCA analysis for harmony or RPCA \n")
    } else if (batch_correction == FALSE) {
      cat("--- Performing SCTransform and PCA analysis without batch correction \n")
    }

    ##---------- normalisation

    seurat_object <- normalise_data(seurat_object = seurat_object,
                                    normalisation_method = normalisation_method)

  }

  ##---------- choosing dimensions

  if (batch_correction == TRUE & correction_method %notin% c("cca", "scvi")) {
    message("Generating elbow plot of PCs")
    print(Seurat::ElbowPlot(seurat_object, ndims = 50))

    if (export_elbow == TRUE) {

      cat("--- Exporting elbow plot to check PC number")
      ggplot2::ggsave("elbow_plot.png", width = 5, height = 5, dpi = 600)

    }

    elbow_value <- readline(prompt = "Choose number of PCs based on elbow plot: ")
    elbow_value <<- as.numeric(elbow_value)
  }

  ##---------- batch correction

  if (batch_correction == TRUE) {

    cat("\n---------- Running pipeline for batch correction \n")

    seurat_object <- batch_correction(seurat_object = seurat_object,
                                      correction_method = correction_method,
                                      batch_correction_group = batch_correction_group,
                                      nintegration_features = nintegration_features,
                                      integration_strength = integration_strength)

  }

  ##---------- umap/tsne and clustering

  cat("\n---------- Downstream clustering and dimensionality reduction pipeline \n")

  if (batch_correction == TRUE & correction_method == "harmony") {

    cat("--- Running downstream harmony pipeline \n")

    seurat_object <- harmony_scvi_clustering(seurat_object = seurat_object,
                                             generate_tsne = generate_tsne,
                                             correction_method = correction_method)

  } else {

    if (batch_correction == TRUE & correction_method %in% c("rpca", "cca")) {
      cat("--- Running RPCA/CCA pipeline \n")
    } else {
      cat("--- Running standard clustering pipeline \n")
    }

    seurat_object <- seurat_clustering(seurat_object = seurat_object,
                                       generate_tsne = generate_tsne,
                                       export_elbow = export_elbow)

  }

  message("\n...Done")
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
    seurat_object <- lapply(seurat_object, function(x) {

      x %>% Seurat::NormalizeData() %>%
        Seurat::FindVariableFeatures(nfeatures = num_var_features) %>%
        Seurat::ScaleData() %>%
        Seurat::RunPCA(verbose = FALSE)

    })

  } else if (normalisation_method == "SCT") {

    ## normalise by sct
    cat("--- Normalizing data using SCT \n")
    seurat_object <- Seurat::SCTransform(seurat_object,
                                         method = "glmGamPoi",
                                         vst.flavor = "v2") %>%
      Seurat::RunPCA(verbose = FALSE)

  } else {

    stop("Choose either LogNormalize or SCT for normalisation_method")

  }

}


#' Batch correction
#'
#' This function takes a Seurat object as an input and performs batch
#' correction using either harmony, RPCA, or CCA
#'
#' @param seurat_object Seurat object with normalised and scaled values, typically with
#' PCA calculated
#' @param correction_method Which batch correction method to use. Either "harmony", "rpca", or "cca"
#' @param batch_correction_group Group used to correct for batch effects. Grouping data should
#'    reference a column in the Seurat object e.g. "donor" or "donor_and_stim"
#' @param nintegration_features Number of features used for integration
#' @param integration_strength Strength of integration in RPCA. Refers to k
#' @param num_sct_features Number of features to use if using SCTransform
#'
#' @examples \dontrun{
#'   data <- batch_correction(seurat_object = seu_obj,
#'                            batch_correction = TRUE,
#'                            correction_method = "harmony",
#'                            batch_correction_group = "sample")
#' }
#'
#' @return A Seurat object that contains batch corrected data
#'
#' @export
#'
#'


batch_correction <- function(seurat_object = NULL,
                             correction_method = NULL,
                             batch_correction_group = NULL,
                             nintegration_features = NULL,
                             num_sct_features = NULL,
                             integration_strength = NULL) {

  if (correction_method == "harmony") {

    cat("--- Batch correcting data using Harmony based on metadata group:", batch_correction_group, "\n")
    ## run harmony on calculated pca values
    seurat_object <- harmony::RunHarmony(seurat_object,
                                         group.by.vars = batch_correction_group,
                                         plot_convergence = T)

  } else if (correction_method == "scvi") {

    cat("--- Batch correcting data using scVI based on metadata group:", batch_correction_group, "\n")
    cat("--- scVI model to be stored in dimred scvi slot", "\n")

    ## import python packages
    sc <- reticulate::import("scanpy", convert = FALSE)
    scvi <- reticulate::import("scvi", convert = FALSE)

    adata <- sceasy::convertFormat(seurat_object,
                           from = "seurat", to = "anndata",
                           main_layer = "counts", drop_single_values = FALSE)

    # run setup_anndata
    scvi$model$SCVI$setup_anndata(adata, batch_key = batch_correction_group)

    # create the model
    model <- scvi$model$SCVI(adata)

    # train the model
    model$train()

    # get the latent representation
    latent <- model$get_latent_representation()

    # store it in the Seurat object
    latent <- as.matrix(latent)
    rownames(latent) <- colnames(seurat_object)
    rm(adata)
    seurat_object[["scvi"]] <- Seurat::CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_object))

  } else if (correction_method == "rpca") {

    cat("--- Batch correcting data using Seurat RPCA based on metadata group:", batch_correction_group, "\n")
    cat("--- Data to be split across batch and normalised independently \n")
    ## reprocess data so rpca can be run
    seurat_list <- Seurat::SplitObject(seurat_object, split.by = batch_correction_group)
    seurat_list <- lapply(seurat_list, function(x) Seurat::SCTransform(x, method = "glmGamPoi"))
    features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list)
    seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
    seurat_list <- lapply(seurat_list, function(x) Seurat::RunPCA(x, features = features, verbose = FALSE))

    anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                              normalization.method = "SCT",
                                              anchor.features = features,
                                              dims = 1:elbow_value,
                                              reduction = "rpca",
                                              k.anchor = integration_strength)

    seurat_object <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:elbow_value)


  } else if (correction_method == "cca") {

    cat("--- Batch correcting data using Seurat CCA based on metadata group:", batch_correction_group, "\n")
    cat("--- Data to be split across batch and normalised independently \n")
    seurat_list <- Seurat::SplitObject(seurat_object, split.by = batch_correction_group)
    seurat_list <- lapply(seurat_list, function(x) Seurat::SCTransform(x, method = "glmGamPoi"))
    features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nintegration_features)
    seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
    anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                              normalization.method = "SCT",
                                              anchor.features = features)
    seurat_object <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")


  }


}


#' Harmony-scVI downstream pipeline
#'
#' This function takes a Seurat object with harmony or scVI dimred matrix
#' and performs downstream clustering and dimensionality reduction
#'
#' @param seurat_object Seurat object with normalised data and harmony corrected PCs
#' @param generate_tsne Should tSNE be calculated?
#'
#' @examples \dontrun{
#'   data <- harmony_scvi_clustering(seurat_object = seu_obj)
#' }
#'
#' @return A Seurat object that contains batch corrected data
#'
#' @export
#'
#'

harmony_scvi_clustering <- function(seurat_object = NULL,
                                    generate_tsne = NULL,
                                    correction_method = NULL) {

  if (generate_tsne == TRUE) {

    message("Calculating tSNE using dims 1:", elbow_value)
    seurat_object <- Seurat::RunTSNE(seurat_object,
                                     dims = 1:elbow_value,
                                     reduction = correction_method)

  }

  message("Performing UMAP using dims 1:", elbow_value)

  seurat_object <- Seurat::RunUMAP(seurat_object,
                                   dims = 1:elbow_value,
                                   reduction = correction_method,
                                   verbose = FALSE) %>%
    Seurat::FindNeighbors(dims = 1:elbow_value, reduction = correction_method, verbose = F) %>%
    Seurat::FindClusters(resolution = 0.5, verbose = F)

}


#' RPCA/CCA downstream pipeline
#'
#' This function takes a Seurat object as an input
#' and performs downstream clustering and dimensionality reduction
#'
#' @param seurat_object Seurat object with normalised data
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
                              export_elbow = NULL,
                              generate_tsne = NULL) {

  seurat_object <- Seurat::RunPCA(seurat_object, verbose = FALSE)
  print(Seurat::ElbowPlot(seurat_object, ndims = 50))

  if (export_elbow == TRUE) {

    message("Exporting elbow plot to working directory")
    ggplot2::ggsave("elbow_plot.png", width = 5, height = 5, dpi = 600)

  }

  elbow_value <- readline(prompt = "Choose number of PCs based on elbow plot: ")
  elbow_value <<- as.numeric(elbow_value)

  if (generate_tsne == TRUE) {

    message("Calculating tSNE using dims 1:", elbow_value)
    seurat_object <- Seurat::RunTSNE(seurat_object, dims = 1:elbow_value)

  }

  message("Performing UMAP using dims 1:", elbow_value)

  seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:elbow_value, verbose = FALSE) %>%
    Seurat::FindNeighbors(dims = 1:elbow_value, verbose = F) %>%
    Seurat::FindClusters(resolution = 0.5, verbose = F)

}







