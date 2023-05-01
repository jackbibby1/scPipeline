#' Read in and filter data based on mt percentages
#'
#' This function takes an input of cellranger, h5, or tab delim expression
#' data, filters based on mt percentages, adds metadata to the samples,
#' and merges all the samples in its output.
#'
#' @param filepath File path to the directory containing the data
#' @param file_type Type of files for input. Either "cellranger", "h5", or "tab"
#' @param filename_pattern Regex for the filenames or folders in the filepath. Used if you only want to process
#'    some files or folders e.g. "*h5" for .h5 files.
#' @param mito_pattern Gene pattern to recognise mitochondrial genes in the input. Defaults to
#'    "^MT-".
#' @param plot_mito Produces violin plots for the mito percentages so you can choose a cutoff
#' @param merge_data Should data be merged into a single Seurat object after processing?
#' @param add_metadata Should metadata be added to the final object?
#' @param metadata Metadata file for samples. This should be formatted as one sample per row, and
#'    information regarding that sample in subsequent columns e.g.
#'    data.frame(sample = c("sample1", "sample2", "sample3", "sample4"),
#'               disease = c("healthy", "healthy", "disease", "disease"),
#'               sorted population = c("t_cell", "b_cell", "t_cell", "b_cell"))
#'
#' @examples \dontrun{
#'   data <- pre_process_scrna(filepath = "raw_data",
#'                             file_type = "cellranger",
#'                             filename_pattern = "filtered",
#'                             metadata = "sample_metadata")
#' }
#'
#' @return Either a list of Seurat objects of the processed data if merge = FALSE, or a merged Seurat
#'     object containing processed and filtered data if merge = TRUE.
#'
#' @export

pre_process_scrna <- function(filepath = NULL,
                              file_type = "cellranger",
                              filename_pattern = NULL,
                              mito_pattern = "^MT-",
                              plot_mito = TRUE,
                              merge_data = TRUE,
                              add_metadata = TRUE,
                              metadata = NULL) {

  ##---------- read in specified file type

  `%notin%` <- Negate(`%in%`)

  if (file_type %notin% c("cellranger", "h5", "tab")) {
    stop("For cell_type, choose between cellranger, h5, or tab")
  }

  if (file_type == "cellranger") {

    message("Processing cellranger data")
    ## define folders to use
    folders <- list.dirs(path = filepath, full.names = T, recursive = F)

    if (!is.null(filename_pattern)) {
      folders <- grep(x = folders, pattern = filename_pattern, ignore.case = T, value = T)
    }

    message("Adding folder names to metadata information")
    metadata <- metadata %>% dplyr::mutate(filename = folders)

    cat("----- Folders to process are: ", folders, sep = "\n")

    ## read the cellranger files
    message("Reading in data with `Read10X`")
    df <- pbapply::pblapply(folders, function(x) Seurat::Read10X(data.dir = x))

  } else if (file_type == "h5") {

    message("Processing .h5 data")
    ## define files to process
    files <- list.files(path = filepath, full.names = T)

    if (!is.null(filename_pattern)) {
      files <- grep(x = files, pattern = filename_pattern, ignore.case = T, value = T)
    }

    message("Adding file names to metadata information")
    metadata <- metadata %>% dplyr::mutate(filename = files)

    cat("----- Files to process are: ", files, sep = "\n")

    message("Reading in data with `Read10X_h5`")
    df <- pbapply::pblapply(files, function(x) Seurat::Read10X_h5(filename = x))

  } else if (file_type == "tab") {

    message("Reading in data with `read.delim`")
    ## define files to process
    files <- list.files(path = filepath, full.names = T)

    if (!is.null(filename_pattern)) {
      files <- grep(x = files, pattern = filename_pattern, ignore.case = T, value = T)
    }

    message("Adding file names to metadata information")
    metadata <- metadata %>% dplyr::mutate(filename = files)

    cat("----- Files to process are: ", files, sep = "\n")

    df <- pbapply::pblapply(folders, function(x) utils::read.delim(file = x, header = T, quote = F, sep = "\t"))

  } else {

    stop("Filetype not recognized. Choose between cellranger, h5, or tab")

  }

  ##---------- process and filter the data

  message("Creating the Seurat object and adding mito percentages")
  ## create the seurat object
  df <- lapply(df, function(x) {

    Seurat::CreateSeuratObject(x) %>%
      Seurat::PercentageFeatureSet(pattern = mito_pattern, col.name = "percent_mt")

  })

  message("Plotting the mito percentages")
  ## plot the mito percentages for all cells
  plots <- lapply(df, function(x) {

    suppressMessages(Seurat::VlnPlot(x, features = "percent_mt", pt.size = 0.1) +
                       ggplot2::geom_hline(yintercept = seq(0, 100, 10),
                                           col = "gray40",
                                           linetype = "dashed",
                                           linewidth = 0.3) +
                       ggplot2::scale_y_continuous(breaks = seq(0, 100, 10), limits = c(-5, 100)) +
                       ggplot2::theme(legend.position = "none",
                                      axis.text.x = ggplot2::element_blank(),
                                      axis.title.x = ggplot2::element_blank(),
                                      plot.title = ggplot2::element_blank()))

  })
  print(patchwork::wrap_plots(plots))

  ## filter based on mito percentages
  mito_percent <- readline(prompt = "Choose % mito cutoff: ")
  mito_percent <- as.numeric(mito_percent)
  message("Using ", mito_percent, "% as the mito cutoff")
  df <- lapply(df, function(x) subset(x, percent_mt < mito_percent))

  ##---------- adding metadata to the samples

  message("Adding metadata from the supplied information: ")
  print(metadata)

  meta_names <- colnames(metadata)
  for (i in 1:nrow(metadata)) {
    for (y in meta_names) {

      df[[i]] <- Seurat::AddMetaData(df[[i]], col.name = y, metadata = metadata[, y][i])

    }
  }

  ##---------- merging the data

  if (merge_data == TRUE) {

    message("Merging the data")

    df <- merge(df[[1]],
                df[2:length(df)])

  }

  return(df)

}


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
#' @param correction_method Which batch correction method to use. Either "harmony", "rpca", or "cca"
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
                          correction_method = "cca",
                          batch_correction_group = NULL,
                          nintegration_features = 3000,
                          integration_strength = 5) {


  if (batch_correction == FALSE |
      batch_correction == TRUE &
      correction_method == "harmony" |
      correction_method == "rpca") {

    cat("----- Initial normalisation and scaling pipeline \n")

    if (batch_correction == TRUE & correction_method == "harmony") {
      message("Performing normalisation and scaling for harmony or RPCA")
    } else if (batch_correction == FALSE) {
      message("Performing normalisation without batch correction")
    }

    ##---------- normalisation

    ## normalise via log
    if (normalisation_method == "LogNormalize") {

      message("Normalizing data using LogNormalize")
      seurat_object <- lapply(seurat_object, function(x) {

        x %>% Seurat::NormalizeData() %>%
          Seurat::FindVariableFeatures(nfeatures = num_var_features) %>%
          Seurat::ScaleData() %>%
          Seurat::RunPCA(verbose = FALSE)

      })

    } else if (normalisation_method == "SCT") {

      ## normalise by sct
      message("Normalizing data using SCT")
      seurat_object <- Seurat::SCTransform(seurat_object, method = "glmGamPoi") %>%
        Seurat::RunPCA(verbose = FALSE)

    } else {

      stop("Choose either LogNormalize or SCT for normalisation_method")

    }

  }


  ##---------- choosing dimensions

  if (batch_correction == TRUE & correction_method != "cca") {
    message("Generating elbow plot of PCs")
    print(Seurat::ElbowPlot(seurat_object, ndims = 50))
    elbow_value <- readline(prompt = "Choose number of PCs based on elbow plot: ")
    elbow_value <- as.numeric(elbow_value)
  }


  ##---------- batch correction

  if (batch_correction == TRUE) {

    cat("----- Batch correction pipeline \n")

    if (correction_method == "harmony") {

      message("Batch correcting data using Harmony based on ", batch_correction_group)
      ## run harmony on calculated pca values
      seurat_object <- harmony::RunHarmony(seurat_object,
                                           group.by.vars = batch_correction_group,
                                           plot_convergence = T)

    } else if (correction_method == "rpca") {

      message("Batch correcting data using Seurat RPCA based on ", batch_correction_group)
      ## reprocess data so rpca can be run
      seurat_list <- Seurat::SplitObject(seurat_object, split.by = batch_correction_group)
      seurat_list <- lapply(seurat_list, function(x) Seurat::SCTransform(x, method = "glmGamPoi"))
      features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = num_sct_features)
      seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
      seurat_list <- lapply(seurat_list, function(x) Seurat::RunPCA(x, features = features, verbose = FALSE))

      anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                                normalization.method = "SCT",
                                                anchor.features = features,
                                                dims = 1:elbow_value,
                                                reduction = "rpca",
                                                k.anchor = integration_strength)

      seurat_object <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:elbow_value)

      ## remove bulky objects
      rm(anchors, seurat_list, features)

    } else if (correction_method == "cca") {

      message("Batch correcting data using Seurat CCA based on ", batch_correction_group)
      seurat_list <- Seurat::SplitObject(seurat_object, split.by = batch_correction_group)
      seurat_list <- lapply(seurat_list, function(x) Seurat::SCTransform(x, method = "glmGamPoi"))
      features <- Seurat::SelectIntegrationFeatures(object.list = seurat_list, nfeatures = nintegration_features)
      seurat_list <- Seurat::PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
      anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_list,
                                                normalization.method = "SCT",
                                                anchor.features = features)
      seurat_object <- Seurat::IntegrateData(anchorset = anchors, normalization.method = "SCT")

      ## remove bulky objects
      rm(anchors, seurat_list, features)

    }

  }

  ##---------- umap/tsne and clustering

  cat("----- Downstream clustering and dimensionality reduction pipeline \n")

  if (batch_correction == TRUE & correction_method == "harmony") {

    message("Running downstream harmony pipeline")

    if (generate_tsne == TRUE) {

      message("Calculating tSNE using dims 1:", elbow_value)
      seurat_object <- Seurat::RunTSNE(seurat_object, dims = 1:elbow_value)

    }

    message("Calculating UMAP using dims 1:", elbow_value)

    seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:elbow_value) %>%
      Seurat::FindNeighbors(dims = 1:elbow_value) %>%
      Seurat::FindClusters(resolution = 0.5)

  } else if (batch_correction == TRUE & correction_method %in% c("rpca", "cca")) {

    message("Running RPCA/CCA pipeline")

    seurat_object <- Seurat::RunPCA(seurat_object, verbose = FALSE)
    print(Seurat::ElbowPlot(seurat_object, ndims = 50))
    elbow_value <- readline(prompt = "Choose number of PCs based on elbow plot: ")
    elbow_value <- as.numeric(elbow_value)

    seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:elbow_value) %>%
      Seurat::FindNeighbors(dims = 1:elbow_value) %>%
      Seurat::FindClusters(resolution = 0.5)

  }

  return(seurat_object)

}



