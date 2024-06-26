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
#'    human "^MT-".
#' @param cores Number of cores to use in mclapply() for reading data
#' @param plot_mito Produces violin plots for the mito percentages so you can choose a cutoff
#' @param merge_data Should data be merged into a single Seurat object after processing?
#' @param metadata Metadata file for samples. This should be formatted as one sample per row, and
#'    information regarding that sample in subsequent columns e.g.
#'    data.frame(sample = c("sample1", "sample2", "sample3", "sample4"),
#'               disease = c("healthy", "healthy", "disease", "disease"),
#'               sorted population = c("t_cell", "b_cell", "t_cell", "b_cell"))
#'
#' @examples \dontrun{
#'   data <- read_scrna(filepath = "raw_data",
#'                             file_type = "cellranger",
#'                             filename_pattern = "filtered",
#'                             metadata = "sample_metadata")
#' }
#'
#' @return Either a list of Seurat objects of the processed data if merge = FALSE, or a merged Seurat
#'     object containing processed and filtered data if merge = TRUE.
#'
#' @export

read_scrna <- function(filepath = NULL,
                       file_type = "cellranger",
                       filename_pattern = NULL,
                       mito_pattern = "^MT-",
                       cores = 8,
                       cite_seq = FALSE,
                       plot_mito = TRUE,
                       filter_mito = TRUE,
                       merge_data = TRUE,
                       metadata = NULL,
                       min_cells = 10,
                       min_features = 200) {

  ##---------- read in specified file type

  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  cat("\n---------- Step 1: Reading in data \n")
  `%notin%` <- Negate(`%in%`)

  if (file_type %notin% c("cellranger", "h5", "tab")) {
    stop("For file_type, choose between cellranger, h5, or tab", call. = F)
  }

  df <- reading_data(filepath = filepath,
                     file_type = file_type,
                     filename_pattern = filename_pattern,
                     cores = cores)

  if (file_type == "cellranger") {

    cat("  ...Adding folder names to metadata information \n")
    updated_metadata <- metadata %>% dplyr::mutate(filepath = folders)

  } else if (file_type %in% c("h5", "tab")) {

    message("Adding file names to metadata information")
    updated_metadata <- metadata %>% dplyr::mutate(filepath = files)

  }

  ##---------- process and filter the data

  cat("\n---------- Step 2: Processing data \n")
  cat("--- Creating the Seurat object and adding mito percentages \n")

  if (cite_seq == TRUE) {

    cat("--- Expecting CITE-Seq data \n")

    df <- parallel::mclapply(df, function(x) {

      combined_data <- Seurat::CreateSeuratObject(x[["Gene Expression"]], min.cells = min_cells, min.features = min_features)
      combined_data[["ADT"]] <- Seurat::CreateAssayObject(x[["Antibody Capture"]][, colnames(combined_data)])
      combined_data <- Seurat::PercentageFeatureSet(combined_data, pattern = mito_pattern, col.name = "percent_mt")
      return(combined_data)

    }, mc.cores = cores)

  } else {

    ## create the seurat object
    df <- parallel::mclapply(df, function(x) {

      Seurat::CreateSeuratObject(x, min.cells = min_cells, min.features = min_features) %>%
        Seurat::PercentageFeatureSet(pattern = mito_pattern, col.name = "percent_mt")

    }, mc.cores = cores)

  }

  if (plot_mito == TRUE) {

    cat("--- Plotting the mito percentages \n")
    ## plot the mito percentages for all cells
    plots <- lapply(df, function(x) {

      suppressMessages(Seurat::VlnPlot(x, features = "percent_mt", pt.size = 0.1, slot = "counts") +
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
    cat("--- Filtering based on mito percentages \n")
    mito_percent <- readline(prompt = "Choose % mito cutoff: ")
    mito_percent <- as.numeric(mito_percent)
    cat("Using", mito_percent, "% as the mito cutoff \n")
    df <- lapply(df, function(x) subset(x, percent_mt < mito_percent))

  }

  ##---------- adding metadata to the samples

  cat("\n---------- Step 3: Adding metadata to the Seurat object based on: \n")
  print(updated_metadata)

  meta_names <- colnames(updated_metadata)
  for (i in 1:nrow(updated_metadata)) {
    for (y in meta_names) {

      df[[i]] <- Seurat::AddMetaData(df[[i]], col.name = y, metadata = updated_metadata[, y][i])

    }
  }

  ##---------- merging the data

  cat("\n---------- Step 4: Merging the data \n")

  if (merge_data == TRUE) {

    cat("--- Merging samples 1:", length(df), sep = "")
    df <- suppressWarnings(merge(x = df[[1]], y = df[2:length(df)]))

  } else {

    cat("--- Not merging samples \n")

    }

  cat("\n...Done \n")
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")

  return(df)

}


#' Reading in single cell data
#'
#' This function processes single cell data that is either formatted as
#' cellranger output, h5, or tab delim files.
#'
#' @param filepath File path to the directory containing the data
#' @param file_type Type of files for input. Either "cellranger", "h5", or "tab"
#' @param filename_pattern Regex for the filenames or folders in the filepath. Used if you only want to process
#'    some files or folders e.g. "*h5" for .h5 files.
#'
#' @examples \dontrun{
#'   data <- reading_data(seurat_object = seu_obj,
#'                        file_type = "cellranger",
#'                        filename_pattern = "*h5")
#' }
#'
#' @return A Seurat object that contains batch corrected data
#'
#' @export
#'
#'


reading_data <- function(filepath = NULL,
                         file_type = NULL,
                         filename_pattern = NULL,
                         cores = NULL) {

  cat("--- Using", cores, "cores \n")

  if (file_type == "cellranger") {

    ## define folders to use
    folders <<- list.dirs(path = filepath, full.names = T, recursive = F)

    if (!is.null(filename_pattern)) {
      folders <<- grep(x = folders, pattern = filename_pattern, ignore.case = T, value = T)
    }

    cat("--- Folders to process are: ", folders, sep = "\n")
    cat("  ...Processing cellranger data\n")

    ## read the cellranger files
    cat("  ...Reading in data with Read10X()\n")
    df <- parallel::mclapply(folders, function(x) Seurat::Read10X(data.dir = x), mc.cores = cores)

  } else if (file_type == "h5") {

    ## define files to process
    files <<- list.files(path = filepath, full.names = T)

    if (!is.null(filename_pattern)) {
      files <<- grep(x = files, pattern = filename_pattern, ignore.case = T, value = T)
    }

    cat("--- Files to process are: ", files, sep = "\n")
    message("  ...Processing .h5 data")

    message("  ...Reading in data with `Read10X_h5`")
    df <- parallel::mclapply(files, function(x) Seurat::Read10X_h5(filename = x), mc.cores = cores)

  } else if (file_type == "tab") {

    ## define files to process
    files <<- list.files(path = filepath, full.names = T)

    if (!is.null(filename_pattern)) {
      files <<- grep(x = files, pattern = filename_pattern, ignore.case = T, value = T)
    }

    cat("--- Files to process are: ", files, sep = "\n")
    message("  ...Reading in data with `fread`")

    df <- parallel::mclapply(files, function(x) as.matrix(data.table::fread(file = x), rownames = 1), mc.cores = cores)

  } else {

    stop("Filetype not recognized. Choose between cellranger, h5, or tab")

  }

}






