% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{seurat_clustering}
\alias{seurat_clustering}
\title{Clustering and dimensional reduction}
\usage{
seurat_clustering(
  seurat_object = NULL,
  reduction = "integrated",
  export_elbow = FALSE,
  generate_tsne = FALSE
)
}
\arguments{
\item{seurat_object}{Seurat object with normalised data}

\item{export_elbow}{Should the elbow plot be exported? Mainly for use in non-GUI sessions}

\item{generate_tsne}{Should tSNE be calculated?}
}
\value{
A Seurat object that contains batch-corrected data
}
\description{
This function takes a Seurat object as an input
and performs downstream clustering and dimensionality reduction
}
\examples{
\dontrun{
  data <- seurat_clustering(seurat_object = seurat_object)
}

}
