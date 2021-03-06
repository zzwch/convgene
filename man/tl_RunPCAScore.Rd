% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tl.R
\name{tl_RunPCAScore}
\alias{tl_RunPCAScore}
\title{a robust PCA for Seurat}
\usage{
tl_RunPCAScore(
  object,
  slot = "data",
  assay = NULL,
  reduction = "pca",
  use_all_genes = T,
  topn = 30,
  do.score.scale = T,
  reduction.name = "rpca"
)
}
\arguments{
\item{object}{Seurat3 object}

\item{slot}{use which data to calc rpca}

\item{assay}{use which assay. NULL to use DefaultAssay(object)}

\item{reduction}{use which reduction to obtain top features.}

\item{use_all_genes}{get top n genes from all genes loadings, not just typically highly variable genes. This will call ProjectDim if inexist.}

\item{topn}{robust top n genes of loadings}

\item{do.score.scale}{scale robust scores}

\item{reduction.name}{new reduction name}
}
\value{
updated Seurat3 object with a new pca_score reduction added
}
\description{
Based on a pre-computed dimensional reduction (typically calculated on a subset of genes and projects this onto the entire dataset (all genes).
Then generate a new rPCA reduction.
}
\examples{

}
