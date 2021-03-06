% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pl.R
\name{pl_dimplot}
\alias{pl_dimplot}
\title{Dimensional reduction plot}
\usage{
pl_dimplot(
  object,
  group_by,
  shape_by = NULL,
  ncol = 2,
  slot = "data",
  subset_by = NULL,
  subset_groups = NULL,
  subset_cells = NULL,
  reduction = "umap",
  dims = c(1, 2),
  pt.size = 1,
  cols = scanpy_colors$default_64
)
}
\arguments{
\item{object}{Seurat object}

\item{group_by}{Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
You can also plot features with gradient cols.}

\item{shape_by}{one column name to control point shape}

\item{ncol}{Number of columns for display when faceting plots}

\item{slot}{FetchData of features from which slot.}

\item{subset_by}{enhanced option, subset cells by the colname.}

\item{subset_groups}{enhanced option, subset cells of these categories.}

\item{subset_cells}{A vector of cell names. if not NULL, subset these cells, otherwise use subset_by and subset_groups}

\item{reduction}{Which dimensionality reduction to use. e.g. umap, tsne, pca}

\item{dims}{Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions}

\item{pt.size}{Adjust point size for plotting}

\item{cols}{Vector of colors, each color corresponds to an identity class.}
}
\value{
a ggplot object
}
\description{
Enhanced DimPlot of Seurat v3
}
\examples{

}
