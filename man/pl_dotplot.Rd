% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pl.R
\name{pl_dotplot}
\alias{pl_dotplot}
\title{Use dotplot/bubbleplot to visualize average/median and percentage expression of features (eg. genes) in groups}
\usage{
pl_dotplot(
  data,
  features,
  group_by,
  groups = NULL,
  method = base::mean,
  by_log = F,
  log_base = base::exp(1),
  exclude_zero = F,
  outlier_cutoff = 1,
  cap_value = NULL,
  swap_axes = F,
  do_scale = F,
  plot_gene_force = F,
  duplicate_features = T,
  gradientn_colors = c("grey90", "grey", "orange", "red", "brown"),
  title = NULL
)
}
\arguments{
\item{data}{expression matrix of feature by cell, ususally normalized (in log scale) but not feature scaled.}

\item{features}{features to be visualized}

\item{group_by}{a vector including cell annotations}

\item{groups}{only these selected groups to be visualized}

\item{method}{using base::mean or stats::median or your custom function for each feature vector;}

\item{by_log}{calculate by log scale ; If TRUE log_base will be ignored.}

\item{log_base}{base::exp(1) or 2 or other else, depends on your normalization method.}

\item{exclude_zero}{exclude zeros when calculate average/median feature value. Note exclude_zero first, outlier_cutoff second.}

\item{outlier_cutoff}{sometimes outliers (several extremely high cells) should be excluded when do summarise. Set 0.99 to exclude top 1 percent cells. (default: 1)}

\item{cap_value}{the max value to show in dotplot. Any value larger than it will be capped to this value. (default: NA)}

\item{swap_axes}{swap x-axis and y-axis or not. (default: False)}

\item{do_scale}{scale the mean value across groups. Useful to demonstrate the differences among groups.}

\item{plot_gene_force}{Force plot gene, which is missing in data.}

\item{duplicate_features}{whether display duplicate features or not.}

\item{gradientn_colors}{the color vector for visualization.}

\item{title}{title to display in the top.}
}
\value{
a ggplot object
}
\description{
Use dotplot/bubbleplot to visualize average/median and percentage expression of features (eg. genes) in groups
}
\examples{

}
