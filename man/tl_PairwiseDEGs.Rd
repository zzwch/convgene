% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tl.R
\name{tl_PairwiseDEGs}
\alias{tl_PairwiseDEGs}
\title{clac DEGs among multiple groups}
\usage{
tl_PairwiseDEGs(
  data,
  group_by,
  groups = NULL,
  features = NULL,
  p_adjust = T,
  multi_test = kruskal.test,
  do_pairwise = F,
  pairs = NULL,
  pairwise_test = wilcox.test,
  ...
)
}
\arguments{
\item{data}{feature by cell matrix}

\item{group_by}{a vector of cell annotations}

\item{groups}{subset cells belonging to these groups}

\item{features}{features to calc}

\item{p_adjust}{adjust p value}

\item{multi_test}{function to test across multiple groups}

\item{do_pairwise}{do pairwised test}

\item{pairs}{a list of pairs to test, NULL for all pairs}

\item{pairwise_test}{function to test two groups included in pairs.}

\item{...}{ohter params passed to pairwise_test.}
}
\value{
a data.frame of gene and pvalues
}
\description{
clac DEGs among multiple groups
}
\examples{

}
