% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pl.R
\name{pl_vioboxplot}
\alias{pl_vioboxplot}
\title{Voinlin and Boxplot in one plot}
\usage{
pl_vioboxplot(
  meta_data,
  x,
  y,
  fill_colors = NULL,
  box_color = "black",
  box_fill_color = "white"
)
}
\arguments{
\item{meta_data}{data.frame}

\item{x}{a string in colnames of meta_data, to be mapped as x-axis}

\item{y}{a string in colnames of meta_data, to be mapped as y-axis}

\item{fill_colors}{fill colors of violins}

\item{box_color}{border color of boxplot}

\item{box_fill_color}{fill color of boxplot}
}
\value{
a ggplot object
}
\description{
Voinlin and Boxplot in one plot
}
\examples{

}
