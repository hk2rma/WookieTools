% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WookieTools.R
\name{wookiePlotQCHistogram}
\alias{wookiePlotQCHistogram}
\title{Plot Count v Features with percentage mito}
\usage{
wookiePlotQCHistogram(
  seurat_obj,
  title = "Histogram",
  fi = 0,
  ci = 0,
  mi = 0,
  ri = 0,
  bins = 300,
  log.counts = FALSE,
  silentwookie = FALSE
)
}
\arguments{
\item{seurat_obj}{Seurat Object}

\item{title}{title of the plot}

\item{fi}{features threshold}

\item{ci}{counts threshold}

\item{bins}{bins}

\item{log.counts}{Boolean to transform nCount to log scale}

\item{silentwookie}{stop wookie from printing puns, default is FALSE}
}
\value{
plot
}
\description{
Get histograms of nFeature and nCounts and possible thresholds
}
