% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WookieTools.R
\name{wookiePlotQCMetrics}
\alias{wookiePlotQCMetrics}
\title{Plot Count v Features with percentage mito}
\usage{
wookiePlotQCMetrics(
  seu,
  limits = c(50, 20000),
  pmt = 15,
  size = 0.2,
  alpha = 0.9,
  silentwookie = FALSE
)
}
\arguments{
\item{seu}{Seurat Object}

\item{limits}{milimits for ggplot}

\item{pmt}{mito cutoff}

\item{size}{size geom}

\item{alpha}{geom alpha}

\item{silentwookie}{stop wookie from printing puns, default is FALSE}
}
\value{
plot
}
\description{
Plot Count v Features with percentage mito
}
