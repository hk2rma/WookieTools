% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WookieTools.R
\name{wookieJaccardPlot}
\alias{wookieJaccardPlot}
\title{Function to plot Jaccard similarity index}
\usage{
wookieJaccardPlot(
  seurat_obj,
  clusters = "seurat_clusters",
  logfc.threshold = 0.95,
  min.pct = 0.25,
  test.use = "wilcox",
  fdr.threshold = 0.05,
  p_val_bonf.threshold = 0.05,
  avg.log2fc.threshold = 0.25,
  title = "Wookie_Jaccard",
  saveplot = FALSE,
  height = 10,
  width = 10,
  dpi = 700,
  units = "cm",
  limitsize = FALSE,
  silentwookie = FALSE
)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{clusters}{Name of the column in the Seurat object containing cluster assignments, default is 'seurat_clusters'}

\item{logfc.threshold}{Log fold-change threshold, default is 0.95}

\item{min.pct}{Minimum percentage of cells expressing a gene, default is 0.25}

\item{test.use}{Statistical test to use for finding differentially expressed genes, default is "wilcox"}

\item{fdr.threshold}{False discovery rate threshold, default is 0.05}

\item{p_val_bonf.threshold}{Bonferroni-corrected p-value threshold, default is 0.05}

\item{avg.log2fc.threshold}{Average log2 fold-change threshold, default is 0.25}

\item{title}{Title for the plot, default is 'Wookie_Jaccard'}

\item{saveplot}{Boolean indicating whether to save the plot, default is FALSE}

\item{height}{Height of the plot in the specified units, default is 10}

\item{width}{Width of the plot in the specified units, default is 10}

\item{dpi}{Resolution of the saved plot, default is 700}

\item{units}{Units for the plot dimensions, default is 'cm'}

\item{limitsize}{Boolean indicating whether to limit the size of the plot, default is FALSE}

\item{silentwookie}{stop wookie from printing puns, default is FALSE}
}
\value{
plot Jaccard similarity index
}
\description{
Function to plot Jaccard similarity index
}
