Wookie Tools
Version 0.5.1

Hand Tools for Single Cell RNA-Seq Analysis

Installation: 
```{r, install Wookie Tools}

# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("WookieTools", quietly = TRUE)) {
  devtools::install_github("hk2rma/WookieTools")
}
library(Seurat)
library(WookieTools)

```