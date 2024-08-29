Wookie Tools
Version 0.1.24A

Handy Tools for Single Cell RNA-Seq Analysis

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

Note:
To Install the development version, which has the latest features and fixes:
```{r, install dev branch}
if (!requireNamespace("WookieTools", quietly = TRUE)) {
  devtools::install_github("hk2rma/WookieTools@babywookie")
}
```
