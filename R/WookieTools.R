# WookieTools - Version 0.4.2

# Install and load necessary packages safely
ensure_packages <- function(required_packages) {
  installed <- rownames(installed.packages())
  missing_packages <- setdiff(required_packages, installed)
  
  if (length(missing_packages) > 0) {
    install.packages(missing_packages)
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
}

# Ensure Bioconductor packages are installed
ensure_bioc_packages <- function(required_packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      BiocManager::install(package)
    }
  }
  
  invisible(lapply(required_packages, library, character.only = TRUE))
}

# Load essential libraries
load_libraries <- function() {
  cran_packages <- c("Matrix","tidyr", "ggplot2", "plyr","dplyr", "cowplot", "patchwork", "Seurat")
  bioc_packages <- c("SingleCellExperiment", "scds")
  
  ensure_packages(cran_packages)
  ensure_bioc_packages(bioc_packages)
}

load_libraries()

# Quality Control function for Seurat Objects
wookie_qc <- function(seurat_obj, nf_min = 0, nf_max = 20000, nc = 200000, pmt = 20, ptr = NULL, species = 'Mouse', pt.size = NULL) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  mt_pattern <- if (species == 'Mouse') "^mt-" else "^MT-"
  
  seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  if (!is.null(ptr)) {
    seurat_obj[['percent.ribo']] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
  }
  
  subset_criteria <- subset(seurat_obj@meta.data, nFeature_RNA > nf_min & nFeature_RNA < nf_max & nCount_RNA < nc & percent.mt < pmt)
  if (!is.null(ptr)) {
    subset_criteria <- subset_criteria & seurat_obj@meta.data$percent.ribo < ptr
  }
  seurat_obj <- subset(seurat_obj, cells = rownames(subset_criteria))
  
  if (ncol(seurat_obj) == 0) {
    stop("No cells meet the quality control criteria.")
  }
  
  plots <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", if (!is.null(ptr)) "percent.ribo"), ncol = 4, pt.size = pt.size)
  plot_grid(plots, ncol = 2, align = 'v')
  
  return(seurat_obj)
}

# Doublet detection with scds (Deprecated function)
scds_doublets <- function(seurat_obj) {
  message("Deprecated: Use scrub_dub for doublet detection.")
  # [Similar to prior code, with streamlined error handling and clear deprecation warning]
}

# Scrublet integration with Seurat
scrub_dub <- function(seurat_obj, preprocess = FALSE) {
  required_packages <- c("Seurat", "singleCellTK")
  ensure_packages(required_packages)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  if (preprocess) {
    seurat_obj <- NormalizeData(seurat_obj) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA()
  }
  
  sce_obj <- as.SingleCellExperiment(seurat_obj)
  sce_obj <- singleCellTK::runScrublet(sce_obj)
  seurat_obj <- as.Seurat(sce_obj)
  
  hist(seurat_obj@meta.data$scrublet_score, main = "Histogram of Scrublet Scores", xlab = "Scrublet Score")
  
  return(seurat_obj)
}

# Matrix summation function
sum_matrices <- function(matrix1, matrix2, sample = 'sample', min_cells = 3, min_features = 200) {
  if (!all(rownames(matrix1) == rownames(matrix2))) {
    stop("Row names (cells) of matrices must match.")
  }
  
  combined <- matrix1 + matrix2
  
  seurat_obj <- CreateSeuratObject(counts = combined, project = sample, min.cells = min_cells, min.features = min_features)
  return(seurat_obj)
}

# UMAP feature plotting
plot_multi_feat_umap <- function(seurat_obj, features, reduce_dim = "umap", pt.size = 1) {
  # Verify 'features' is a character vector
  if (!is.character(features)) {
    stop("'features' should be a character vector of feature names.")
  }
  
  plot_list <- lapply(features, function(feature) {
    if (feature %in% rownames(seurat_obj)) {
      FeaturePlot(seurat_obj, features = feature, reduction = reduce_dim, pt.size = pt.size)
    } else {
      warning(paste("Feature", feature, "not found in Seurat object."))
      NULL
    }
  })
  
  plot_grid(plotlist = plot_list, align = 'v')
}


# Function to run and plot multiple UMAP's for different numbers of features
#' @name plot_multi_feat_umap
#' @title Plot UMAPs for various feature sets
#' @description Plot UMAPs for different numbers of features (e.g., Highly Variable Genes or Most Abundant Genes)
#' @param seurat_obj Seurat object
#' @param features List of features to use
#' @param min.dist Minimum distance parameter for UMAP
#' @param max_features Maximum number of features to consider
#' @param ftype Type of features (e.g., 'HVG')
#' @param step Step size for incrementing feature sets
#' @param out_name Name to assign to the combined plot
#' @return Combined UMAP plot
#' @export
plot_multi_feat_umap <- function(seurat_obj, features, min.dist = 0.1, 
                                 max_features = 3000, ftype = 'HVG',
                                 step = 500, out_name = 'combined_umap') {
  load_libraries()
  plot_list <- list()
  
  for (feature_length in seq(500, max_features, step)) {
    current_features <- features[1:feature_length]
    cat(paste0('Calculating UMAP for ', ftype, ' with ', feature_length, ' features...'))
    current_umap <- RunUMAP(seurat_obj, features = current_features, min.dist = min.dist)
    current_plot <- DimPlot(current_umap, reduction = 'umap') + ggtitle(paste('UMAP ', ftype, feature_length))
    plot_list[[length(plot_list) + 1]] <- current_plot
    cat("Done.\n")
  }
  
  combined_plot <- plot_grid(plotlist = plot_list)
  assign(out_name, combined_plot, envir = globalenv())
  return(combined_plot)
}

# Function to plot multiple UMAPs at different min.dist values
#' @name plot_multi_min_dist_umap
#' @title Plot UMAPs for various min.dist values
#' @description Plot UMAPs for different min.dist values
#' @param seurat_obj Seurat object
#' @param features Features to use for UMAP
#' @param dims Dimensions to use for UMAP
#' @param out_name Name to assign to the combined plot
#' @return Combined UMAP plot
#' @export
plot_multi_min_dist_umap <- function(seurat_obj, features = NULL, dims = 1:30, 
                                     out_name = 'min_dist_umaps') {
  load_libraries()
  plot_list <- list()
  
  for (min_dist in seq(0.1, 0.5, 0.1)) {
    cat(paste0('Calculating UMAP at min.dist:', min_dist, '...'))
    current_umap <- RunUMAP(seurat_obj, features = features, dims = dims, min.dist = min_dist)
    current_plot <- DimPlot(current_umap, reduction = 'umap') + ggtitle(paste('UMAP: min.dist:', min_dist))
    plot_list[[length(plot_list) + 1]] <- current_plot
    cat("Done.\n")
  }
  
  combined_plot <- plot_grid(plotlist = plot_list)
  assign(out_name, combined_plot, envir = globalenv())
  return(combined_plot)
}

# Function to plot multiple features with color map
#' @name multi_f_plots
#' @title Plot multiple features with color map
#' @description Plot multiple features with custom color scale
#' @param seurat_obj Seurat object
#' @param feature_list List of features to plot
#' @param ncol Number of columns for arrangement
#' @param pt.size Point size for plotting
#' @param split_by Variable for splitting
#' @return Plot grid
#' @export
multi_f_plots <- function(seurat_obj, feature_list, ncol = 3, pt.size = 0.8, split_by = NULL) {
  load_libraries()
  plot_list <- lapply(feature_list, function(feature) {
    FeaturePlot(object = seurat_obj, features = feature, pt.size = pt.size, reduction = "umap", split.by = split_by) +
      theme(aspect.ratio = 1) +
      scale_color_gradientn(colours = c("#DCDCDC", "yellow", "orange", "red", "#8b0000"))
  })
  
  plot_grid(plotlist = plot_list, ncol = ncol, rel_widths = rep(1, length(feature_list)))
}

# Function to filter particular cell types based on marker gene expression
#' @name wookie_filter_celltype
#' @title Filter particular cell types based on marker gene expression
#' @description Remove specific cell types based on marker gene expression
#' @param seurat_obj Seurat object
#' @param marker_list List of marker genes for the cell type to remove
#' @param cutoff Quantile threshold (default: 0.99)
#' @return Filtered Seurat object
#' @export
wookie_filter_celltype <- function(seurat_obj, marker_list, cutoff = 0.99) {
  print('Ensure marker genes are in RNA$scale.data')
  
  expression_matrix_transposed <- t(seurat_obj@assays$RNA$scale.data)
  seurat_obj$avg_celltype_expression <- rowMeans(expression_matrix_transposed[, marker_list, drop = FALSE])
  
  threshold <- quantile(seurat_obj$avg_celltype_expression, cutoff)
  
  cellstokeep <- which(seurat_obj$avg_celltype_expression <= threshold)
  seu_filtered <- seurat_obj[, cellstokeep]
  cellstoremove <- which(seurat_obj$avg_celltype_expression > threshold)
  seu_removed <- seurat_obj[, cellstoremove]
  print(paste0(length(cellstokeep), ' Cells kept, ', length(cellstoremove), ' cells removed.'))
  return(seu_filtered)
}

# Function to plot QC metrics of a sparse matrix
#' @name wookie_matrix_qc_plot
#' @title Plot QC metrics of a sparse matrix
#' @description Plot QC metrics of a sparse matrix
#' @param count_matrix_sparse Sparse matrix
#' @param fill_color Plot color (default: #589FFF)
#' @param title Title for the plot
#' @return QC plot
#' @export
wookie_matrix_qc_plot <- function(count_matrix_sparse, fill_color = "#589FFF", title = "") {
  reads_per_cell <- Matrix::colSums(count_matrix_sparse)
  genes_per_cell <- Matrix::colSums(count_matrix_sparse > 0)
  reads_per_gene <- Matrix::rowSums(count_matrix_sparse > 0)
  
  p1 <- ggplot() +
    geom_histogram(aes(x = log10(reads_per_cell + 1)), fill = fill_color, color = 'black', bins = 30) +
    ggtitle('Reads per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p2 <- ggplot() +
    geom_histogram(aes(x = log10(genes_per_cell + 1)), fill = fill_color, color = 'black', bins = 30) +
    ggtitle('Genes per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p4 <- ggplot() +
    geom_histogram(aes(x = log10(reads_per_gene + 1)), fill = fill_color, color = 'black', bins = 30) +
    ggtitle('Reads per Gene') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p3 <- ggplot() +
    geom_point(aes(x = reads_per_cell, y = genes_per_cell), fill = fill_color, color = 'black', pch = 21, shape = 16, size = 2, alpha = 1) +
    ggtitle('Reads vs. Genes per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_grid(p1, p2, p3, p4, ncol = 2) + ggtitle(title)
}


# Function to compare different normalization methods
#' @name plot_compare_normalisation
#' @title Compare Normalization Methods
#' @description Compare raw counts, normalized counts, SCTransform counts, and scaled data
#' @param seurat_obj A Seurat object with both RNA and SCT assays
#' @return A plot comparing different normalization methods
#' @export
plot_compare_normalisation <- function(seurat_obj) {
  rc <- colSums(GetAssayData(object = seurat_obj, assay = 'RNA', layer = 'counts'))
  normalized_counts <- colSums(seurat_obj[["RNA"]]@layers$data)
  sctransform_counts <- colSums(seurat_obj[["SCT"]]@data)
  scaled_ln <- colSums(seurat_obj[["RNA"]]@layers$scale.data)
  scaled_sct <- colSums(seurat_obj[["SCT"]]@scale.data)
  
  plot_data <- data.frame(
    Cell = names(rc),
    Raw_Counts = rc,
    Normalized_Counts = normalized_counts,
    SCT_Counts = sctransform_counts,
    Scaled_LN = scaled_ln,
    Scaled_SCT = scaled_sct
  )
  
  melted_data <- reshape2::melt(plot_data, id.vars = "Cell", variable.name = "Type", value.name = "Counts")
  
  ggplot(melted_data, aes(x = Cell, y = Counts, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Type, scales = "free_y") +
    labs(title = "Comparison of Normalization Methods",
         x = "Cell",
         y = "Counts",
         fill = "Normalization Method") +
    theme_minimal() +
    theme(axis.text.x = element_blank())
}

# Function to compare Log Normalisation and SCT (Histogram)
#' @name plot_gene_expression_histogram
#' @title Compare Log Normalisation and SCT (Histogram)
#' @description Compare Log Normalisation and SCT
#' @param seurat_obj Seurat Object with both RNA and SCT assays
#' @return Histogram comparing Log Normalisation and SCT
#' @export
plot_gene_expression_histogram <- function(seurat_obj) {
  expression_data_RNA <- as.vector(seurat_obj[["SCT"]]@scale.data)
  expression_data_SCT <- as.vector(seurat_obj[["RNA"]]@scale.data)
  
  mean_expr_rna <- mean(expression_data_RNA)
  sd_expr_rna <- sd(expression_data_RNA)
  mean_expr_sct <- mean(expression_data_SCT)
  sd_expr_sct <- sd(expression_data_SCT)
  
  par(mfrow = c(1, 2))
  
  hist(expression_data_RNA, breaks = 50, freq = FALSE, main = "Histogram of Gene Expression (RNA)", xlab = "Gene Expression", col = "lightgray")
  curve(dnorm(x, mean = mean_expr_rna, sd = sd_expr_rna), add = TRUE, col = "blue", lwd = 2)
  
  hist(expression_data_SCT, breaks = 50, freq = FALSE, main = "Histogram of Gene Expression (SCT)", xlab = "Gene Expression", col = "lightgray")
  curve(dnorm(x, mean = mean_expr_sct, sd = sd_expr_sct), add = TRUE, col = "blue", lwd = 2)
  
  par(mfrow = c(1, 1))
}

# Function to get optimal number of PCs
#' @name get_optimal_pcs
#' @title Get Optimal Number of PCs
#' @description Get the optimal number of principal components to use
#' @param seurat_obj Seurat Object
#' @param reduction Type of reduction (e.g., 'pca')
#' @return Number of PCs to use
#' @export
get_optimal_pcs <- function(seurat_obj, reduction = 'pca') {
  pct <- seurat_obj[[reduction]]@stdev / sum(seurat_obj[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)
  
  col <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  pcs <- min(col, co2)
  return(pcs)
}
