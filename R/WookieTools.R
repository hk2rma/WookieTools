# WookieTools - Version 0.5.3

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

# Seurat Object Quality Control function
#' @name wookieqc
#' @title Seurat Object Quality Control function
#' @description Run Iteratively for the necessary QC...
#' @param matrix Seurat object ...
#' @param nf_min Minimum number of features ...
#' @param nf_max Maximum number of features ...
#' @param nc Maximum number of counts ...
#' @param pmt Percentage of mitochondrial genes ...
#' @param ptr Percentage of ribosomal genes ...
#' @param group Grouping variable ...
#' @param species species in dataset Mouse or Human only ...
#' @param colors Colors for facetting ...
#' @param pt.size data points in violin plot
#' @return Seurat object after quality control
#' @export
wookie_qc <- function(seurat_obj, nf_min = 0, nf_max = 20000,
                      nc_max = 200000,nc_min = 0, pmt = 20,
                      ptr = NULL, species = 'Mouse', 
                      pt.size = NULL,legend = TRUE) {
 
  
   if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  mt_pattern <- if (species == 'Mouse') "^mt-" else "^MT-"
  
  seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  if (!is.null(ptr)) {
    seurat_obj[['percent.ribo']] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
  }
  
  subset_criteria <- subset(seurat_obj@meta.data, nFeature_RNA > nf_min &
                              nFeature_RNA < nf_max & nCount_RNA < nc_max &
                              nCount_RNA > nc_min &percent.mt < pmt)
  
  if (!is.null(ptr)) {
    ribo_indices <- which(seurat_obj@meta.data$percent.ribo < ptr)
    subset_criteria <- subset_criteria[ribo_indices, ]
  }
  
  seurat_obj <- subset(seurat_obj, cells = rownames(subset_criteria))
  
  if (ncol(seurat_obj) == 0) {
    stop("No cells meet the quality control criteria.")
  }
  
  
  # Visualizations
  vl_plot <- VlnPlot(seurat_obj,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                  if (!is.null(ptr)) "percent.ribo"),
                                  ncol = 4, pt.size = pt.size)
  
  plot1 <- if (legend) {
    FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                   feature2 = "nFeature_RNA")
  } else {
    FeatureScatter(seurat_obj, feature1 = "nCount_RNA",
                   feature2 = "nFeature_RNA") + NoLegend()
  }
  
  plot2 <- if (legend) {
    FeatureScatter(seurat_obj, feature1 = "percent.mt",
                   feature2 = "nFeature_RNA")
  } else {
    FeatureScatter(seurat_obj, feature1 = "percent.mt",
                   feature2 = "nFeature_RNA") + NoLegend()
  }
  
  plots_list <- list(vl_plot, plot1, plot2)
  
  if (!is.null(ptr)) {
    plot3 <- if (legend) {
      FeatureScatter(seurat_obj, feature1 = 'percent.ribo',
                     feature2 = 'nFeature_RNA')
    }        else {
      FeatureScatter(seurat_obj, feature1 = 'percent.ribo',
                     feature2 = 'nFeature_RNA') + NoLegend()
    }
    plots_list <- c(plots_list, list(plot3))
  }
  
  # Combine all plots into a single plot
  combined_plot <- cowplot::plot_grid(plotlist = plots_list,
                                      ncol = 2, align = 'v')
  
  # Display the combined plot
  
  print(combined_plot)
  return(seurat_obj)
  
}

# Doublet finder using Scrublet
#' @name wookie_scrub
#' @title Scrublet for Seurat ...
#' @description Run Scrublet on a Seurat Object ...
#' @param seu_obj Seurat object ...
#' @param preprocess specify if the object has been preprocces and runPCA was done  ...
#' @return Seurat object with scrublet scores and call
#' @export
wookie_scrub <- function(seu_obj, preprocess = FALSE) {
  # Load required packages
  required_packages <- c("Seurat", "SingleCellExperiment", "singleCellTK")
  for (package in required_packages) {
    if (!require(package, character.only = TRUE)) {
      install.packages(package)
      library(package, character.only = TRUE)
    }
  }
  
  # Input validation
  if (!inherits(seu_obj, "Seurat")) {
    stop("The 'seu_obj' parameter must be a valid Seurat object.")
  }
  
  # Preprocess the Seurat object if required
  if (preprocess) {
    message("Preprocessing the Seurat object...")
    seu_obj <- Seurat::NormalizeData(seu_obj)
    seu_obj <- Seurat::FindVariableFeatures(seu_obj)
    seu_obj <- Seurat::ScaleData(seu_obj)
    seu_obj <- Seurat::RunPCA(seu_obj)
  }
  
  # Convert Seurat object to SingleCellExperiment
  message("Running Scrublet...")
  sce_obj <- as.SingleCellExperiment(seu_obj)
  
  # Run Scrublet
  sce_scrub <- singleCellTK::runScrublet(sce_obj)
  
  # Convert SingleCellExperiment back to Seurat
  seu_obj_scrubbed <- as.Seurat(sce_scrub)
  
  # Extract Scrublet scores and cell type calls
  scrub_scores <- seu_obj_scrubbed@meta.data$scrublet_score
  scrub_type <- seu_obj_scrubbed@meta.data$scrublet_call
  
  # Add Scrublet scores and cell type calls to the original Seurat object
  message("Adding Scrublet results to the Seurat object...")
  seu_obj@meta.data$scrublet_score <- scrub_scores
  seu_obj@meta.data$scrublet_call <- scrub_type
  
  # Plot histogram of Scrublet scores
  hist(scrub_scores, main = "Histogram of Scrublet Scores",
       xlab = "Scrublet Score")
  
  return(seu_obj)
}

# Function to sum the counts of two matrices containing the same cells
# Input: Count Matrices | Output: Seurat Object
#' @name wookie_matrix_sum
#' @title Matrix Sum function
#' @description Merge two count matrices, where the cells are the same, to obtain a single seurat object with the counts from two matrices summed for each cell ...
#' @param matrix1 count matrix ...
#' @param matrix2 count matrix ...
#' @param sample sample/library name ...
#' @param min_cells minimum number cells a gene is found in ...
#' @param min_features minimum number of features found in a cell ...
#' @return Summed and merged Seurat object
#' @export
wookie_matrix_sum <- function(matrix1, matrix2, sample = 'sample',
                              min_cells = 3, min_features = 200) {
  load_libraries()
  
  # Check if row names are identical
  if (!identical(rownames(matrix1), rownames(matrix2))) {
    stop('Error: Row names are not identical.')
  }
  # Check if Column names are identical
  if (!identical(rownames(matrix1), rownames(matrix2))) {
    print(paste0('Warning: Column names are not identical.'))
  }
  # Identify columns not common to both matrices
  extra_cols_matrix1 <- setdiff(colnames(matrix1), colnames(matrix2))
  extra_cols_matrix2 <- setdiff(colnames(matrix2), colnames(matrix1))
  
  common_rows <- intersect(rownames(matrix1), rownames(matrix2))
  common_cols <- intersect(colnames(matrix1), colnames(matrix2))
  
  # Subset matrices to common rows and columns
  matrix1_common <- matrix1[which(rownames(matrix1) %in% common_rows),
                            which(colnames(matrix1) %in% common_cols)]
  matrix2_common <- matrix2[which(rownames(matrix2) %in% common_rows), 
                            which(colnames(matrix2) %in% common_cols)]
  
  # Sum the matrices
  result_matrix <- matrix1_common + matrix2_common
  
  # Concatenate extra columns to the right of the result matrix
  if (length(extra_cols_matrix1) > 0) {
    matrix1_uncommon <- matrix1[common_rows,extra_cols_matrix1]
    result_matrix <- cbind(result_matrix, matrix1_uncommon)
  }
  
  if (length(extra_cols_matrix2) > 0) {
    matrix2_uncommon <- matrix2[common_rows,extra_cols_matrix2]
    result_matrix <- cbind(result_matrix, matrix2_uncommon)
  }
  
  original_col_order <- c(colnames(matrix1), extra_cols_matrix2)
  result_matrix <- result_matrix[, original_col_order]
  
  # Create Seurat object
  seu_obj <- CreateSeuratObject(result_matrix, min.cells = min_cells,
                                min.features = min_features, project = sample)
  
  return(seu_obj)
}

# Function to run and plot multiple UMAP's for different numbers of a feature(Highly variable genes or Most Abundant genes)
# Features must be obtained and given as input
#' @name wookie_multifeatureumap
#' @title Plot UMAPs to test features
#' @description plot multiple UMAP's for different numbers of a feature i.e Highly variable genes or Most Abundant genes ...
#' @return plot saved to global environment
#' @export
wookie_multifeatureumap <- function(object = seu_obj, features = features,
                                    min.dist = 0.1, max_features = 3000,
                                    ftype='HVG',step = 500,
                                    out_name = 'combined_umap') {
  load_libraries()
  plot_list <- list()
  
  for (feature_length in seq(500, max_features, step)) {
    current_features <- features[1:feature_length]
    cat(paste0('Calculating UMAP at ',ftype,':',feature_length))
    current_umap <- RunUMAP(object, features = current_features, min.dist = min.dist)
    current_plot <- DimPlot(current_umap, reduction = 'umap') + 
                    ggtitle(paste('UMAP ',ftype, feature_length))
    plot_list[[length(plot_list) + 1]] <- current_plot
    cat(paste0('UMAP done for ',ftype,':',feature_length))
  }
  
  # Combine plots into a grid
  combined_plot <- plot_grid(plotlist = plot_list)
  
  # Assign the combined plot to a variable in the global environment
  assign(out_name, combined_plot, envir = globalenv())
  
  # Return the combined plot
  return(combined_plot)
}


# Function to plot multiple UMAPs at different min.dist values
#' @name wookie_Mindist
#' @title Plot UMAPs for various min.dist values
#' @description Plot UMAPs for different min.dist values
#' @param seurat_obj Seurat object
#' @param features Features to use for UMAP
#' @param dims Dimensions to use for UMAP
#' @param out_name Name to assign to the combined plot
#' @return Combined UMAP plot
#' @export
wookie_Mindist <- function(seurat_obj, features = NULL, dims = 1:30, 
                                     out_name = 'min_dist_umaps') {
  load_libraries()
  plot_list <- list()
  
  for (min_dist in seq(0.1, 0.5, 0.1)) {
    cat(paste0('Calculating UMAP at min.dist:', min_dist, '...'))
    current_umap <- RunUMAP(seurat_obj, features = features,
                            dims = dims, min.dist = min_dist)
    current_plot <- DimPlot(current_umap, reduction = 'umap') + 
                          ggtitle(paste('UMAP: min.dist:', min_dist))
    plot_list[[length(plot_list) + 1]] <- current_plot
    cat("Done.\n")
  }
  
  combined_plot <- plot_grid(plotlist = plot_list)
  assign(out_name, combined_plot, envir = globalenv())
  return(combined_plot)
}

# Function to plot multiple features with color map
#' @name wookie_featureplot
#' @title Plot multiple features with color map
#' @description Plot multiple features with custom color scale
#' @param seurat_obj Seurat object
#' @param feature_list List of features to plot
#' @param ncol Number of columns for arrangement
#' @param pt.size Point size for plotting
#' @param split_by Variable for splitting
#' @return Plot grid
#' @export
wookie_featureplot <- function(seuratObject, featureList, ncol = 3,
                               pt.size = 0.8,split_by = NULL) {
  load_libraries()
  plotList <- lapply(featureList, function(feature) {
    FeaturePlot(object = seuratObject, features = feature, 
                pt.size = pt.size, reduction = "umap",
                split.by = split_by) +
      theme(aspect.ratio = 1) +
      scale_color_gradientn(colours = c("#DCDCDC", "yellow", "orange", "red", "#8b0000"))
  })
  
  plotGrid <- plot_grid(plotlist = plotList, ncol = ncol, 
                        rel_widths = rep(1, length(featureList)))
  
  return(plotGrid)
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
  seurat_obj$avg_celltype_expression <- rowMeans(expression_matrix_transposed[, marker_list,
                                                                              drop = FALSE])
  
  threshold <- quantile(seurat_obj$avg_celltype_expression, cutoff)
  
  cellstokeep <- which(seurat_obj$avg_celltype_expression <= threshold)
  seu_filtered <- seurat_obj[, cellstokeep]
  cellstoremove <- which(seurat_obj$avg_celltype_expression > threshold)
  seu_removed <- seurat_obj[, cellstoremove]
  print(paste0(length(cellstokeep), ' Cells kept, ', 
               length(cellstoremove), ' cells removed.'))
  return(seu_filtered)
}

# Function to plot QC metrics of a sparse matrix
#' @name wookie_matrix_qc
#' @title Plot QC metrics of a sparse matrix
#' @description Plot QC metrics of a sparse matrix
#' @param count_matrix_sparse Sparse matrix
#' @param fill_color Plot color (default: #589FFF)
#' @param title Title for the plot
#' @return QC plot
#' @export
wookie_matrix_qc <- function(count_matrix_sparse, fill_color = "#589FFF", title = "") {
  reads_per_cell <- Matrix::colSums(count_matrix_sparse)
  genes_per_cell <- Matrix::colSums(count_matrix_sparse > 0)
  reads_per_gene <- Matrix::rowSums(count_matrix_sparse > 0)
  
  p1 <- ggplot() +
    geom_histogram(aes(x = log10(reads_per_cell + 1)), fill = fill_color,
                   color = 'black', bins = 30) +
    ggtitle('Reads per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p2 <- ggplot() +
    geom_histogram(aes(x = log10(genes_per_cell + 1)), fill = fill_color,
                   color = 'black', bins = 30) +
    ggtitle('Genes per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p4 <- ggplot() +
    geom_histogram(aes(x = log10(reads_per_gene + 1)), fill = fill_color,
                   color = 'black', bins = 30) +
    ggtitle('Reads per Gene') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p3 <- ggplot() +
    geom_point(aes(x = reads_per_cell, y = genes_per_cell), fill = fill_color,
               color = 'black', pch = 21, shape = 16, size = 2, alpha = 1) +
    ggtitle('Reads vs. Genes per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_grid(p1, p2, p3, p4, ncol = 2) + ggtitle(title)
}

# Function to compare Log Normalisation and SCT (Histogram)
#' @name wookie_ge_histogram
#' @title Compare Log Normalisation and SCT (Histogram)
#' @description Compare Log Normalisation and SCT
#' @param seurat_obj Seurat Object with both RNA and SCT assays
#' @return Histogram comparing Log Normalisation and SCT
#' @export
wookie_ge_histogram <- function(seurat_obj) {
  expression_data_RNA <- as.vector(seurat_obj[["SCT"]]@scale.data)
  expression_data_SCT <- as.vector(seurat_obj[["RNA"]]@scale.data)
  
  mean_expr_rna <- mean(expression_data_RNA)
  sd_expr_rna <- sd(expression_data_RNA)
  mean_expr_sct <- mean(expression_data_SCT)
  sd_expr_sct <- sd(expression_data_SCT)
  
  par(mfrow = c(1, 2))
  
  hist(expression_data_RNA, breaks = 50, freq = FALSE,
       main = "Histogram of Gene Expression (RNA)",
       xlab = "Gene Expression", col = "lightgray")
  curve(dnorm(x, mean = mean_expr_rna, sd = sd_expr_rna),
        add = TRUE, col = "blue", lwd = 2)
  
  hist(expression_data_SCT, breaks = 50, freq = FALSE,
       main = "Histogram of Gene Expression (SCT)", 
       xlab = "Gene Expression", col = "lightgray")
  curve(dnorm(x, mean = mean_expr_sct, sd = sd_expr_sct),
        add = TRUE, col = "blue", lwd = 2)
  
  par(mfrow = c(1, 1))
}

# Function to get optimal number of PCs
#' @name wookie_get_pc 
#' @title Get Optimal Number of PCs
#' @description Get the optimal number of principal components to use
#' @param seurat_obj Seurat Object
#' @param reduction Type of reduction (e.g., 'pca')
#' @return Number of PCs to use
#' @export
wookie_get_pc <- function(seurat_obj, reduction = 'pca') {
  pct <- seurat_obj[[reduction]]@stdev / sum(seurat_obj[[reduction]]@stdev) * 100
  cumu <- cumsum(pct)
  
  col <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
              decreasing = TRUE)[1] + 1
  pcs <- min(col, co2)
  return(pcs)
}


# Function to create a plot for a given Seurat object
#' @name wookie_fc_hist 
#' @title histograms of nFeature and nCounts
#' @description Get histograms of nFeature and nCounts and possible thresholds
#' @param seurat_obj Seurat Object
#' @param title title of the plot
#' @param fi features threshold
#' @param ci counts threshold 
#' @return plot
#' @export
wookie_fc_hist <- function(seurat_obj, title = 'Histogram', fi = 0, ci = 0) {
  # Extract data
  data <- FetchData(seurat_obj, vars = c("nFeature_RNA", "nCount_RNA"))
  
  # Create histograms using ggplot2
  p1 <- ggplot(data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100, fill = "#06125F") +
    geom_vline(xintercept = fi, color = "#FF0909", linetype = "dashed") +
    ggtitle(paste(title, "- Features"))
  
  p2 <- ggplot(data, aes(x = nCount_RNA)) +
    geom_histogram(bins = 100, fill = "#06125F") +
    geom_vline(xintercept = ci, color = "#FF0909", linetype = "dashed") +
    ggtitle(paste(title, "- Counts"))
  
  # Return combined plot for each seurat object
  return(p1 + p2 + plot_layout(ncol = 2))
}