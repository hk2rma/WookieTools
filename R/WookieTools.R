# WookieTools - Version 0.6.5.1

# Seurat Object Quality Control function
#' @name wookieqc
#' @import Seurat 
#' @import ggplot2
#' @import cowplot
#' @title Seurat Object Quality Control function
#' @description Run Iteratively for the necessary QC...
#' @param matrix Seurat object ...
#' @param nf_min Minimum number of features ...
#' @param nf_max Maximum number of features ...
#' @param nc Maximum number of counts ...
#' @param pmt Percentage of mitochondrial genes ...
#' @param ptr_max maximum percentage of ribosomal genes ...
#' @param ptr_min minimum percentage of ribosomal genes ...
#' @param group Grouping variable ...
#' @param species species in dataset Mouse or Human only ...
#' @param colors Colors for facetting ...
#' @param pt.size data points in violin plot
#' @return Seurat object after quality control
#' @export
wookie_qc <- function(seurat_obj, nf_min = 0, nf_max = 20000,
                      nc_max = 200000, nc_min = 0, pmt = 20,
                      ptr_max = NULL, ptr_min = NULL, species = 'Mouse',
                      pt.size = NULL, legend = TRUE) {
  
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  
  mt_pattern <- if (species == 'Mouse') "^mt-" else "^MT-"
  
  seurat_obj[['percent.mt']] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  seurat_obj[['percent.ribo']] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")
  
  seurat_obj <- subset(seurat_obj, nFeature_RNA > nf_min &
                         nFeature_RNA < nf_max & nCount_RNA < nc_max &
                         nCount_RNA > nc_min & percent.mt < pmt)
  
  if (!is.null(ptr_min)) {
    seurat_obj <- subset(seurat_obj, percent.ribo > ptr_min) 
  }
  if (!is.null(ptr_max) ) {
    seurat_obj <- subset(seurat_obj, percent.ribo < ptr_max) 
  }  
  
  if (ncol(seurat_obj) == 0) {
    stop("No cells meet the quality control criteria.")
  }
  
  
  # Visualizations
  vl_plot <- VlnPlot(seurat_obj,
                     features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                                  if (!is.null(ptr_max) || !is.null(ptr_min)) "percent.ribo"),
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
  
  if (!is.null(ptr_max) || !is.null(ptr_min)) {
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
  wookieSay()
  return(seurat_obj)
  
  
}

# Doublet finder using Scrublet
#' @name wookie_scrub
#' @title Scrublet for Seurat ...
#' @import Seurat 
#' @import ggplot2
#' @importFrom singleCellTK runScrublet
#' @import cowplot
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
  sce_scrub <- runScrublet(sce_obj)
  
  # Convert SingleCellExperiment back to Seurat
  seu_obj_scrubbed <- as.Seurat(sce_scrub)
  
  # Extract Scrublet scores and cell type calls
  scrub_scores <- seu_obj_scrubbed@meta.data$scrublet_score
  scrub_type <- seu_obj_scrubbed@meta.data$scrublet_call
  
  # Add Scrublet scores and cell type calls to the original Seurat object
  message("Adding Scrublet results to the Seurat object...")
  seu_obj@meta.data$scrublet_score <- scrub_scores
  seu_obj@meta.data$scrublet_call <- scrub_type
  wookieSay()
  return(seu_obj)
}

# Function to sum the counts of two matrices containing the same cells
# Input: Count Matrices | Output: Seurat Object
#' @name wookie_matrix_sum
#' @title Matrix Sum function
#' @import Matrix
#' @import dplyr
#' @import tidyr
#' @import patchwork
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
  wookieSay()
  return(seu_obj)
}

# Function to run and plot multiple UMAP's for different numbers of a feature(Highly variable genes or Most Abundant genes)
# Features must be obtained and given as input
#' @name wookie_multifeatureumap
#' @title Plot UMAPs to test features
#' @import Seurat
#' @import ggplot2
#' @import cowplot
#' @description plot multiple UMAP's for different numbers of a feature i.e Highly variable genes or Most Abundant genes ...
#' @return plot saved to global environment
#' @export
wookie_multifeatureumap <- function(object = seu_obj, features = features,
                                    min.dist = 0.1, max_features = 3000,
                                    ftype='HVG',step = 500,
                                    out_name = 'combined_umap') {
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
  wookieSay()
  # Return the combined plot
  return(combined_plot)
}


# Function to plot multiple UMAPs at different min.dist values
#' @name wookie_Mindist
#' @title Plot UMAPs for various min.dist values
#' @import Seurat
#' @import cowplot
#' @import ggplot2
#' @description Plot UMAPs for different min.dist values
#' @param seurat_obj Seurat object
#' @param features Features to use for UMAP
#' @param dims Dimensions to use for UMAP
#' @param out_name Name to assign to the combined plot
#' @return Combined UMAP plot
#' @export
wookie_Mindist <- function(seurat_obj, features = NULL, dims = 1:30, 
                                     out_name = 'min_dist_umaps') {
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
  wookieSay()
  return(combined_plot)
}

# Function to plot multiple features with color map
#' @name wookie_featureplot
#' @title Plot multiple features with color map
#' @import Seurat
#' @import ggplot2
#' @import cowplot
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
  plotList <- lapply(featureList, function(feature) {
    FeaturePlot(object = seuratObject, features = feature, 
                pt.size = pt.size, reduction = "umap",
                split.by = split_by) +
      theme(aspect.ratio = 1) +
      scale_color_gradientn(colours = c("#DCDCDC", "yellow", "orange", "red", "#8b0000"))
  })
  
  plotGrid <- plot_grid(plotlist = plotList, ncol = ncol, 
                        rel_widths = rep(1, length(featureList)))
  wookieSay()
  return(plotGrid)
}




# Function to filter particular cell types based on marker gene expression
#' @name wookie_filter_celltype
#' @title Filter particular cell types based on marker gene expression
#' @import Seurat
#' @import patchwork
#' @import tidyr
#' @import Matrix
#' @import dplyr
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
  wookieSay()
  return(seu_filtered)
}

# Function to plot QC metrics of a sparse matrix
#' @name wookie_matrix_qc
#' @title Plot QC metrics of a sparse matrix
#' @import Seurat
#' @import patchwork
#' @import ggplot2
#' @import cowplot
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
  wookieSay()
}

# Function to compare Log Normalisation and SCT (Histogram)
#' @name wookie_ge_histogram
#' @title Compare Log Normalisation and SCT (Histogram)
#' @import Seurat
#' @import patchwork
#' @import ggplot2
#' @import cowplot
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
  wookieSay()
}

# Function to get optimal number of PCs
#' @name wookie_get_pc 
#' @title Get Optimal Number of PCs
#' @import Seurat
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
  cat(paste0('No. of Optimal PCs: ',pcs))
  wookieSay()
}


# Function to create a plot for a given Seurat object
#' @name wookie_fc_hist 
#' @title histograms of nFeature and nCounts
#' @import Seurat
#' @import patchwork
#' @import ggplot2
#' @import cowplot
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
  
  wookieSay()
  # Return combined plot for each seurat object
  return(p1 + p2 + plot_layout(ncol = 2))
  
}

# Function to plot multiple features with color map
#' @name wookie_dotplot
#' @title Plot dotplots with color map
#' @import Seurat
#' @import ggplot2
#' @import cowplot
#' @description dotplot
#' @param seurat_obj Seurat object
#' @param assay default RNA
#' @param feature_list List of markers
#' @param tag plot title
#' @param scale.by scale.by
#' @return Pdotplot
#' @export
wookie_dotplot <- function(seurat_obj, feature_list , assay = 'RNA',scale.by = 'size', tag = 'DotPlot'){
  dp <- DotPlot(seurat_obj, features = feature_list, assay = assay,scale.by = scale.by) +
    coord_flip()+
    scale_color_gradientn(colours = c("#DCDCDC", "yellow", "orange", "red", "#8b0000"))+ 
    ggtitle(tag)
  wookieSay()
  return(dp)
}


# Function to plot Cluster relationship Tree
#' @name wookie_clustertree
#' @title Function to plot Cluster relationship Tree
#' @import Seurat
#' @import ggplot2
#' @description Function to plot Cluster relationship Tree
#' @param seurat_obj Seurat object
#' @param dims dimensions to use
#' @param reorder default is FALSE
#' @param assay default is 'RNA'
#' @param reduction default is 'pca'
#' @param slot default is 'data'
#' @param reorder.numeric default is FALSE
#' @param saveplot save a jpeg image to working directory
#' @param dpi default is 1080
#' @param height default is 10
#' @param width default is 10
#' @param units deafults is 'cm'
#' @return p
#' @export
wookie_clustertree <- function(seurat_obj, dims = NULL, features = NULL, reorder = FALSE ,
                               reorder.numeric = FALSE, saveplot = FALSE,assay = 'RNA',
                               reduction = 'pca', slot = 'data',
                               dpi = 1080, height = 10, width = 10, units = 'cm'){
  seurat <- BuildClusterTree(
    seurat_obj,
    dims = dims,
    assay = assay,
    reduction = reduction,
    slot = slot,
    features = features,
    reorder = reorder,
    reorder.numeric = reorder.numeric
  )
  
  tree <- seurat@tools$BuildClusterTree
  tree$tip.label <- paste0("Cluster ", tree$tip.label)
  num_tips <- length(tree$tip.label)
  p <- ggtree::ggtree(tree, aes(x, y)) +
    scale_y_reverse() +
    ggtree::geom_tree() +
    ggtree::theme_tree() +
    ggtree::geom_tiplab(offset = 1) +
    ggtree::geom_tippoint(color = rainbow(num_tips), shape = 16, size = 5) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))
  if(saveplot == TRUE){
    ggsave('Cluster_tree.jpeg',p,width = width,height = height,dpi = dpi,units = units)
  }
  wookieSay()
  return(p)
  
}


# Function to plot expression per cluster
#' @name wookie_pcePlot
#' @title Function to plot expression per cluster
#' @import Seurat
#' @import ggplot2
#' @import gghalves
#' @description Function to plot Cluster relationship Tree
#' @param seurat_obj Seurat object
#' @param seurat_clusters default is 'seurat_clusters'
#' @param saveplot save a jpeg image to working directory
#' @param dpi default is 1080
#' @param height default is 10
#' @param width default is 10
#' @param units deafults is 'cm'
#' @return pce_plot
#' @export
wookie_pcePlot <- function(seurat, seurat_clusters= 'seurat_clusters' ,
  saveplot = FALSE,dpi = 1080, height = 10, width = 10, units = 'cm'){
  
  temp_labels <- seurat@meta.data %>%
    group_by(seurat_clusters) %>%
    tally()

  p1 <- ggplot() +
    geom_half_violin(
      data = seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
      side = 'l', show.legend = FALSE, trim = FALSE
      ) +
    geom_half_boxplot(
      data = seurat@meta.data, aes(seurat_clusters, nCount_RNA, fill = seurat_clusters),
      side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
      ) +
    geom_text(
      data = temp_labels,
      aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
      ) +
    scale_color_manual(values = custom_colors$discrete) +
    scale_fill_manual(values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of transcripts', labels = scales::comma, expand = c(0.08, 0)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank()
    )

  p2 <- ggplot() +
    geom_half_violin(
      data = seurat@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
      side = 'l', show.legend = FALSE, trim = FALSE
      ) +
    geom_half_boxplot(
      data = seurat@meta.data, aes(seurat_clusters, nFeature_RNA, fill = seurat_clusters),
      side = 'r', outlier.color = NA, width = 0.4, show.legend = FALSE
      ) +
    geom_text(
      data = temp_labels,
      aes(x = seurat_clusters, y = -Inf, label = paste0('n = ', format(n, big.mark = ',', trim = TRUE)), vjust = -1),
      color = 'black', size = 2.8
      ) +
    scale_color_manual(values = custom_colors$discrete) +
    scale_fill_manual(values = custom_colors$discrete) +
    scale_y_continuous(name = 'Number of expressed genes', labels = scales::comma, expand = c(0.08, 0)) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.title.x = element_blank()
    )
   
    pce_plot <- p1 + p2 + plot_layout(ncol = 1)
    if(saveplot == TRUE){
    ggsave(
      'plots/ncount_nfeature_by_cluster.jpeg',
      p1 + p2 + plot_layout(ncol = 1),
      height = height, width = width,
      dpi = dpi,units = units)
    }
    wookieSay()
    return(pce_plot)
    
}


# Function to plot silhoutte scores
#' @name wookie_silhouttePlot
#' @title Function to plot silhoutte scores for each cluster
#' @import Seurat
#' @import ggplot2
#' @import cluster
#' @description Function to plot silhoutte scores for each cluster
#' @param seurat Seurat object
#' @param cluster default is 'seurat_clusters'
#' @param dims dimension, default is 1:30
#' @param reduction default is 'pca'
#' @return silhoutte plot
#' @export
wookie_silhouttePlot <- function(seurat,cluster = 'seurat_clusters',dims = 1:30,reduction = 'pca'){
  seurat$seurat_clusters <- seurat[[cluster]]
  distance_matrix <- dist(Embeddings(seurat[[reduction]])[, dims])
  clusters <- seurat@meta.data$seurat_clusters
  silhouette <- silhouette(as.numeric(clusters), dist = distance_matrix)
  seurat@meta.data$silhouette_score <- silhouette[,3]
  
  mean_silhouette_score <- mean(seurat@meta.data$silhouette_score)
  
  silhoutte_plot <- seurat@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(seurat_clusters,-silhouette_score) %>%
    mutate(barcode = factor(barcode, levels = barcode)) %>%
    ggplot() +
    geom_col(aes(barcode, silhouette_score, fill = seurat_clusters), show.legend = FALSE) +
    geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette score') +
    scale_fill_manual(values = custom_colors$discrete) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )+
    geom_text(
      data = . %>% group_by(seurat_clusters) %>% slice(1),
      aes(x = barcode, y = silhouette_score, label = seurat_clusters),
      hjust = -0.2,
      vjust = 0.5,
      angle = 90,
      size = 3
    )
  wookieSay()
  return(silhoutte_plot)
  
}


# Function to plot cluster Similarity
#' @name wookie_clusterSimilarityPlot
#' @title Function to plot cluster Similarity
#' @import Seurat
#' @import scran
#' @import bluster
#' @import ggplot2
#' @import cluster
#' @description Function to plot cluster Similarity
#' @param seurat_obj Seurat object
#' @param clusters default is 'seurat_clusters'
#' @param dims dimension, default is 1:30
#' @param reduction default is 'PCA'
#' @return cluster Similarity heatmap
#' @export
wookie_clusterSimilarityPlot <- function(seurat_obj,dims = 1:30,clusters = 'seurat_clusters',reduction = 'PCA'){
  seurat_obj$seurat_clusters <- seurat_obj[[clusters]]
  sce <- as.SingleCellExperiment(seurat_obj)
  reducedDim(sce, 'PCA_sub') <- reducedDim(sce, 'PCA')[,dims, drop = FALSE]
  g <- scran::buildSNNGraph(sce, use.dimred = 'PCA_sub')
  ratio <- bluster::pairwiseModularity(g, seurat_obj@meta.data$seurat_clusters, as.ratio = TRUE)
  ratio_to_plot <- log10(ratio+1)
  clus_similarity_Plot <- ratio_to_plot %>%
    as_tibble() %>%
    rownames_to_column(var = 'cluster_1') %>%
    pivot_longer(
      cols = 2:ncol(.),
      names_to = 'cluster_2',
      values_to = 'probability'
    ) %>%
    mutate(
      cluster_1 = as.character(as.numeric(cluster_1) - 1),
      cluster_1 = factor(cluster_1, levels = rev(unique(cluster_1))),
      cluster_2 = factor(cluster_2, levels = unique(cluster_2))
    ) %>%
    ggplot(aes(cluster_2, cluster_1, fill = probability)) +
    geom_tile(color = 'white') +
    geom_text(aes(label = round(probability, digits = 2)), size = 2.5) +
    scale_x_discrete(name = 'Cluster', position = 'top') +
    scale_y_discrete(name = 'Cluster') +
    scale_fill_gradient(
      name = 'log10(ratio)', low = 'white', high = '#c0392b', na.value = '#bdc3c7',
      guide = guide_colorbar(
        frame.colour = 'black', ticks.colour = 'black', title.position = 'left',
        title.theme = element_text(hjust = 1, angle = 90),
        barwidth = 0.75, barheight = 10
      )
    ) +
    coord_fixed() +
    theme_bw() +
    theme(
      legend.position = 'right',
      panel.grid.major = element_blank()
   )
  wookieSay()
  return(clus_similarity_Plot)
}








wookieSay <- function() {
  messages <- c(
    "May the midichlorians be with you!",
    "This data analysis is strong with the Force.",
    "Your single-cell journey to the Endor system has begun.",
    "You must unlearn what you have learned about normalization.",
    "Do. Or do not. There is no try when it comes to dimensionality reduction.",
    "I find your lack of quality control disturbing.",
    "Great shot, kid! That's one in a million reads!",
    "Laugh it up, fuzzball! Your p-value adjustment is impressive.",
    "Hokey religions and ancient data analysis tools are no match for a good Seurat workflow.",
    "I've got a bad feeling about this cluster.",
    "The Force flows through every cell.",
    "Sacredible! Your single-cell data is a true Masterwork.",
    "Your analysis is an elegant weapon for a more civilized age.",
    "You don't know the power of the dark side of batch effects.",
    "Your mind powers will have doubled since the last time we met, count.",
    "This is the way to handle single-cell data.",
    "Search your feelings, you know it to be true... that normalization is necessary.",
    "I am one with the Force, and the Force is with me... and your single-cell data.",
    "You must learn the ways of the Force if you're to come with me to a higher dimension.",
    "Strap yourselves in, we're in for some fancy data integration!",
    "The ability to identify rare cell types is insignificant next to the power of the Force.",
    "You will never find a more wretched hive of scum and villainy than batch effects.",
    "These aren't the cells you're looking for. *waves hand*",
    "Your cells have paid the price for your lack of vision regarding normalization.",
    "I am altering the analysis. Pray I don't alter it any further.",
    "The circle is now complete. Your analysis has begun.",
    "You have controlled your Cell Ranger data. But you have allowed this Seurat... this Seurat to twist your mind.",
    "The data is strong with this one.",
    "I find your lack of faith in single-cell analysis disturbing.",
    "The Force will be with you. Always.",
    "I'll never join you in your pursuit of batch effects!",
    "You've failed, Your Highness. I am a Jedi, like my father before me. I will not turn to the dark side of poor quality control.",
    "The Force is what gives a Jedi their power. It's an energy field created by all living cells that surrounds us, penetrates us, and binds the galaxy together.",
    "I've been waiting for you, Obi-Wan. We meet again, at last. The circle is now complete. When I left you, I was but the learner; now I am the master of single-cell analysis.",
    "Your focus determines your reality. Focus on the data, not on your fears of misinterpreting it.",
    "You must confront your fear of batch effects. Only then will you become a true master of single-cell analysis.",
    "The fear of batch effects is the path to the dark side. Fear leads to anger, anger leads to hate, hate leads to suffering... in your data interpretation.",
    "May the Force guide your analysis, and your analysis guide the Force.",
    "In my experience, there's no such thing as too much quality control.",
    "Normalization is the path to the dark side. Normalization leads to hate. Hate leads to suffering.",
    "Remember, a Jedi's strength flows from the cell cycle.",
    "Impressive. Most impressive. But you are not a master of single-cell analysis yet.",
    "You have paid the price for your lack of vision regarding batch effects.",
    "The Force is strong with this cluster, but you are not a Jedi yet.",
    "You were the chosen one! You were meant to balance the data, not leave it in darkness!",
    "I have brought peace, freedom, justice, and security to my single-cell experiment.",
    "You underestimate the power of the dark side... of batch effects.",
    "I am a master of single-cell analysis, Darth Vader. I must not let myself be defeated.",
    "Your hate has made you powerful. But I sense there is still good in you, young analyst.",
    "I'll try normalization. That's a good trick!",
    "The Fear of Loss is a path to the dark side of data analysis.",
    "You're going to find that many of the truths we cling to depend greatly on our own point of view of the data.",
    "I've just about had enough of these Star Wars puns!",
    "I don't believe what I'm hearing. Normalization leading to hate? That's insane!",
    "You were supposed to bring balance to the data, not leave it in darkness!",
    "I fear nothing. For my data is all-powerful and you have lost my faith in single-cell analysis.",
    "Your journey towards the dark side of batch effects will be your undoing.",
    "Help me, Obi-Wan Kenobi. You're my only hope for proper normalization.",
    "The Force is strong in your family of cell types.",
    "I have a bad feeling about this... quality of your data.",
    "Use the Force, Luke... and Seurat's integration tools.",
    "Your data has a good motivation... for normalization.",
    "You don't know the power of the dark side... of overclustering!",
    "I find your lack of faith in single-cell analysis disturbing... and dangerous.",
    "You are unwise to lower your quality control thresholds.",
    "Ah, yes, the single-cell analyst... we have been expecting you.",
    "The data will decide your fate, not my bias.",
    "The ability to identify rare cell types is insignificant next to the power of the Force... and good experimental design.",
    "I sense a disturbance in the data... as if millions of cells cried out in terror and were suddenly silenced.",
    "You have controlled your Cell Ranger data, but you have allowed this Seurat... this Seurat to twist your mind.",
    "I am a master of single-cell analysis, Anakin. I must not let myself be defeated by batch effects.",
    "The Force is what gives a Jedi their power. It's an energy field created by all living cells that surrounds us, penetrates us, and binds the galaxy together... and your data.",
    "Your focus determines your reality. Focus on the data, not on your fears of misinterpreting it... or batch effects.",
    "The fear of batch effects is the path to the dark side. Fear leads to anger, anger leads to hate, hate leads to suffering... in your data interpretation and experimental design.",
    "Remember, a Jedi's strength flows from the cell cycle... and proper normalization.",
    "Impressive. Most impressive. But you are not a master of single-cell analysis yet... you still have much to learn.",
    "You have paid the price for your lack of vision regarding batch effects... and quality control.",
    "The Force is strong with this cluster, but you are not a Jedi yet... you must learn to trust your data.",
    "I have brought peace, freedom, justice, and security to my single-cell experiment... through the power of the Force and good experimental design.",
    "You underestimate the power of the dark side... of batch effects and overclustering.",
    "I am a master of single-cell analysis, Darth Vader. I must not let myself be defeated by your lack of faith in the data.",
    "Your hate has made you powerful. But I sense there is still good in you, young analyst... you must embrace the light side of proper normalization and quality control.",
    "I'll try normalization. That's a good trick... but you must also master batch effect correction and data integration.",
    "The Fear of Loss is a path to the dark side of data analysis... and poor experimental design.",
    "You're going to find that many of the truths we cling to depend greatly on our own point of view of the data... and our biases.",
    "I've just about had enough of these Star Wars puns... let's get back to the serious business of single-cell analysis!",
    "I don't believe what I'm hearing. Normalization leading to hate? That's insane! Proper normalization is the path to enlightenment and accurate data interpretation.",
    "Cells! Thousands of them! I sense a disturbance in the data.",
    "You don't know the power of the dark side of dropout events.",
    "I find your lack of faith in UMAP disturbing.",
    "The Force is strong with this cluster, but it's no Jedi yet.",
    "I've got a bad feeling about this... lack of quality control.",
    "These aren't the cell types you're looking for. *waves hand*",
    "Your hate for batch effects has made you powerful.",
    "Use the Force, Luke... and Seurat's data integration tools.",
    "The circle is now complete. Your single-cell analysis has begun.",
    "I am one with the Force, and the Force is with me... and your single-cell data.",
    "Normalization is the path to the light side of data analysis.",
    "You underestimate the power of the dark side... of overclustering!",
    "I sense a disturbance in the data... as if millions of cells cried out in terror and were suddenly silenced.",
    "Impressive. Most impressive. But you are not a master of single-cell analysis yet... you still have much to learn.",
    "The Fear of Loss is a path to the dark side of data analysis... and poor experimental design.",
    "I don't believe what I'm hearing. Normalization leading to hate? That's insane! Proper normalization is the path to enlightenment and accurate data interpretation.",
    "You were supposed to bring balance to the data, not leave it in darkness!",
    "The Force is what gives a Jedi their power. It's an energy field created by all living cells that surrounds us, penetrates us, and binds the galaxy together... and your data.",
    "Your focus determines your reality. Focus on the data, not on your fears of misinterpreting it... or batch effects.",
    "The fear of batch effects is the path to the dark side. Fear leads to anger, anger leads to hate, hate leads to suffering... in your data interpretation and experimental design.",
    "Remember, a Jedi's strength flows from the cell cycle... and proper normalization.",
    "You have paid the price for your lack of vision regarding batch effects... and quality control.",
    "I have brought peace, freedom, justice, and security to my single-cell experiment... through the power of the Force and good experimental design.",
    "I am a master of single-cell analysis, Darth Vader. I must not let myself be defeated by your lack of faith in the data.",
    "Your hate has made you powerful. But I sense there is still good in you, young analyst... you must embrace the light side of proper normalization and quality control.",
    "I'll try normalization. That's a good trick... but you must also master batch effect correction and data integration.",
    "You're going to find that many of the truths we cling to depend greatly on our own point of view of the data... and our biases.",
    "I've just about had enough of these Star Wars puns... let's get back to the serious business of single-cell analysis!",
    "The ability to identify rare cell types is insignificant next to the power of the Force... and good experimental design.",
    "You have controlled your Cell Ranger data, but you have allowed this Seurat... this Seurat to twist your mind.",
    "I am a master of single-cell analysis, Anakin. I must not let myself be defeated by batch effects.",
    "The fear of batch effects is the path to the dark side. Fear leads to anger, anger leads to hate, hate leads to suffering... in your data interpretation and experimental design.",
    "You underestimate the power of the dark side... of batch effects and overclustering.",
    "I sense a disturbance in the data... as if millions of cells cried out in terror and were suddenly silenced... by poor quality control.",
    "The Force is strong in your family of cell types... but you must learn to control your bias.",
    "You don't know the power of the dark side... of overclustering and improper normalization!",
    "Use the Force, Luke... and Seurat's powerful clustering algorithms.",
    "I find your lack of faith in single-cell analysis disturbing... and dangerous for your scientific career.",
    "You are unwise to lower your quality control thresholds... it leads to the dark side of data interpretation.",
    "Ah, yes, the single-cell analyst... we have been expecting you to join the light side of proper normalization.",
    "The data will decide your fate, not my bias... or my lack of proper experimental design.",
    "I sense a disturbance in the data... as if millions of cells cried out in terror and were suddenly silenced... by batch effects.",
    "You have controlled your Cell Ranger data, but you have allowed this Seurat... this Seurat to twist your mind... with improper clustering.",
    "I am a master of single-cell analysis, Anakin. I must not let myself be defeated by batch effects... or your lack of faith in quality control.",
    "The Force is what gives a Jedi their power. It's an energy field created by all living cells that surrounds us, penetrates us, and binds the galaxy together... and your data... and your experimental design.",
    "Your focus determines your reality. Focus on the data, not on your fears of misinterpreting it... or batch effects... or improper normalization.",
    "The fear of batch effects is the path to the dark side. Fear leads to anger, anger leads to hate, hate leads to suffering... in your data interpretation and experimental design... and your scientific career.",
    "Remember, a Jedi's strength flows from the cell cycle... and proper normalization... and a good experimental design.",
    "Impressive. Most impressive. But you are not a master of single-cell analysis yet... you still have much to learn... about batch effects, normalization, and quality control.",
    "You have paid the price for your lack of vision regarding batch effects... and quality control... and proper experimental design.",
    "The Force is strong with this cluster, but you are not a Jedi yet... you must learn to trust your data... and your experimental design.",
    "I have brought peace, freedom, justice, and security to my single-cell experiment... through the power of the Force and good experimental design... and proper normalization.",
    "You underestimate the power of the dark side... of batch effects and overclustering... and improper normalization.",
    "I am a master of single-cell analysis, Darth Vader. I must not let myself be defeated by your lack of faith in the data... or your poor experimental design.",
    "Your hate has made you powerful. But I sense there is still good in you, young analyst... you must embrace the light side of proper normalization and quality control... and a good experimental design."
  )
    
 
  
  message <- sample(messages, 1)
  cat("\n")
  cat('WookieSays:', message, "\n")
}
custom_colors <- list()
colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)
colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)
custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)