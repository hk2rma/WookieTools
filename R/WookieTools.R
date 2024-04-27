# WookieTools - Version 0.7

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

#' @name wookie_featureplot
#' @title Plot multiple features with color map
#' @import Seurat
#' @import ggplot2
#' @import cowplot
#' @description Plot multiple features with custom color scale
#' @param seuratObject Seurat object
#' @param featureList List of features to plot
#' @param ncol Number of columns for arrangement
#' @param pt.size Point size for plotting
#' @param split_by Variable for splitting
#' @param alpha Alpha transparency value
#' @param order Whether to order cells based on expression
#' @param min.cutoff Minimum cutoff value for expression
#' @param max.cutoff Maximum cutoff value for expression
#' @param reduction Dimensional reduction technique to use
#' @param keep.scale Whether to keep the same scale for all features or scale each individually
#' @param shape.by Variable to represent groups by shape
#' @param slot slot to use for plotting, default is 'data'
#' @param blend Whether to blend colors for continuous data
#' @param blend.threshold Threshold for blending colors
#' @param label Whether to label cells
#' @param label.size Size of cell labels
#' @param label.color Color of cell labels
#' @param repel Whether to repel labels to avoid overlapping
#' @param coord.fixed Whether to fix the coordinate aspect ratio
#' @param by.col Whether to group cells by color
#' @param sort.cell Deprecated argument
#' @param interactive Whether to make the plot interactive
#' @param combine Whether to combine plots into a single grid
#' @param raster Optional raster object for plotting
#' @param raster.dpi Resolution for raster plotting
#' @return Plot grid
#' @export
wookie_featureplot <- function(seuratObject, featureList, ncol = 3,
                  pt.size = 0.8,split_by = NULL,alpha = 1,
                  order = FALSE,
                  min.cutoff = NA,
                  max.cutoff = NA,
                  reduction = NULL,
                  keep.scale = "feature",
                  shape.by = NULL,
                  slot = "data",
                  blend = FALSE,
                  blend.threshold = 0.5,
                  label = FALSE,
                  label.size = 4,
                  label.color = "black",
                  repel = FALSE,
                  coord.fixed = FALSE,
                  by.col = TRUE,
                  sort.cell = deprecated(),
                  interactive = FALSE,
                  combine = TRUE,
                  raster = NULL,
                  raster.dpi = c(512, 512)) {
  plotList <- lapply(featureList, function(feature) {
    FeaturePlot(object = seuratObject, features = feature, 
                pt.size = pt.size,
                alpha = alpha,
                order = order,
                min.cutoff = min.cutoff,
                max.cutoff = max.cutoff,
                reduction = reduction,
                split.by = split_by,
                keep.scale = keep.scale,
                shape.by = shape.by ,
                slot = slot,
                blend = blend,
                blend.threshold = blend.threshold,
                label = label,
                label.size = label.size,
                label.color = label.color,
                repel = repel,
                ncol = ncol,
                coord.fixed = coord.fixed,
                by.col = by.col,
                interactive = interactive,
                combine = combine,
                raster = raster,
                raster.dpi = raster.dpi,) +
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

# Function to plot Jaccard similarity index
#' @name wookie_jaccardPlot
#' @title Function to plot Jaccard similarity index
#' @import Seurat
#' @import reshape2
#' @import viridis
#' @import ggplot2
#' @import cluster
#' @description Function to plot Jaccard similarity index
#' @param seurat_obj Seurat object
#' @param clusters Name of the column in the Seurat object containing cluster assignments, default is 'seurat_clusters'
#' @param logfc.threshold Log fold-change threshold, default is 0.95
#' @param min.pct Minimum percentage of cells expressing a gene, default is 0.25
#' @param test.use Statistical test to use for finding differentially expressed genes, default is "wilcox"
#' @param fdr.threshold False discovery rate threshold, default is 0.05
#' @param p_val_bonf.threshold Bonferroni-corrected p-value threshold, default is 0.05
#' @param avg.log2fc.threshold Average log2 fold-change threshold, default is 0.25
#' @param title Title for the plot, default is 'Wookie_Jaccard'
#' @param saveplot Boolean indicating whether to save the plot, default is FALSE
#' @param height Height of the plot in the specified units, default is 10
#' @param width Width of the plot in the specified units, default is 10
#' @param dpi Resolution of the saved plot, default is 700
#' @param units Units for the plot dimensions, default is 'cm'
#' @param limitsize Boolean indicating whether to limit the size of the plot, default is FALSE
#' @return plot Jaccard similarity index
#' @export
wookie_jaccardPlot <- function(seurat_obj,clusters = 'seurat_clusters',
                               logfc.threshold = 0.95,
                               min.pct = 0.25,
                               test.use = "wilcox",
                               fdr.threshold = 0.05,
                               p_val_bonf.threshold = 0.05,
                               avg.log2fc.threshold = 0.25,
                               title = 'Wookie_Jaccard',
                               saveplot = FALSE,
                               height = 10, width = 10,
                               dpi = 700,units = 'cm',limitsize = FALSE){
  
  seurat_obj$seurat_clusters <- seurat_obj[[clusters]]
  Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
  clusters <- levels(seurat_obj)
  # Initialize empty lists to store markers and DEGs for each cluster
  markers_list <- list()
  degs_list <- list()
  
  # Iterate over each cluster
  for (i in clusters) {
    # Define clusters to compare with
    compare_me <- clusters[clusters != i]
    
    # Find markers for the current cluster
    markers_i <- FindMarkers(seurat_obj, ident.1 = i, ident.2 = compare_me,
                             logfc.threshold = logfc.threshold, min.pct = min.pct,
                             test.use = test.use, verbose = FALSE)
    markers_i$fdr <- p.adjust(markers_i$p_val, method = "fdr") 
    markers_i <- markers_i[markers_i$fdr < fdr.threshold & 
                             markers_i$p_val_adj < p_val_bonf.threshold &
                             markers_i$avg_log2FC > avg.log2fc.threshold , ]
    
    # Store markers in the markers_list
    markers_list[[i]] <- markers_i
    
    # Store DEGs in the degs_list
    degs_list[[i]] <- rownames(markers_i)
  }
  
  # Calculate pairwise Jaccard similarities
  jaccard_similarities <- sapply(degs_list, function(x) {
    sapply(degs_list, function(y) {
      length(intersect(x, y)) / length(union(x, y))
    })
  })
  
  # Convert the similarity matrix to a data frame
  jaccard_df <- as.data.frame(jaccard_similarities)
  jaccard_df$Cluster1 <- rownames(jaccard_df)
  jaccard_df <- melt(jaccard_df, id.vars = "Cluster1", variable.name = "Cluster2",
                     value.name = "Jaccard_Similarity")
  
  # Convert Cluster1 and Cluster2 to factors with specified levels
  jaccard_df$Cluster1 <- factor(jaccard_df$Cluster1, levels = clusters)
  jaccard_df$Cluster2 <- factor(jaccard_df$Cluster2, levels = clusters)
  
  # Plot the Jaccard similarities
  js_plot <- ggplot(jaccard_df, aes(x = Cluster1, y = Cluster2, fill = Jaccard_Similarity)) +
    geom_tile() +
    geom_text(aes(label = round(Jaccard_Similarity, 2)), color = "white", size = 3) +
    scale_fill_gradientn(colours = viridis(20,direction = 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
    labs(x = "Cluster", y = "Cluster", fill = "Jaccard Similarity",
         title = "Jaccard Similarity of DEGs between Clusters",
         caption = title) +
    coord_fixed(ratio = 1)
  
  if(saveplot){
    ggsave('jaccard.jpeg', js_plot,
           height = height, width = width,
           dpi = dpi,units = units,limitsize = limitsize)
  }
  wookieSay()
  return(js_plot)
  
}

# Function to get qc stats
#' @name wookie_stats
#' @title Function to get qc stats
#' @import Seurat
#' @description Function to get qc stats
#' @param seurat_obj Seurat object
#' @param counts default is 'nCount_RNA'
#' @param features default is 'nFeature_RNA'
#' @param mito default is 'percent.mt'
#' @param ribo default is NULL
#' @return stats df
#' @export
wookie_stats <- function(seurat_obj, counts = 'nCount_RNA', features = 'nFeature_RNA', mito = 'percent.mt', ribo = NULL) {
  nc <- obj[[counts]]
  nf <- obj[[features]]
  pmt <- obj[[mito]]
  
  # Calculate mean, standard deviation, and variance for each value
  nc_mean <- mean(nc)
  nc_sd <- sd(nc)
  nc_var <- var(nc)
  
  nf_mean <- mean(nf)
  nf_sd <- sd(nf)
  nf_var <- var(nf)
  
  pmt_mean <- mean(pmt)
  pmt_sd <- sd(pmt)
  pmt_var <- var(pmt)
  
  # Handle the case when ribo is NULL
  if (!is.null(ribo)) {
    ptr <- obj[[ribo]]
    ptr_mean <- mean(ptr)
    ptr_sd <- sd(ptr)
    ptr_var <- var(ptr)
    results_df <- data.frame(
      Metric = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
      Mean = c(nc_mean, nf_mean, pmt_mean, ptr_mean),
      SD = c(nc_sd, nf_sd, pmt_sd, ptr_sd),
      Variance = c(nc_var, nf_var, pmt_var, ptr_var)
    )
  } else {
    results_df <- data.frame(
      Metric = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
      Mean = c(nc_mean, nf_mean, pmt_mean),
      SD = c(nc_sd, nf_sd, pmt_sd),
      Variance = c(nc_var, nf_var, pmt_var)
    )
  }
  
  return(results_df)
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
    "You've failed, Your Highness. I am a Jedi, like my father before me. I will not turn to the dark side of poor quality control.",
    "The Force is what gives a Jedi their power. It's an energy field created by all living cells that surrounds us, penetrates us, and binds the galaxy together.",
    "I've been waiting for you, Obi-Wan. We meet again, at last. The circle is now complete. When I left you, I was but the learner; now I am the master of single-cell analysis.",
    "Your focus determines your reality. Focus on the data, not on your fears of misinterpreting it.",
    "You must confront your fear of batch effects. Only then will you become a true master of single-cell analysis.",
    "The fear of batch effects is the path to the dark side. Fear leads to anger, anger leads to hate, hate leads to suffering... in your data interpretation.",
    "May the Force guide your analysis, and your analysis guide the Force.",
    "In my experience, there's no such thing as too much quality control.",
    "Remember, a Jedi's strength flows from the cell cycle.",
    "Impressive. Most impressive. But you are not a master of single-cell analysis yet.",
    "You were the chosen one! You were meant to balance the data, not leave it in darkness!",
    "I have brought peace, freedom, justice, and security to my single-cell experiment.",
    "You underestimate the power of the dark side... of batch effects.",
    "I'll try normalization. That's a good trick!",
    "You're going to find that many of the truths we cling to depend greatly on our own point of view of the data.",
    "I've just about had enough of these Star Wars puns!",
    "I fear nothing. For my data is all-powerful and you have lost my faith in single-cell analysis.",
    "Your journey towards the dark side of batch effects will be your undoing.",
    "The Force is strong in your family of cell types.",
    "I have a bad feeling about this... quality of your data.",
    "Use the Force, Luke... and Seurat's integration tools.",
    "Your data has a good motivation... for normalization.",
    "You don't know the power of the dark side... of overclustering!",
    "You are unwise to lower your quality control thresholds.",
    "Ah, yes, the single-cell analyst... we have been expecting you.",
    "The data will decide your fate, not my bias.",
    "I sense a disturbance in the data... as if millions of cells cried out in terror and were suddenly silenced.",
    "The Force is what gives a Jedi their power. It's an energy field created by all living cells that surrounds us, penetrates us, and binds the galaxy together... and your data.",
    "Your focus determines your reality. Focus on the data, not on your fears of misinterpreting it... or batch effects.",
    "Cells! Thousands of them! I sense a disturbance in the data.",
    "You don't know the power of the dark side of dropout events.",
    "These aren't the cell types you're looking for. *waves hand*",
    "Your hate for batch effects has made you powerful.",
    "Normalization is the path to the light side of data analysis.",
    "I've just about had enough of these Star Wars puns... let's get back to the serious business of single-cell analysis!",
    "This is not another star wars pun. *waves hand*",
    "May the Force of high-throughput sequencing be with you.",
    "This is your destiny. Join me, and together, we can rule the galaxy of single-cell analysis.",
    "In my experience, there’s no such thing as luck in single-cell transcriptomics.",
    "Remember... the Force will be with you, always, and so will your flow cytometer.",
    "Help me, Single-Cell Analyzer. You’re my only hope.",
    "It's a trap! Beware of overfitting in your analysis.",
    "The ability to destroy a planet is insignificant next to the power of the single-cell.",
    "Use the Force, think creatively for novel cell type discovery.",
    "A Jedi uses the Force for knowledge and defense, never for missing data.",
    "Fear is the path to the dark side. Fear leads to bias, bias leads to skew, skew leads to misinterpretation in data analysis.",
    "Your overclustering will be your undoing.",
    "You were the chosen one! It was said that you would destroy the batch effects, not join them!",
    "I sense much sequencing in you.",
    "When nine hundred years old you reach, look as good your longitudinal data will not.",
    "Much to learn, you still have. This is just the beginning of your single-cell journey.",
    "The Force will be with you. Always. Especially in high-dimensional data analysis.",
    "Size matters not. Look at me. Judge me by my cell size, do you?",
    "Never tell me the odds of this stochastic analysis!",
    "Let go of your redundant variables, and only then, a true model of the data you will see.",
    "I have a good feeling about this hypothesis.",
    "Difficult to see. Always in motion is the future of single-cell technologies.",
    "Search your data. You know it to be clean.",
    "Only a master of alignment can understand the true structure of the genome.",
    "Let's keep the clonotypes balanced, as all things should be.",
    "A Jedi seeks not these false discoveries.",
    "Patience you must have, my young Padawan, to interpret these single-cell experiments.",
    "A disturbance in the data, there is. Prepare to visualize!",
    "Impressive. Most impressive. Your model has now conquered the dimensionality reduction challenge.",
    "Your faith in your pipeline is your strength... and your weakness.",
    "Heed the prophecy of cross-validation. Do not underestimate its power.",
    "May your cluster resolution be fine and your cell identities clear.",
    "Escape the dark side of inadequate sample sizes.",
    "Yoda warned of the dark side: Overfitting, Batch Effects, and the dreaded Technical Variability.",
    "Come to the single-cell side; we have the definitive resolution.",
    "Rebellions are built on hope—and robust data analysis.",
    "To defeat batch effects, into the normalization you must go.",
    "The time to publish, it is. Ready your manuscripts we must!",
    "Chewbacca might not understand RNA-seq, but he’d appreciate a good co-expression network.",
    "Like the Millennium Falcon, your data’s journey through the analysis pipeline must be swift and sure.",
    "MIDI-chlorians: the bioinformatics tools you are looking for.",
    "Trust in the Force, but also in your control samples.",
    "These are the significant results you are looking for. *waves hand*",
    "Deploy the scatter plots and PCA, and prepare for the defense of your findings!",
    "Not as clumsy or random as bulk sequencing; an elegant tool for a more civilized age—single-cell RNA-seq.",
    "Rise, my friend. The data does not betray you—it makes you powerful!",
    "If you strike this outlier down, I shall become more powerful than you could possibly imagine."
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