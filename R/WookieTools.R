# WookieTools

# Version 0.3.2
# All input objects are Seurat Objects unless mentioned otherwise
#' @export
load_libraries<- function(){
  suppressMessages(library(Matrix))
  suppressMessages(library(dplyr))
  suppressMessages(library(tidyr))
  suppressMessages(library(tidyverse))
  suppressMessages(library(scds))
  suppressMessages(library(SingleCellExperiment))
  suppressMessages(library(Seurat))
  suppressMessages(library(ggplot2))
  suppressMessages(library(plyr))
  suppressMessages(library(cowplot))
  suppressMessages(library(patchwork))
}

load_libraries()

# Seurat Object Quality Control function
# Run Iteratively
# change 'mt' to 'MT' depending on Mouse/Human dataset
# Potential error with feature name, change as needed
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
#' @param colors Colors for facetting ...
#' @return Seurat object after quality control
#' @export
wookie_qc <- function(matrix, nf_min = 0, nf_max = 20000, nc = 200000, pmt = 20, ptr = 100,group = 'orig.ident', colors = NULL) {
  load_libraries()
  matrix[['percent.ribo']] <- PercentageFeatureSet(matrix, pattern = "^Rp[sl]")
  matrix[["percent.mt"]] <- PercentageFeatureSet(matrix, pattern = "^mt-")
  matrix <- subset(matrix, subset = nf_min < nFeature_RNA & nFeature_RNA < nf_max & nCount_RNA < nc & percent.mt < pmt & percent.ribo < ptr)
  options(repr.plot.width = 16, repr.plot.height = 30) 
  vl_plot <- VlnPlot(matrix, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
                     ncol = 4)
  options(repr.plot.width = 16, repr.plot.height = 35)
  g <- FeatureScatter(matrix,
                      feature1 = "nCount_RNA",
                      feature2 = "nFeature_RNA",
                      group.by = group, shape = 1) + ggtitle('nCount_vs_nFeature') +
    geom_point(shape = 1, alpha = 0.3)
  
  plot2 <- FeatureScatter(matrix,
                          feature1 = "percent.mt",
                          feature2 = "nFeature_RNA") + ggtitle('percent.mt_vs_nFeature') +
    geom_point(shape = 1, alpha = 0.3) + facet_wrap(~colors)
  
  plot3 <- FeatureScatter(matrix,
                          feature1 = "percent.ribo",
                          feature2 = "nFeature_RNA") + ggtitle('percent.ribo_vs_nFeature') +
    geom_point(shape = 1, alpha = 0.3) + facet_wrap(~colors)
  
  plot4 <- FeatureScatter(matrix,
                          feature1 = "percent.mt",
                          feature2 = "percent.ribo") + ggtitle('percent.mt_vs_percent.ribo') +
    geom_point(shape = 1, alpha = 0.3) + facet_wrap(~colors)
  
  # Combine all plots into a single plot
  c_plot<- plot_grid(g, plot2, plot3, plot4, ncol = 2,label_size = 10)
  combined_plot<-  plot_grid(vl_plot,c_plot,ncol = 1,nrow = 2,rel_widths = c(2,1))
  # Display the combined plot
  print(combined_plot)
  
  return(matrix)
}


# Doublet finder using scds
# Not recommended, use scrubdub (scrublet) instead
#' @name scds_doublets
#' @title scds Doublet finder for Seurat
#' @description Use not recommended ...
#' @param matrix Seurat object ...
#' @return Seurat object with doublet scores and call
#' @export
scds_doublets <- function(matrix){
  load_libraries()
  print('Not recommended!, use scrubdub (Scrublet) instead')
  suppressMessages(a_matrix <- NormalizeData(matrix))
  suppressMessages(a_matrix <- FindVariableFeatures(a_matrix, selection.method = "vst", nfeatures = 3000))
  suppressMessages(a_matrix <- ScaleData(a_matrix))
  suppressMessages(a_matrix <- RunPCA(a_matrix))
  suppressMessages(a_matrix <- RunUMAP(a_matrix, dims = 1:10))  
  sce <- as.SingleCellExperiment(a_matrix)
  sce = bcds(sce, retRes = TRUE, estNdbl=TRUE)
  ## Annotate doublet using co-expression based doublet scoring:
  try({
    sce = cxds(sce, retRes = TRUE, estNdbl=TRUE)
  })
  ### If cxds worked, run hybrid, otherwise use bcds annotations
  if ("cxds_score" %in% colnames(colData(sce))) {
    ## Combine both annotations into a hybrid annotation
    sce = cxds_bcds_hybrid(sce, estNdbl=TRUE)
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$hybrid_score, colData(sce)$hybrid_call))
  } else {
    print("this pool failed cxds so results are just the bcds calls")
    Doublets <- as.data.frame(cbind(rownames(colData(sce)), colData(sce)$bcds_score, colData(sce)$bcds_call))
  }
  ## Doublet scores are now available via colData:
  colnames(Doublets) <- c("Barcode","scds_score","scds_DropletType")
  Doublets$scds_DropletType <- gsub("FALSE","singlet",Doublets$scds_DropletType)
  Doublets$scds_DropletType <- gsub("TRUE","doublet",Doublets$scds_DropletType)
  a_matrix@meta.data$Doublet_score = as.numeric(Doublets$scds_score)
  matrix@meta.data$Doublet_score = as.numeric(Doublets$scds_score)
  a_matrix@meta.data$Doublet_score
  matrix@meta.data$Doublet_type = Doublets$scds_DropletType
  DimPlot(a_matrix, reduction = 'umap')
  hist(matrix@meta.data$Doublet_score)
  return(matrix)
}

# Doublet finder using Scrublet
#' @name scrub_dub
#' @title Scrublet for Seurat ...
#' @description Run Scrublet on a Seurat Object ...
#' @param seu_obj Seurat object ...
#' @return Seurat object with scrublet scores and call
#' @export
scrub_dub <- function(seu_obj){
  load_libraries()
  suppressMessages(seu_obj <- NormalizeData(seu_obj))
  suppressMessages(seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 3000))
  suppressMessages(seu_obj <- ScaleData(seu_obj))
  suppressMessages(seu_obj <- RunPCA(seu_obj))
  sce_obj <- as.SingleCellExperiment(seu_obj)
  sce_scrub <- runScrublet(sce_obj)
  seu_obj_scrubbed <- as.Seurat(sce_scrub)
  scrub_scores <- seu_obj_scrubbed@meta.data$scrublet_score
  scrub_type <- seu_obj_scrubbed@meta.data$scrublet_call
  seu_obj@meta.data$scrublet_score = scrub_scores
  seu_obj@meta.data$scrublet_call = scrub_type
  hist(scrub_scores)
  return(seu_obj)
}

# Function to sum the counts of two matrices containing the same cells
# Input: Count Matrices | Output: Seurat Object
#' @name sum_matrices
#' @title Matrix Sum function
#' @description Merge two count matrices, where the cells are the same, to obtain a single seurat object with the counts from two matrices summed for each cell ...
#' @param matrix1 count matrix ...
#' @param matrix2 count matrix ...
#' @param sample sample/library name ...
#' @param min_cells minimum number cells a gene is found in ...
#' @param min_features minimum number of features found in a cell ...
#' @return Summed and merged Seurat object
#' @export
sum_matrices <- function(matrix1, matrix2, sample = 'sample', min_cells = 3, min_features = 200) {
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
  matrix1_common <- matrix1[which(rownames(matrix1) %in% common_rows), which(colnames(matrix1) %in% common_cols)]
  matrix2_common <- matrix2[which(rownames(matrix2) %in% common_rows), which(colnames(matrix2) %in% common_cols)]
  
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
  seu_obj <- CreateSeuratObject(result_matrix, min.cells = min_cells, min.features = min_features, project = sample)
  
  return(seu_obj)
}

# Function to run and plot multiple UMAP's for different numbers of a feature(Highly variable genes or Most Abundant genes)
# Features must be obtained and given as input
#' @name plot_multi_feat_umap
#' @title Plot UMAPs to test features
#' @description plot multiple UMAP's for different numbers of a feature i.e Highly variable genes or Most Abundant genes ...
#' @return plot saved to global environment
#' @export
plot_multi_feat_umap <- function(object = seu_obj, features = features, min.dist = 0.1, 
                                 max_features = 3000,ftype='HVG',
                                 step = 500,out_name = 'combined_umap') {
  load_libraries()
  plot_list <- list()
  
  for (feature_length in seq(500, max_features, step)) {
    current_features <- features[1:feature_length]
    cat(paste0('Calculating UMAP at ',ftype,':',feature_length))
    current_umap <- RunUMAP(object, features = current_features, min.dist = min.dist)
    current_plot <- DimPlot(current_umap, reduction = 'umap') + ggtitle(paste('UMAP ',ftype, feature_length))
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

# Function plots multiple UMAP's at different min.dist values
#' @name plot_multi_min_dist_umap
#' @title Plot UMAPs to test min.dist
#' @description plots multiple UMAP's at different min.dist values ...
#' @return plot saved to global environment
#' @export
plot_multi_min_dist_umap <- function(object = seu_obj, features = features, 
                                     out_name = 'min_dist_umaps'){
  load_libraries()
  plot_list <- list()
  
  for (min_dist in seq(0.1, 0.5, 0.1)) {
    current_md <- min_dist
    cat(paste0('Calculating UMAP at min.dist:',current_md))
    current_umap <- RunUMAP(object, features = features, min.dist = current_md)
    current_plot <- DimPlot(current_umap, reduction = 'umap') + ggtitle(paste('UMAP: min.dist:', current_md))
    plot_list[[length(plot_list) + 1]] <- current_plot
    cat(paste0('UMAP done at min.dist:',current_md))
  }
  
  # Combine plots into a grid
  combined_plot <- plot_grid(plotlist = plot_list)
  
  # Assign the combined plot to a variable in the global environment
  assign(out_name, combined_plot, envir = globalenv())
  
  # Return the combined plot
  return(combined_plot)
}

# Function to plot multiple features with color map
# Input seurat object and a list of features to plot
#' @name multi_f_plots
#' @title Plot multiple features 
#' @description plots n number features given in a list with custom colour scale ...
#' @param featureList list of features to plot ...
#' @return plot saved to global environment
#' @export
multi_f_plots <- function(seuratObject, featureList, ncol = 3, pt.size = 0.8) {
  load_libraries()
  plotList <- lapply(featureList, function(feature) {
    FeaturePlot(object = seuratObject, features = feature, pt.size = pt.size, reduction = "umap") +
      theme(aspect.ratio = 1) +
      scale_color_gradientn(colours = c("#DCDCDC", "yellow", "orange", "red", "#8b0000"))
  })
  
  plotGrid <- plot_grid(plotlist = plotList, ncol = ncol, rel_widths = rep(1, length(featureList)))
  
  return(plotGrid)
}



# Function to remove particular cell types based on marker gene expression
#' @name wookie_filter_celltypes
#' @title wookie_filter_celltypes
#' @description Function to remove particular cell types based on marker gene expression ...
#' @param seurat_object object to filter
#' @param marker_List list of marker genes of celltype to remove ...
#' @param cutoff quantile threshold, default is 0.99
#' @return filtered seurat object
#' @export
wookie_filter_celltype <- function(seurat_object,marker_list,cutoff = 0.99){
  
  print('Ensure marker genes are in RNA$scale.data')
  
  # Transpose the expression matrix
  expression_matrix_transposed <- t(seurat_object@assays$RNA$scale.data)
  
  # Calculate average normalized values for fibroblast marker genes
  seurat_object$avg_celltype_expression <- rowMeans(expression_matrix_transposed[, marker_list, drop = FALSE])
  
  # Set a threshold for filtering
  threshold <- quantile(seurat_object$avg_celltype_expression, cutoff)
  
  # Filter out non-epithelial cells
  cellstokeep <- which(seurat_object$avg_celltype_expression <= threshold)
  seu_filtered <- seurat_object[, cellstokeep]
  cellstoremove <- which(seurat_object$avg_celltype_expression > threshold)
  seu_removed <- seurat_object[, cellstoremove]
  print(paste0(length(cellstokeep), ' Cells kept, ', length(cellstoremove), ' cells removed.'))
  return(seu_filtered)
}



# Function to plot qc metrics of a sparce matrix
#' @name wookie_matrix_qc_plot
#' @title wookie_matrix_qc_plot
#' @description Function to to plot qc metrics of a sparce matrix ...
#' @param count_matrix_sparse sparce matrix
#' @param fill_color plot color, default is #589FFF ...
#' @param title title for the plot
#' @return qc plot
#' @export
wookie_matrix_qc_plot <- function(count_matrix_sparse, fill_color = "#589FFF",title="") {
  reads_per_cell <- Matrix::colSums(count_matrix_sparse)
  genes_per_cell <- Matrix::colSums(count_matrix_sparse > 0)
  reads_per_gene <- Matrix::rowSums(count_matrix_sparse > 0)
  
  p1 <- ggplot() +
    geom_histogram(aes(x = log10(reads_per_cell + 1)), fill = fill_color, color = 'black', bins = 30) +
    ggtitle('reads per cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p2 <- ggplot() +
    geom_histogram(aes(x = log10(genes_per_cell + 1)), fill = fill_color, color = 'black', bins = 30) +
    ggtitle('genes per cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p4 <- ggplot() +
    geom_histogram(aes(x = log10(reads_per_gene + 1)), fill = fill_color, color = 'black', bins = 30) +
    ggtitle('reads per gene') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  p3 <- ggplot() +
    geom_point(aes(x = reads_per_cell, y = genes_per_cell), fill = fill_color,color='black',pch=21, shape = 16, size = 2, alpha = 1) +
    ggtitle('Reads vs. Genes per Cell') +
    theme_minimal() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  qplot <- plot_grid(p1, p2, p3, p4, ncol = 2) + ggtitle(title)
  
  
  return(qplot)
}


#' Function to compare different normalization methods
#' @title plot_compare_normalisation
#' @name plot_compare_normalisation
#' @param seuratObj A Seurat object with both RNA and SCT assays
#' @return A plot comparing raw counts, normalized counts, SCTransform counts, and scaled data
#' @export
plot_compare_normalisation <- function(seuratObj) {
  
  # Extract count data
  rc <- GetAssayData(object = seuratObj,assay = 'RNA',layer = 'counts')
  rawCounts <- colSums(rc)
  normalizedCounts <- colSums(seuratObj[['RNA']]@layers$data)
  sctransformCounts <- colSums(seuratObj[['SCT']]@data)
  scaledLN <- colSums(seuratObj[['RNA']]@layers$scale.data)
  scaledSCT <- colSums(seuratObj[['SCT']]@scale.data)
  
  # Create a data frame for plotting
  plotDataCounts <- data.frame(
    Cell = names(rawCounts),
    Raw_Counts = rawCounts,
    Normalized_Counts = normalizedCounts,
    SCT_Counts = sctransformCounts,
    Scaled_LN = scaledLN,
    Scaled_SCT = scaledSCT
  )
  
  # Melt data for easier plotting
  plotDataMelted <- reshape2::melt(plotDataCounts, id.vars = "Cell", variable.name = "Type", value.name = "Counts")
  
  # Plot
  plot <- ggplot(plotDataMelted, aes(x = Cell, y = Counts, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ Type, scales = "free_y") +
    labs(title = "Comparison of Normalization Methods",
         x = "Cell",
         y = "Counts",
         fill = "Normalization Method") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) # Remove x-axis labels
  
  return(plot)
}

# Function to compare Log Normalisation and SCT (Histogram)
#' @name plot_gene_expression_histogram
#' @title plot_gene_expression_histogram
#' @description Function to compare Log Normalisation and SCT, Ensure Seurat Obj is SCTransformed and both scale.data contains all genes ...
#' @param seurat_obj Seurat Object with both RNA and SCT assays
#' @return histogram
#' @export
plot_gene_expression_histogram <- function(seurat_obj) {
  expression_data_RNA <- as.vector(seurat_obj@assays$SCT$scale.data)
  expression_data_SCT <- as.vector(seurat_obj@assays$RNA$scale.data)
  
  mean_expr_rna <- mean(expression_data_RNA)
  sd_expr_rna <- sd(expression_data_RNA)
  mean_expr_sct <- mean(expression_data_SCT)
  sd_expr_sct <- sd(expression_data_SCT)
  
  par(mfrow = c(1, 2))
  
  hist(expression_data_RNA, breaks = 50, freq = FALSE, main = "Histogram of Gene Expression RNA", xlab = "Gene Expression", col = "lightgray")
  curve(dnorm(x, mean = mean_expr_rna, sd = sd_expr_rna), add = TRUE, col = "blue", lwd = 2)
  
  hist(expression_data_SCT, breaks = 50, freq = FALSE, main = "Histogram of Gene Expression SCT", xlab = "Gene Expression", col = "lightgray")
  curve(dnorm(x, mean = mean_expr_sct, sd = sd_expr_sct), add = TRUE, col = "blue", lwd = 2)
  
  par(mfrow = c(1, 1))
}
