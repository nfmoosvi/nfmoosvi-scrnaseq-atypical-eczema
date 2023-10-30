### Mouse dataset: XXX
# Box folders: https://ucsf.box.com/s/jeqp1aix9hrp41xgz1af0agnojzanuvl (atypical pre-treatment)
# and https://ucsf.box.com/s/j16rueszrmcghmvplbrjd3zzsk92pipx (atypical mid-treatment)
# Link to FindConservedMarkers output of this analysis: https://docs.google.com/spreadsheets/d/1YhsLmJe_Ygwn8B8iU7h5SS32Jllk9OrASTpg2f2J3QU/edit?usp=sharing
# Link to FindMarkers output of this analysis: https://docs.google.com/spreadsheets/d/1o7PZRBuoB0ySXrtaflQ4rM-jpfE0Kew92ocZyK9lZ9Q/edit?usp=sharing

### Follow-Ups/Questions/TODOs:
#
#  1) Code is currently not modular. Need to clean up, refactor, and add comments. Also look
# into switching over to use lab's newer pipeline.
# 
#  2) Need to modify preprocessing/filtering logic.
#
#  3) Need to determine labels for unlabeled clusters.

# Load libraries

library(cowplot)
library(devtools)
library(dplyr)
library(future)
library(ggExtra)
library(ggplotify)
library(gplots)
library(harmony)
library(here)
library(Matrix)
library(metap)
library(mvtnorm)
library(openxlsx)
library(patchwork)
library(rcartocolor)
library(readr)
library(Seurat)
library(stringi)
library(tidyverse)
library(uwot)
library(wutilities)

# Enable parallel processing

plan("multisession")
options(future.globals.maxSize = 8000 * 1024^2)

# Define helper functions

trim.feature.names <- function(inmat){
  newnames <- sapply(strsplit(rownames(inmat), "-"),
                     function(x) {
                       if(length(x) == 1) return(x)
                       else return(paste(x[-length(x)], collapse = "-"))
                     })
  rownames(inmat) <- newnames
  return(inmat)
}

# Define paths

input_data_root_dir <- "./analysis/atypical_eczema/input/nfmoosvi_data_20231023"
output_data_root_dir <- "./analysis/atypical_eczema/output/nfmoosvi_data_20231023"

# Paths for atypical eczema (AE)

ae_input_data_root_dir <- paste(input_data_root_dir, "/ae_samples/", sep = "")

ae_folders <- c()
ae_folders[[1]] <- "skin_357_pre_treatment"
ae_folders[[2]] <- "skin_385_mid_treatment"

cur_output_data_dir <- paste(output_data_root_dir, "/pv_ad_hc_ae/", sep = "")
cur_output_data_dir_eda <- paste(cur_output_data_dir, "eda/", sep = "")

# Load data

for (i in 1:2){
  # AE
  cur_dir <- paste(ae_input_data_root_dir, ae_folders[i], "/sample_filtered_feature_bc_matrix", sep = "")
  cur_obj_name <- paste0("ae", i, ".data")
  assign(cur_obj_name, Read10X(data.dir = cur_dir, unique.features = TRUE))
}

# Create seurat objects from input data

ae1 <- CreateSeuratObject(counts = ae1.data, project = "ae1", min.cells = 3, min.features = 200)
ae2 <- CreateSeuratObject(counts = ae2.data, project = "ae2", min.cells = 3, min.features = 200)
ae1[["percent.mt"]] <- PercentageFeatureSet(ae1, pattern = "^mt-")
ae2[["percent.mt"]] <- PercentageFeatureSet(ae2, pattern = "^mt-")

# Visualize QC metrics and save plots

# Add metrics for % of genes that are mitochondrial and ribosomal

ae1 <- PercentageFeatureSet(ae1, pattern = "^MT-|mt-", col.name = "pct_mito")
ae1 <- PercentageFeatureSet(ae1, pattern = "^RP[SL]|^MRP[SL]|Rp[sl]|Mrp[sl]", col.name = "pct_ribo")

ae2 <- PercentageFeatureSet(ae2, pattern = "^MT-|mt-", col.name = "pct_mito")
ae2 <- PercentageFeatureSet(ae2, pattern = "^RP[SL]|^MRP[SL]|Rp[sl]|Mrp[sl]", col.name = "pct_ribo")

ae1_metadata <- ae1@meta.data %>% rownames_to_column() %>% as_tibble()
ae2_metadata <- ae2@meta.data %>% rownames_to_column() %>% as_tibble()

# ae1 nFeature_RNA vs nCount_RNA

p1 <-  ae1@meta.data %>% 
  ggplot(aes(x = nCount_RNA, 
             y = nFeature_RNA)) +
  ggtitle("QC Metrics (Sample: ae1): nFeature_RNA and nCount_RNA") + 
  geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
  geom_hline(yintercept=2500, 
             linetype="dashed", 
             color = "black", 
             size=1) + 
  geom_vline(xintercept = 10000, 
             linetype = "dashed", 
             color = "black", 
             size=1)
p1 <- ggMarginal(p1, type = "densigram", color="lightblue", fill="darkblue")
p1
ggsave(paste(cur_output_data_dir_eda, "ae1_nFeature_RNA_nCount_RNA_scatter.png", sep = ""), plot = p1)

# ae2 nFeature_RNA vs nCount_RNA

p1 <-  ae2@meta.data %>% 
  ggplot(aes(x = nCount_RNA, 
             y = nFeature_RNA)) +
  ggtitle("QC Metrics (Sample: ae2): nFeature_RNA and nCount_RNA") + 
  geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
  geom_hline(yintercept=2500, 
             linetype="dashed", 
             color = "black", 
             size=1) + 
  geom_vline(xintercept = 10000, 
             linetype = "dashed", 
             color = "black", 
             size=1)
p1 <- ggMarginal(p1, type = "densigram", color="lightblue", fill="darkblue")
p1
ggsave(paste(cur_output_data_dir_eda, "ae2_nFeature_RNA_nCount_RNA_scatter.png", sep = ""), plot = p1)

# ae1 pct_mito and pct_ribo plots

p2 <- ae1@meta.data %>% 
  ggplot(aes(x = pct_ribo, 
             y = pct_mito)) +
  ggtitle("QC Metrics (Sample: ae1): pct_mito and pct_ribo") + 
  geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
  geom_hline(yintercept = 20, linetype = "dashed", color = "black", size = 1) + 
  geom_vline(xintercept = 50, linetype="dashed", color = "black", size = 1)
p2 <- ggMarginal(p2, type = "densigram", color="lightblue", fill="darkblue")
p2
ggsave(paste(cur_output_data_dir_eda, "ae1_pct_mito_pct_ribo_scatter.png", sep = ""), plot = p2)

# ae2 pct_mito and pct_ribo plots

p2 <- ae2@meta.data %>% 
  ggplot(aes(x = pct_ribo, 
             y = pct_mito)) +
  ggtitle("QC Metrics (Sample: ae2): pct_mito and pct_ribo") + 
  geom_point(color = "blue", alpha = 0.1, size = 0.5) + 
  geom_hline(yintercept = 20, linetype = "dashed", color = "black", size = 1) + 
  geom_vline(xintercept = 50, linetype="dashed", color = "black", size = 1)
p2 <- ggMarginal(p2, type = "densigram", color="lightblue", fill="darkblue")
p2
ggsave(paste(cur_output_data_dir_eda, "ae2_pct_mito_pct_ribo_scatter.png", sep = ""), plot = p2)

# Generate histogram of RNA counts per cell 

# ae1

cell_count_hist <- qplot(ae1$nCount_RNA, geom="histogram", binwidth=500, main="QC Metrics (Sample: ae1): nCount_RNA Histogram", xlab = "Count depth (Counts/cell)", fill=I("darkblue"), col=I("black"))
cell_count_hist
ggsave(paste(cur_output_data_dir_eda, "ae1_nCount_RNA_histogram.png", sep = ""), plot = cell_count_hist)

# ae2

cell_count_hist <- qplot(ae2$nCount_RNA, geom="histogram", binwidth=500, main="QC Metrics (Sample: ae2): nCount_RNA Histogram", xlab = "Count depth (Counts/cell)", fill=I("darkblue"), col=I("black"))
cell_count_hist
ggsave(paste(cur_output_data_dir_eda, "ae2_nCount_RNA_histogram.png", sep = ""), plot = cell_count_hist)

#Generate histogram of feature counts per cell

# ae1

feature_count_hist <- qplot(ae1$nFeature_RNA, geom = "histogram", binwidth = 100, main = "QC Metrics (Sample: ae1): nFeature_RNA Histogram", xlab = "Number of genes per cell", fill = I("darkblue"), col = I("black"))
feature_count_hist
ggsave(paste(cur_output_data_dir_eda, "ae1_nFeature_RNA_histogram.png", sep = ""), plot = feature_count_hist)

# ae2

feature_count_hist <- qplot(ae2$nFeature_RNA, geom = "histogram", binwidth = 100, main = "QC Metrics (Sample: ae2): nFeature_RNA Histogram", xlab = "Number of genes per cell", fill = I("darkblue"), col = I("black"))
feature_count_hist
ggsave(paste(cur_output_data_dir_eda, "ae2_nFeature_RNA_histogram.png", sep = ""), plot = feature_count_hist)

# Filter data using QC metrics to prep for integration

ae1 <- subset(ae1, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & percent.mt < 25)
ae2 <- subset(ae2, subset = nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 25)

# Set tx conditions in dataset

ae1$STIM <- "PRE"
ae2$STIM <- "MID"

# Merge duplicates

ae.merge <- merge(x = ae1, y = ae2, add.cell.ids = c("ae1", "ae2"), project = "ae")

# Generate list from data

skin.list <- list(ae.merge)

# Normalize data and find variable features

for (i in 1:length(skin.list)) {
  skin.list[[i]] <- NormalizeData(skin.list[[i]], verbose = TRUE)
  skin.list[[i]] <- FindVariableFeatures(skin.list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = TRUE)
}

# Perform dataset integration

features <- SelectIntegrationFeatures(object.list = skin.list)

skin.list <- lapply(X = skin.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# TODO: How do we know to use 30 dims, 2k anchor features here?
# Note that using reciprocal PCA (RPCA) in the below function time significantly cuts down on runtime

anchors <- FindIntegrationAnchors(object.list = skin.list, anchor.features = 2000, reduction = "rpca", 
                                  dims = 1:30)

# Save anchors prior to integrating.
# Note that this is primarily being done because of memory issues.

saveRDS(anchors, paste(cur_output_data_dir, "anchors.rds", sep = ""))

# Load anchors (NOTE: only needed for memory issues)

anchors <- readRDS(paste(cur_output_data_dir, "anchors.rds", sep = ""))

# Integrate dataset

skin.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Save integrated dataset (unprocessed)

saveRDS(skin.integrated, paste(cur_output_data_dir, "skin_integrated_unprocessed.rds", sep = ""))

# Load unprocessed integrated dataset (NOTE: only needed for memory issues)

skin.integrated <- readRDS(paste(cur_output_data_dir, "skin_integrated_unprocessed.rds", sep = ""))

# Set default assay to integrated

DefaultAssay(skin.integrated) <- "integrated"

skin.integrated <- ScaleData(skin.integrated, verbose = FALSE)
skin.integrated <- RunPCA(skin.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(skin.integrated, ndims = 30)
skin.integrated <- RunUMAP(skin.integrated, dims = 1:30) # Default reduction is PCA
skin.integrated <- RunTSNE(skin.integrated, reduction = "pca", check_duplicates = FALSE, dims = 1:30)
skin.integrated <- FindNeighbors(skin.integrated, reduction = "pca", dims = 1:30)

# Visualize the top genes associated with PCs

pc_viz <- VizDimLoadings(skin.integrated, dims = 1:30, nfeatures = 10, reduction = "pca")
pc_viz
ggsave(paste(cur_output_data_dir, "PC_top_genes.png", sep = ""), plot = pc_viz, width = 40, height = 60, units = "cm")

# Create scree plot and manually determine elbow in order to know which PCs are most important

elbow_plot <- ElbowPlot(skin.integrated, ndims = 30) + ggtitle("Principal Component Elbow Plot")
ggsave(paste(cur_output_data_dir, "PC_elbow_plot.png", sep = ""), plot = elbow_plot)

skin.integrated <- FindClusters(skin.integrated, resolution = 0.35)

# Visualize cluster results

# Note that in the paper, they use the TSNE reduction for their DimPlot, not UMAP
cluster_plot <- DimPlot(skin.integrated, reduction = "tsne")
ggsave(paste(cur_output_data_dir, "Clusters.png", sep = ""), plot = cluster_plot)

# Cluster plot but split by tx condition

cluster_plot_split <- DimPlot(skin.integrated, reduction = "tsne", split.by = "STIM")
ggsave(paste(cur_output_data_dir, "Clusters_split.png", sep = ""), plot = cluster_plot_split)

# Save integrated dataset

saveRDS(skin.integrated, paste(cur_output_data_dir, "skin_integrated.rds", sep = ""))

# Load integrated dataset

skin.integrated <- readRDS(paste(cur_output_data_dir, "skin_integrated.rds", sep = ""))

# Get list of cluster IDs

clusters <- sort(unique(skin.integrated@meta.data$seurat_clusters))

# Print number of cells per cluster

table(Idents(skin.integrated))

# Print proportion of cells in each cluster

prop.table(table(Idents(skin.integrated)))

# Iterate through each cluster, find differentially expressed genes, and save markers to Excel file
# TODO: we need to run both FindConservedMarkers (for finding genes expressed in the same
# cluster and conserved across each tx condition, i.e. genes that are constant for a fixed cluster 
# and variable condition), and FindAllMarkers (for finding genes expressed in 
# a particular cluster overall vs other clusters, irrespective of tx condition, i.e.
# genes that are upregulated in a fixed cluster versus all other clusters
# when we do NOT split by tx condition). See table S3 for more info.

OUT <- createWorkbook()

# Switch default assay back to RNA for differential analysis after integration

DefaultAssay(skin.integrated) <- "RNA"

for(cluster in clusters){
  print(cluster)
  
  # Conserved markers, use FindConservedMarkers (same cluster, highly expressed across both tx)
  
  cur_markers <- FindConservedMarkers(
    skin.integrated,
    ident.1 = paste0(cluster),
    grouping.var = "STIM",
    verbose = TRUE,
    test.use = "MAST"
  ) %>%
    mutate(., gene = rownames(.))
  tab  <- paste0("Cluster_conserved_", cluster)
  addWorksheet(OUT, tab)
  writeData(OUT, sheet = tab, x = cur_markers)
}

saveWorkbook(OUT, paste(cur_output_data_dir, "FindConservedMarkers.xlsx", sep = ""), overwrite = TRUE)

# # All markers, use FindAllMarkers (same cluster, highly expressed for all tx combined)
# Meaning enriched in given cluster but upregulated only in certain tx condition
# # Note: this does not need to be in the loop, computes for all clusters automatically

# Differentially expressed markers for pre vs mid, use FindMarkers (same cluster, differentially expressed across tx)

skin.integrated$clusternumber.STIM <- paste(Idents(object = skin.integrated),sep = "_", skin.integrated$STIM)
skin.integrated$clusternumber <- Idents(object = skin.integrated)
Idents(object = skin.integrated) <- "clusternumber.STIM"

OUT3 <- createWorkbook()

for(cluster in clusters){
  if (is.element(paste0(cluster, "_MID"),unique(skin.integrated$clusternumber.STIM))
      & is.element(paste0(cluster, "_PRE"),unique(skin.integrated$clusternumber.STIM)))
    {
      print(cluster)
      cur_markers <- FindMarkers(
        skin.integrated, 
        ident.1 = paste0(cluster, "_MID"),
        ident.2 = paste0(cluster, "_PRE"),
        verbose = TRUE, 
        assay = "RNA", 
        slot = "data", 
        test.use = "MAST"
      ) %>%
        mutate(., gene = rownames(.))
      tab  <- paste0("Cluster_deg_", cluster)
      addWorksheet(OUT3, tab)
      writeData(OUT3, sheet = tab, x = cur_markers)
  }
}

saveWorkbook(OUT3, paste(cur_output_data_dir, "FindMarkers_pre_mid.xlsx", sep = ""), overwrite = TRUE)

table(skin.integrated$clusternumber.STIM) # Checking cell count in Seurat obj

# After identifying cell types by inspecting marker genes (outputted to FindConservedMarkers
# Excel file), update labels.
# FIXME the below fn call may need to be updated if the clusters selected are inaccurate.

Idents(object = skin.integrated) <- "clusternumber"
retained_clusters <- subset(x = skin.integrated, idents = c(0, 
                                                            1, 
                                                            2,
                                                            3,
                                                            4,
                                                            5,
                                                            6,
                                                            7,
                                                            8
)
)

DefaultAssay(retained_clusters) <- "RNA"

retained_clusters <- RenameIdents(retained_clusters, 
                                  '0' = "ILC2", 
                                  '1' = "Unknown_1", 
                                  '2' = "NK",
                                  '3' = "Unknown_3", #"dγδT", 
                                  '4' = "Unknown_4",
                                  '5' = "DETC",
                                  '6' = "Unknown_6", 
                                  '7' = "M/MdM", 
                                  '8' = "M/B"
)

labeled_clusters <- DimPlot(retained_clusters, reduction = "tsne", label = TRUE)
ggsave(paste(cur_output_data_dir, "clustered_labeled.png", sep = ""), plot = labeled_clusters)

cell_populations <- c("ILC2",
                      "Unknown_1",
                      "NK",
                      "Unknown_3",
                      "Unknown_4",
                      "DETC",
                      "Unknown_6",
                      "M/MdM",
                      "M/B"
)

# TODO: used same markers below as in existing Github code, update to assign dynamically?

markers.to.plot <- rev(c("Adgre1", "Itgam", "Cd68", "Lyz2", "Ms4a7", "C1qc", "Plac8", "Ly6c2", "Clec9a", "Xcr1", "Mgl2", "Cd209d",
                         "H2-Ab1", "Sirpa", "Fscn1", "Cacnb3","Apol7c", "Flt3", "Ccr7", "Cd207", "Epcam", "Ly75", "Itgax", "Itgae", 
                         "Cd3g", "Cd3e", "Thy1", "Trbc1", "Trdc", "Tcrg-C1", "Cd163l1", "X5830411N06Rik","Icos", "Trac", "Cd4", "Cd8a", "Il2ra", "Foxp3", "Gata3", "Rora", 
                         "Gzma", "Gzmb", "Klra7", "Klra8", "Ncr1", "Gata2", "Ms4a2", "Cpa3",	"Mcpt4", "Mcpt8",
                         "Cd14", "Csf3r", "S100a9", "S100a8", "Sell",
                         "FKBP5", "CTLA4", "IL7R", "GZMH", "CD8B", "CD8A", "GZMA", "PRF1",
                         "CD3D", "CD3G", "HAVCR2", "EOMES", "CD3E", "TMEM123", "KIT",
                         "PLAC8", "FTH1", "GATA2"))

conserved_marker_plot <- DotPlot(retained_clusters, features = rev(markers.to.plot), cols=c("lightgrey", "darkslateblue"), col.min = 0, col.max = 0.2, dot.min=0) + 
  ggtitle("Relative Gene Expression by Cell Type (All Tx Conditions)") + 
  theme(axis.text.x = element_text(size = 8, angle =90))+ scale_size(range = c(0, 6)) + 
  theme(legend.title = element_text(color = "black", size = 12))

ggsave(paste(cur_output_data_dir, "Relative Gene Expression by Cell Type (All Tx Conditions).png", sep = ""), 
       plot = conserved_marker_plot,
       width = 40, 
       height = 15, 
       units = "cm"
)

# Now we should plot DEGs between OXA vs EtOH tx conditions
# average log fold change > 0.5 and p.adj < 0.05

idx <- 1
FindMarkers_sheetnames <- getSheetNames(paste(cur_output_data_dir, "FindMarkers_pre_mid.xlsx", sep = ""))
for (sheet in FindMarkers_sheetnames) {
  cur_cluster_ident <- parse_number(sheet) #stri_sub(sheet,-1,-1)
  if (cur_cluster_ident %in% sort(unique(retained_clusters@meta.data$seurat_clusters))) {
    cur_markers <- read.xlsx(paste(cur_output_data_dir, "FindMarkers_pre_mid.xlsx", sep = ""), sheet = sheet)
    theme_set(theme_cowplot())
    cur_cell_population <- cell_populations[idx]
    genes_to_label <- pull(slice_min(cur_markers, n = 15, order_by = p_val_adj, with_ties = FALSE), var = gene)
    cur_cluster_cells <- subset(skin.integrated, idents = cur_cluster_ident)
    Idents(cur_cluster_cells) <- "STIM"
    avg_cur_cluster_cells <- log1p(AverageExpression(cur_cluster_cells, verbose = FALSE)$RNA)
    p1 <- ggplot(as.data.frame(avg_cur_cluster_cells), aes(MID, PRE)) + geom_point(color='darkblue') + ggtitle(paste0(cur_cell_population, ": Differentially Expressed Genes (PRE vs MID)"))
    p1 <- LabelPoints(plot = p1, points = genes_to_label)
    ggsave(paste(cur_output_data_dir, paste0("pre_vs_mid_deg_dotplot_cluster_", cur_cluster_ident, ".png"), sep = ""), 
           plot = p1)
  }
  idx <- idx + 1
}  

DefaultAssay(retained_clusters) <- "RNA"

# Create new column to store cell type

retained_clusters$celltype <- Idents(retained_clusters)

plots <- VlnPlot(retained_clusters, features = markers.to.plot, split.by = "STIM", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

idx <- 1
for (plot in plots) {
  cur_feature <- markers.to.plot[idx]
  ggsave(paste(cur_output_data_dir, paste0(cur_feature, " Expression by Cell Type.png"), sep = ""), 
         plot = plot)
  idx <- idx + 1
}