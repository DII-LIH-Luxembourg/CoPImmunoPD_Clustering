
###############################################################################
## Load RDS sub-type
###############################################################################
pht <- "m_m" # select the desired sub-type
# pht <- "p_p"
# pht <- "m_p"
# pht <- "p_m"

alldata.int <- readRDS(paste0("./data/sce_",pht,".rds"))

###############################################################################
## Sub-type m_p and p_m will have merged clusters
###############################################################################
if (pht == "p_m") {
  sel.clust = "CCA_snn_res.0.1"
  alldata.int <- SetIdent(alldata.int, value = sel.clust)
  
  new.cluster.ids <- c("Merged_0_1", "Merged_0_1", "Original_2")
  names(new.cluster.ids) <- levels(alldata.int)
  alldata.int <- RenameIdents(alldata.int, new.cluster.ids)
  
} else if (pht == "m_p") {
  sel.clust = "CCA_snn_res.0.1"
  alldata.int <- SetIdent(alldata.int, value = sel.clust)
  
  new.cluster.ids <- c("Merged_0_1", "Merged_0_1", "Original_2", "Original_3")
  names(new.cluster.ids) <- levels(alldata.int)
  alldata.int <- RenameIdents(alldata.int, new.cluster.ids)
} else if (pht == "p_p") {
  sel.clust = "CCA_snn_res.0.2"
  alldata.int <- SetIdent(alldata.int, value = sel.clust)
  
} else if (pht == "m_m") {
  sel.clust = "CCA_snn_res.0.1"
  alldata.int <- SetIdent(alldata.int, value = sel.clust)
  
}

###############################################################################
## Downsample the number of cells per identity class for UMAP visualization
###############################################################################
Idents(alldata.int) <- "sample"
head(alldata.int@meta.data)

t <- table(alldata.int@meta.data$sample)
sub_sce <- subset(x = alldata.int, downsample = min(t))

table(sub_sce@meta.data$sample)

clustree(alldata.int@meta.data, prefix = "CCA_snn_res.")
# select cluster resolution according to sub-type
Idents(alldata.int) <- sel.clust

mycolors <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2")
DimPlot(sub_sce, split.by = "condition", pt.size = 1, cols = mycolors, group.by = sel.clust)              
         

###############################################################################
## DGE
###############################################################################

# plot this clustering
DimPlot(alldata.int)

markers_genes <- FindAllMarkers(alldata.int, log2FC.threshold = 0.5, test.use = "wilcox",
                                min.pct = 0.3, 
                                only.pos = TRUE, assay = "RNA")
dim(markers_genes)

write.csv(markers_genes, paste0("./results/",
                                pht,"/marker_genes_",pht,"_resolution_",
                                sel.clust,".csv"))

markers_genes %>%
  group_by(cluster) %>%
  top_n(-25, p_val_adj) -> top25
top25

mypar(2, 4, mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  barplot(sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], F),
          horiz = T, las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
  abline(v = c(0, 0.25), lty = c(1, 2))
}

markers_genes %>%
  group_by(cluster) %>%
  top_n(-5, p_val_adj) -> top5

###############################################################################
## 
# topVlnPlot <-  markers_genes %>% group_by(cluster) %>% 
#   arrange(desc(avg_log2FC), .by_group = T) %>%
#   top_n(10, avg_log2FC)
# 
# VlnPlot(alldata.int, features = as.character(unique(topVlnPlot$gene)), ncol = 10, 
#         group.by = sel.clust, split.by = "condition",
#         assay = "RNA", pt.size = 0, cols = c("black", "red"))

###############################################################################
## display the top genes from the same top25 selection
## with custom threshold settings

topVlnPlot <-  top25 %>% group_by(cluster) %>% 
  arrange(desc(avg_log2FC), .by_group = T) %>%
  top_n(5, avg_log2FC)

alldata.int@meta.data$condition <- as.factor(alldata.int@meta.data$condition)
# set pt.size to zero if you do not want all the points to hide the violin
# shapes, or to a small value like 0.1
VlnPlot(alldata.int, features = as.character(unique(topVlnPlot$gene)), ncol = 5, 
        group.by = sel.clust, split.by = "condition",
        assay = "RNA", pt.size = 0, cols = c("black", "red"))

###############################################################################
## define features to be plotted
## adapt to the sub-type accordingly 
myfeatures <- c("THEMIS", "CAMK4", "GZMK", "PAG1", "COTL1",
                "IKZF2", "TYROBP", "TRIO", "TRDC", "KLRC2", "KLRC3", "TIGIT", "NCR1", "FCGR3A")
###############################################################################
# shapes, or to a small value like 0.1
VlnPlot(alldata.int, features = myfeatures, ncol = 5, split.by = "condition",
        assay = "RNA", pt.size = 0, cols = c("black", "red"))

###############################################################################
## Differential expression across conditions
###############################################################################
# select all cells in a cluster
for (l in levels(Idents(alldata.int))) {
  
  # select all cells in a cluster
  selCluster <- l
  
  cell_selection <- subset(alldata.int, cells = colnames(alldata.int)[alldata.int@meta.data[, sel.clust] ==
                                                                        selCluster])
  
  cell_selection <- SetIdent(cell_selection, value = "condition")
  
  DGE_cell_selection <- FindAllMarkers(cell_selection, log2FC.threshold = 0.2, 
                                       test.use = "wilcox", min.pct = 0.1, 
                                       only.pos = T,
                                       assay = "RNA")
  
  avg_expr <- AverageExpression(cell_selection, assays = "RNA")
  head(avg_expr$RNA)
  avg_expr <- as.data.frame(avg_expr$RNA)
  avg_expr$gene <- rownames(avg_expr)
  
  DGE_complete <- merge(x=DGE_cell_selection, y=avg_expr, by="gene")
  write.csv(DGE_complete, paste0("./results/",pht,"/DGE_",pht,"_cluster_",
                                 selCluster,".csv"))
}


###############################################################################
## Heatmap m_m
###############################################################################
glist <- c("PRF1", "HLA-DRA", "GZMB", "KLRD1", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
           "FCGR3A", "ITGB2", "PLCG2", "PRF1", "GZMB", "KLRD1", "ITGAL",
           "TNFAIP3", "PIK3R1", "JUNB", "NFKB1")

glist <- unique(glist)

heatmap_df <- data.frame()

sce_scale <- alldata.int
Idents(sce_scale) <- sel.clust

# Loop through conditions and clusters to create heatmaps
for (c in levels(as.factor(alldata.int@meta.data$condition))) {
  
  subset_obj <- subset(sce_scale, subset = condition == c)
  
  tmp <- AverageExpression(subset_obj)
  mydf <- as.data.frame(tmp$RNA)
  
  colnames(mydf) <- paste0(colnames(mydf), "_",c)
  mydf <- mydf[rownames(mydf) %in% glist,]
  
  # Append the z-scores to the heatmap data frame
  if (dim(heatmap_df)[2] == 0) {
    heatmap_df <- mydf  
  } else {
    heatmap_df <- cbind(heatmap_df, mydf)  
  }
}

# check the umap cluster colors
mycolors <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2")

my_sample_col <- data.frame(condition = as.factor(c("HC","PD","HC","PD")),
                            cluster = as.factor(c("0","0","1","1")))

# re-arrange the column order
heatmap_df <- heatmap_df[, c(1,3,2,4)]
row.names(my_sample_col) <- colnames(heatmap_df)

ann_colors = list(
  cluster = c("0"="#E69F00", "1"="#56B4E9"),
  condition = c("HC" = "#DDCC77","PD"="#CC79A7"))

# Set heatmap colors
col_palette <- colorRampPalette(c("blue", "white", "red"))(100)

pheatmap(heatmap_df,legend = T,annotation_legend = T,
         scale = "row", 
         cluster_rows = T,
         gaps_col = 2,
         # scale = "column",
         annotation_col = my_sample_col,
         annotation_colors = ann_colors,
         cluster_cols = FALSE, 
         color = col_palette, 
         show_colnames = FALSE)

###############################################################################
# length(glist)
VlnPlot(alldata.int, features = glist[22:34], ncol = 5, 
        group.by = sel.clust, split.by = "condition",
        assay = "RNA", pt.size = 0, cols = c("black", "red"))














