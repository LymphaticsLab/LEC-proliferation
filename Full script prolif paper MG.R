knitr::opts_chunk$set(message=FALSE, warning=FALSE, result='hold',fig.width=10, fig.height = 8)
library(SeuratObject)
library(Seurat)
suppressPackageStartupMessages(require(Matrix))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(gridExtra))
suppressPackageStartupMessages(require(pheatmap))
suppressPackageStartupMessages(require(scales))
suppressPackageStartupMessages(require(rlang))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(grid))
suppressPackageStartupMessages(require(topGO))
suppressPackageStartupMessages(require(org.Dr.eg.db))
suppressPackageStartupMessages(require(knitr))
datadir = "/Users/margl916/Downloads/v3_SCT/"
adata_5d <- readRDS(file.path(datadir,"data_5days.Rds"))
adata_7d <- readRDS(file.path(datadir,"data_7days.Rds"))
adata_36h <- readRDS(file.path(datadir,"data_36h.Rds"))
adata_48h <- readRDS(file.path(datadir,"data_48h.Rds"))


#CLEANING Fig3 and FigS3

wholedata_unclean <- readRDS(file.path(datadir,"seurat_integration_SCT-2.Rds"))

wholedata_unclean@active.ident = as.factor(wholedata_unclean$Timepoint)
DimPlot(wholedata_unclean, pt.size=1, label.size=10, cols = c("36h" = "#8c510a", "48h" = "#5ab4ac", "5days" = "#d8b365", "7days" = "#01665e")) & theme(text = element_text(face = "bold", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

wholedata_unclean@active.ident = as.factor(wholedata_unclean$seurat_clusters)
DimPlot(wholedata_unclean, pt.size=1, label.size=10, cols = c("0" = "#d8b365","1"="#5ab4ac", "2" = "dimgrey", "3" = "#c7eae5", "4" = "#f6e8c3", "5" = "black", "6" = "darkgrey", "7" = "#762a83", "8" = "#af8dc3", "9" = "#e7d4e8", "10" = "#d9f0d3", "11" = "#7fbf7b", "12" = "#1b7837", "13" = "#8c510a", "14" = "#01665e")) & theme(text = element_text(face = "bold", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

wholedata_unclean@active.ident = as.factor(wholedata_unclean$Part)
DimPlot(wholedata_unclean, pt.size=1, label.size=10, cols = c("Trunk" = "#d8b365", "Head" = "#5ab4ac", "Whole" = "dimgrey")) & theme(text = element_text(face = "bold", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

datadir2 = "/Users/margl916/Downloads/"
adata_5d_unclean <- readRDS(file.path(datadir2,"dataI.SCT.seurat_5days.Rds"))
adata_7d_unclean <- readRDS(file.path(datadir2,"dataI.SCT.seurat_7days.Rds"))
adata_48h_unclean <- readRDS(file.path(datadir2,"dataI.SCT.seurat_48h.Rds"))
adata_36h_unclean <- readRDS(file.path(datadir2,"dataI.SCT.seurat_36h.Rds"))

adata_5d_unclean@active.ident = as.factor(adata_5d_unclean$seurat_clusters)
adata_7d_unclean@active.ident = as.factor(adata_7d_unclean$seurat_clusters)
adata_48h_unclean@active.ident = as.factor(adata_48h_unclean$seurat_clusters)
adata_36h_unclean@active.ident = as.factor(adata_36h_unclean$seurat_clusters)

DimPlot(adata_5d_unclean, pt.size=1, label.size=10, cols = c("0" = "#d8b365","1"="#5ab4ac", "2" = "dimgrey", "3" = "#c7eae5", "4" = "#f6e8c3", "5" = "#302c2d", "6" = "darkgrey", "7" = "#762a83", "8" = "#af8dc3", "9" = "#e7d4e8", "10" = "#d9f0d3", "11" = "#7fbf7b", "12" = "#1b7837", "13" = "#8c510a", "14" = "#01665e", "15"="lightgrey")) & theme(text = element_text(face = "bold", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DimPlot(adata_7d_unclean, pt.size=1, label.size=10, cols = c("0" = "#d8b365","1"="#5ab4ac", "2" = "dimgrey", "3" = "#c7eae5", "4" = "#f6e8c3", "5" = "#302c2d", "6" = "darkgrey", "7" = "#762a83", "8" = "#af8dc3", "9" = "#e7d4e8", "10" = "#d9f0d3", "11" = "#7fbf7b")) & theme(text = element_text(face = "bold", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DimPlot(adata_48h_unclean, pt.size=1, label.size=10, cols = c("0" = "#d8b365","1"="#5ab4ac", "2" = "dimgrey", "3" = "#c7eae5", "4" = "#f6e8c3", "5" = "#302c2d", "6" = "darkgrey", "7" = "#762a83", "8" = "#af8dc3", "9" = "#e7d4e8", "10" = "#d9f0d3", "11" = "#7fbf7b", "12" = "#1b7837", "13" = "#8c510a", "14" = "#01665e", "15"="lightgrey", "16"="darkgreen", "17"="purple", "18"="turquoise")) & theme(text = element_text(face = "bold", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DimPlot(adata_36h_unclean, pt.size=1, label.size=10, cols = c("0" = "#d8b365","1"="#5ab4ac", "2" = "dimgrey", "3" = "#c7eae5", "4" = "#f6e8c3", "5" = "#302c2d", "6" = "darkgrey", "7" = "#762a83", "8" = "#af8dc3", "9" = "#e7d4e8", "10" = "#d9f0d3", "11" = "#7fbf7b", "12" = "#1b7837", "13" = "#8c510a", "14" = "#01665e")) & theme(text = element_text(face = "bold", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

library(Seurat)
library(ggplot2)
library(viridis)

outdir <- "/Users/margl916/Documents/Unclean_Dotplots"
dir.create(outdir, showWarnings = FALSE)

marker_genes <- rev(c(
  "myb", "gata2b", "runx1", "cebpa",
  "pdgfra", "vim", "twist1a",
  "gata1a", "hbbe3",
  "mpeg1.1", "spi1b"
))

make_dotplot <- function(seurat_obj, timepoint_label, filename_prefix) {
  
  DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj$seurat_clusters <- as.factor(seurat_obj$seurat_clusters)
  
  genes_present <- marker_genes[marker_genes %in% rownames(seurat_obj)]
  
  p <- DotPlot(
    seurat_obj,
    features = genes_present,
    group.by = "seurat_clusters",
    assay = "SCT",
    dot.scale = 5,
    scale = FALSE
  ) +
    coord_flip() +
    scale_color_gradientn(
      colours = viridis::viridis(20)
    ) +
    ggtitle(timepoint_label) +
    theme_classic() +
    theme(
      text = element_text(face = "oblique", size = 10),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(face = "italic", size = 10),
      axis.title = element_text(size = 10),
      axis.title.y.right = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      axis.line = element_line(size = 1),
      plot.title = element_text(hjust = 0.5, size = 18)
    )
  
  ggsave(
    file.path(outdir, paste0(filename_prefix, "_dotplot.pdf")),
    p,
    width = 7,
    height = 5
  )
  
  return(p)
}

p36 <- make_dotplot(adata_36h_unclean, "36 hpf", "36h_unclean")
p48 <- make_dotplot(adata_48h_unclean, "48 hpf", "48h_unclean")
p5  <- make_dotplot(adata_5d_unclean,  "5 dpf",  "5d_unclean")
p7  <- make_dotplot(adata_7d_unclean,  "7 dpf",  "7d_unclean")

p36
p48
p5
p7


#CLUSTER ANALYSIS Fig3 and FigS3

DimPlot(adata_36h, pt.size=2, label.size=20, cols = c("0" = "#109630","1"="#762a83", "2" = "#a4da84", "3" = "#fcf2a7", "4" = "#f49262", "5" = "#c5a9ce"))
DimPlot(adata_48h, pt.size=2, label.size=20, cols = c("0" = "#a4da84","1"="#762a83", "2" = "#109630", "3" = "#96c7dd", "4" = "#f49262"))
DimPlot(adata_5d, pt.size=2, label.size=20, cols = c("0" = "#109630","1"="#762a83", "2" = "#96c7dd", "3" = "#54859b", "4" = "darkblue", "5" = "#c5a9ce", "6" = "black", "7" = "#f49262", "8" = "#a4da84"))
DimPlot(adata_7d, pt.size=2, label.size=20, cols = c("0" = "#762a83","1"="#109630", "2" = "#c5a9ce", "3" = "black", "4" = "#f49262", "5" = "darkblue", "6" = "#96c7dd"))

DoHeatmap(adata_36h, features = c("fli1a", "lyve1b", "dab2", "mafbb", "mafba", "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "prox1a", "mrc1b", "tbx1", "slc7a7", "igfbp5b", "itih5", "dut", "cenpx", "esama", "cdh5", "kdrl", "plvapb", "hlx1", "esm1", "cxcr4a", "mmp25a", "slc2a1a", "pcna", "mki67", "kif23"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DoHeatmap(adata_48h, features = c("fli1a", "lyve1b", "dab2", "mafbb", "mafba", "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "prox1a", "mrc1b", "tbx1", "slc7a7", "igfbp5b", "itih5", "dut", "cenpx", "esama", "cdh5", "kdrl", "plvapb", "hlx1", "esm1", "cxcr4a", "mmp25a", "slc2a1a", "pcna", "mki67", "kif23"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DoHeatmap(adata_5d, features = c("fli1a", "lyve1b", "dab2", "mafbb", "mafba", "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "prox1a", "mrc1b", "tbx1", "slc7a7", "igfbp5b", "itih5", "dut", "cenpx", "esama", "cdh5", "kdrl", "plvapb", "hlx1", "esm1", "cxcr4a", "mmp25a", "slc2a1a", "pcna", "mki67", "kif23"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DoHeatmap(adata_7d, features = c("fli1a", "lyve1b", "dab2", "mafbb", "mafba", "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "prox1a", "mrc1b", "tbx1", "slc7a7", "igfbp5b", "itih5", "dut", "cenpx", "esama", "cdh5", "kdrl", "plvapb", "hlx1", "esm1", "cxcr4a", "mmp25a", "slc2a1a", "pcna", "mki67", "kif23"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

FeaturePlot(adata_36h, features = "cdh6", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))
FeaturePlot(adata_48h, features = "cdh6", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))
FeaturePlot(adata_5d, features = "cdh6", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))
FeaturePlot(adata_7d, features = "cdh6", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))

VlnPlot(adata_36h, features = c("cdh6"))
VlnPlot(adata_48h, features = c("cdh6"))
VlnPlot(adata_5d, features = c("cdh6"))
VlnPlot(adata_7d, features = c("cdh6"))

VlnPlot(adata_36h, features = c("prox1a"))
VlnPlot(adata_48h, features = c("prox1a"))
VlnPlot(adata_5d, features = c("prox1a"))
VlnPlot(adata_7d, features = c("prox1a"))


FeaturePlot(adata_36h, features = "cdh5", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))
FeaturePlot(adata_48h, features = "cdh5", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))
FeaturePlot(adata_5d, features = "cdh5", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))
FeaturePlot(adata_7d, features = "cdh5", pt.size=2, label.size=10, order = TRUE) + scale_color_gradientn(colours = viridis::viridis(12))

VlnPlot(adata_36h, features = c("cdh5"))
VlnPlot(adata_48h, features = c("cdh5"))
VlnPlot(adata_5d, features = c("cdh5"))
VlnPlot(adata_7d, features = c("cdh5"))

VlnPlot(adata_36h, features = c("kdrl"))
VlnPlot(adata_48h, features = c("kdrl"))
VlnPlot(adata_5d, features = c("kdrl"))
VlnPlot(adata_7d, features = c("kdrl"))


source("/Users/margl916/Desktop/shiny/functions.R")
library(ggplot2)

adata_36h@active.ident = as.factor(adata_36h$seurat_clusters)
data36h_cluster_0_1 <- subset(x = adata_36h, idents = c("0", "1"))

plotClusters_marker(data36h_cluster_0_1, "prox1a") + geom_bar(stat = "identity") + scale_fill_manual(values = c("0" = "#109630","1"="#762a83"))
plotClusters_marker(data36h_cluster_0_1, "cdh6") + geom_bar(stat = "identity") + scale_fill_manual(values = c("0" = "#109630","1"="#762a83"))


adata_48h@active.ident = as.factor(adata_48h$seurat_clusters)
data48h_cluster_2_1 <- subset(x = adata_48h, idents = c("2", "1"))

plotClusters_marker(data48h_cluster_2_1, "prox1a") + geom_bar(stat = "identity") + scale_fill_manual(values = c("2" = "#109630","1"="#762a83"))
plotClusters_marker(data48h_cluster_2_1, "cdh6") + geom_bar(stat = "identity") + scale_fill_manual(values = c("2" = "#109630","1"="#762a83"))


adata_5d@active.ident = as.factor(adata_5d$seurat_clusters)
data5d_cluster_0_1 <- subset(x = adata_5d, idents = c("0", "1"))

plotClusters_marker(data5d_cluster_0_1, "prox1a") + geom_bar(stat = "identity") + scale_fill_manual(values = c("0" = "#109630","1"="#762a83"))
plotClusters_marker(data5d_cluster_0_1, "cdh6") + geom_bar(stat = "identity") + scale_fill_manual(values = c("0" = "#109630","1"="#762a83"))


adata_7d@active.ident = as.factor(adata_7d$seurat_clusters)
data7d_cluster_1_0 <- subset(x = adata_7d, idents = c("1", "0"))

plotClusters_marker(data7d_cluster_1_0, "prox1a") + geom_bar(stat = "identity") + scale_fill_manual(values = c("1" = "#109630","0"="#762a83"))
plotClusters_marker(data7d_cluster_1_0, "cdh6") + geom_bar(stat = "identity") + scale_fill_manual(values = c("1" = "#109630","0"="#762a83"))


# DEG cluster analysis

library(Seurat)
library(dplyr)
library(openxlsx)

# make output folder
outdir <- "/Users/margl916/Documents/DEG_lists_all_clusters"
dir.create(outdir, showWarnings = FALSE)

get_DEG_lists_all_clusters <- function(seurat_obj, timepoint_name) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  Idents(seurat_obj) <- as.factor(seurat_obj$seurat_clusters)
  
  # find markers for all clusters
  markers <- FindAllMarkers(
    object = seurat_obj,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25
  )
  
  # keep significant DEGs
  markers_sig <- markers %>%
    dplyr::filter(p_val_adj < 0.01) %>%
    dplyr::arrange(cluster, desc(avg_log2FC))
  
  markers_sig$timepoint <- timepoint_name
  
  # write one gene-symbol list per cluster
  clusters <- sort(unique(markers_sig$cluster))
  
  for (cl in clusters) {
    
    genes_cluster <- markers_sig %>%
      dplyr::filter(cluster == cl) %>%
      dplyr::pull(gene) %>%
      unique()
    
    write.table(
      genes_cluster,
      file = file.path(
        outdir,
        paste0("DEG_", timepoint_name, "_adjpval01_FC_cluster", cl, "_list.txt")
      ),
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  return(markers_sig)
}

DEG_36h <- get_DEG_lists_all_clusters(adata_36h, "36h")
DEG_48h <- get_DEG_lists_all_clusters(adata_48h, "48h")
DEG_5d  <- get_DEG_lists_all_clusters(adata_5d,  "5d")
DEG_7d  <- get_DEG_lists_all_clusters(adata_7d,  "7d")

DEG_all <- dplyr::bind_rows(
  DEG_36h,
  DEG_48h,
  DEG_5d,
  DEG_7d
)

# export full DEG tables as Excel
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "36h")
openxlsx::writeData(wb, "36h", DEG_36h)

openxlsx::addWorksheet(wb, "48h")
openxlsx::writeData(wb, "48h", DEG_48h)

openxlsx::addWorksheet(wb, "5d")
openxlsx::writeData(wb, "5d", DEG_5d)

openxlsx::addWorksheet(wb, "7d")
openxlsx::writeData(wb, "7d", DEG_7d)

openxlsx::addWorksheet(wb, "all_timepoints")
openxlsx::writeData(wb, "all_timepoints", DEG_all)

openxlsx::saveWorkbook(
  wb,
  file = file.path(outdir, "DEG_all_clusters_all_timepoints.xlsx"),
  overwrite = TRUE
)


# GO BP analysis
# convert zebrafish gene symbols to Entrez gene IDs

library(org.Dr.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(viridis)

keytypes(org.Dr.eg.db)

# load DEG gene-symbol list (same for other DEG lists)
DEG_Multi36_adjpval01_FC_cluster0_list <- scan(
  "/Users/margl916/Documents/DEG_Multi36_adjpval01_FC_cluster0_list.txt",
  character(),
  quote = ""
)

t <- c(DEG_Multi36_adjpval01_FC_cluster0_list)

# convert SYMBOL to Entrez ID
et <- bitr(
  t,
  fromType = "SYMBOL",
  toType = c("ENTREZID", "ALIAS", "GENENAME"),
  OrgDb = org.Dr.eg.db
)

head(et)

write.csv(
  et,
  "/Users/margl916/Documents/DEG_Multi36_cluster0_gene_conversion.csv",
  row.names = FALSE
)

# remove duplicate Entrez IDs
entrezID_et_single <- et %>%
  dplyr::select(ENTREZID) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

write.csv(
  entrezID_et_single,
  "/Users/margl916/Documents/DEG_Multi36_cluster0_entrezID_single.csv",
  row.names = FALSE
)

# make clean Entrez vector for enrichGO
entrezID_et_single_list <- entrezID_et_single$ENTREZID

# GO BP enrichment
ans.go <- enrichGO(
  gene          = entrezID_et_single_list,
  ont           = "BP",
  OrgDb         = org.Dr.eg.db,
  keyType       = "ENTREZID",
  readable      = TRUE,
  pvalueCutoff  = 0.5,
  pAdjustMethod = "BH"
)

ans.go

tab.go <- as.data.frame(ans.go)
View(tab.go)

write.csv(
  tab.go,
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_enrichGO.csv",
  row.names = FALSE
)

# visualization: top 20
dotplot(ans.go, showCategory = 20) +
  ggtitle("GO BP") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  viridis::scale_fill_viridis(direction = -1)

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_dotplot_top20.pdf",
  width = 6,
  height = 6
)

# visualization: top 10 by qvalue
dotplot(ans.go, showCategory = 10, orderBy = "qvalue") +
  ggtitle("GO BP") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  viridis::scale_fill_viridis(direction = -1) +
  scale_y_discrete(limits = rev)

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_dotplot_top10_qvalue.pdf",
  width = 6,
  height = 6
)

# selected vascular / proliferation GO terms
selected_GO_terms <- c(
  "blood vessel development",
  "angiogenesis",
  "endothelial cell differentiation",
  "venous blood vessel development",
  "lymph vessel development",
  "lymphangiogenesis",
  "sprouting angiogenesis",
  "endothelial cell migration",
  "vascular endothelial growth factor receptor signaling pathway",
  "branching involved in blood vessel morphogenesis",
  "blood vessel lumenization",
  "endothelial cell proliferation",
  "DNA replication",
  "mitotic cell cycle process",
  "cell population proliferation",
  "otolith development"
)

dotplot(ans.go, showCategory = selected_GO_terms, orderBy = "qvalue") +
  ggtitle("GO BP") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  viridis::scale_fill_viridis(direction = -1) +
  scale_y_discrete(limits = rev)

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_dotplot_selected_terms.pdf",
  width = 8,
  height = 3.54
)

# vascular-focused terms
vascular_GO_terms <- c(
  "blood vessel development",
  "angiogenesis",
  "sprouting angiogenesis",
  "response to growth factor",
  "regulation of catalytic activity",
  "regulation of kinase activity",
  "branching involved in blood vessel morphogenesis",
  "blood vessel lumenization",
  "endothelial cell proliferation",
  "mitotic cell cycle process"
)

dotplot(ans.go, showCategory = vascular_GO_terms) +
  ggtitle("GO BP") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  viridis::scale_fill_viridis(direction = -1)

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_dotplot_vascular_terms.pdf",
  width = 6,
  height = 4
)

# cell-cycle-focused terms
cellcycle_GO_terms <- c(
  "DNA replication",
  "chromosome organization",
  "mitotic cell cycle process",
  "nuclear division",
  "DNA recombination",
  "chromatin organization",
  "chromatin remodeling",
  "regulation of cell cycle process"
)

dotplot(ans.go, showCategory = cellcycle_GO_terms) +
  ggtitle("GO BP") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  ) +
  viridis::scale_fill_viridis(direction = -1)

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_dotplot_cellcycle_terms.pdf",
  width = 5,
  height = 3.5
)

# cnetplot
cnetplot(
  ans.go,
  showCategory = c("blood vessel development", "lymphangiogenesis"),
  circular = TRUE,
  colorEdge = TRUE
) +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 10)
  )

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_cnetplot.pdf",
  width = 7,
  height = 7
)

# upset plot
upsetplot(ans.go, showCategory = 20) +
  ggtitle("GO BP")

ggsave(
  "/Users/margl916/Documents/DEG_Multi36_cluster0_GO_BP_upsetplot.pdf",
  width = 7,
  height = 5
)



#Prolif ANALYSIS Fig4 and FigS4

adata_36h@active.ident = as.factor(adata_36h$Phase)
DimPlot(adata_36h, pt.size=2, label.size=10, cols = c("G1" = "lightgrey","G2M"="#5ab4ac", "S" = "#01665e"), order = TRUE) & theme(text = element_text(face = "bold", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

adata_48h@active.ident = as.factor(adata_48h$Phase)
DimPlot(adata_48h, pt.size=2, label.size=10, cols = c("G1" = "lightgrey","G2M"="#5ab4ac", "S" = "#01665e"), order = TRUE) & theme(text = element_text(face = "bold", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

adata_5d@active.ident = as.factor(adata_5d$Phase)
DimPlot(adata_5d, pt.size=2, label.size=10, cols = c("G1" = "lightgrey","G2M"="#5ab4ac", "S" = "#01665e"), order = TRUE) & theme(text = element_text(face = "bold", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

adata_7d@active.ident = as.factor(adata_7d$Phase)
DimPlot(adata_7d, pt.size=2, label.size=10, cols = c("G1" = "lightgrey","G2M"="#5ab4ac", "S" = "#01665e"), order = TRUE) & theme(text = element_text(face = "bold", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 



#Vulcano Plot LEC Prolif for all 48hpf cluster LECs:

adata_48h@active.ident=as.factor(adata_48h$seurat_clusters)
data48h_clusterLECs <- subset(x = adata_48h, idents = c("0", "2", "4"))
data48h_clusterLECs@active.ident=as.factor(data48h_clusterLECs$CellCycle)
deMarkersLECProl48 = FindMarkers(data48h_clusterLECs, assay = "RNA", ident.1 = "G1", ident.2 = "prol", min.pct = 0.5)  

ggplot(data = deMarkersLECProl48, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') + 
  geom_point()

deMarkersLECProl48$diffexpressed <- "NO"
deMarkersLECProl48$diffexpressed[deMarkersLECProl48$avg_log2FC > 0.5 & deMarkersLECProl48$p_val < 0.1] <- "DOWN"
deMarkersLECProl48$diffexpressed[deMarkersLECProl48$avg_log2FC < -0.5 & deMarkersLECProl48$p_val < 0.1] <- "UP"

head(deMarkersLECProl48[order(deMarkersLECProl48$p_val_adj) & deMarkersLECProl48$diffexpressed == 'DOWN', ])
deMarkersLECProl48$delabel <- row.names(deMarkersLECProl48)

ggplot(data = deMarkersLECProl48,  aes(x = avg_log2FC, y = -log10(p_val), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("lightgrey", "#737373",  "#01665e"), 
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
  geom_text()


library(Seurat)
library(ggplot2)
library(dplyr)
library(openxlsx)

outdir <- "/Users/margl916/Documents/LEC_Prolif_Volcano"
dir.create(outdir, showWarnings = FALSE)

# define LEC subsets

adata_36h@active.ident=as.factor(adata_36h$seurat_clusters)
adata_48h@active.ident=as.factor(adata_48h$seurat_clusters)
adata_5d@active.ident=as.factor(adata_5d$seurat_clusters)
adata_7d@active.ident=as.factor(adata_7d$seurat_clusters)

allLECcluster36 <- subset(adata_36h, idents = c("0", "2", "4"))
allLECcluster48 <- subset(adata_48h, idents = c("0", "2", "4"))
allLECcluster5  <- subset(adata_5d,  idents = c("0", "7", "8"))
allLECcluster7  <- subset(adata_7d,  idents = c("1", "4"))

run_LEC_prolif_DE <- function(seurat_obj, timepoint_label) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj@active.ident <- as.factor(seurat_obj$CellCycle)
  
  deMarkers <- FindMarkers(
    seurat_obj,
    assay = "RNA",
    ident.1 = "G1",
    ident.2 = "prol",
    test.use = "wilcox",
    min.pct = 0.5,
    logfc.threshold = 0
  )
  
  deMarkers$gene <- rownames(deMarkers)
  deMarkers$timepoint <- timepoint_label
  deMarkers$comparison <- "G1_vs_prol"
  
  deMarkers$diffexpressed <- "NO"
  
  deMarkers$diffexpressed[
    deMarkers$avg_log2FC > 0.5 &
      deMarkers$p_val < 0.1
  ] <- "DOWN"
  
  deMarkers$diffexpressed[
    deMarkers$avg_log2FC < -0.5 &
      deMarkers$p_val < 0.1
  ] <- "UP"
  
  deMarkers$delabel <- deMarkers$gene
  
  deMarkers <- deMarkers %>%
    dplyr::select(
      timepoint,
      comparison,
      gene,
      p_val,
      avg_log2FC,
      pct.1,
      pct.2,
      p_val_adj,
      diffexpressed,
      delabel
    )
  
  return(deMarkers)
}

plot_volcano <- function(deMarkers, title_text) {
  
  ggplot(
    data = deMarkers,
    aes(
      x = avg_log2FC,
      y = -log10(p_val),
      col = diffexpressed,
      label = delabel
    )
  ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      col = "gray",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.1),
      col = "gray",
      linetype = "dashed"
    ) +
    geom_point(size = 2) +
    geom_text(size = 2.5, check_overlap = TRUE) +
    scale_color_manual(
      values = c(
        "DOWN" = "lightgrey",
        "NO" = "#737373",
        "UP" = "#01665e"
      ),
      labels = c(
        "DOWN" = "Downregulated",
        "NO" = "Not significant",
        "UP" = "Upregulated"
      )
    ) +
    ggtitle(title_text) +
    theme_classic() +
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      legend.title = element_blank()
    )
}

# run DE
deMarkersLECProl36 <- run_LEC_prolif_DE(
  allLECcluster36,
  "36 hpf"
)

deMarkersLECProl48 <- run_LEC_prolif_DE(
  allLECcluster48,
  "48 hpf"
)

deMarkersLECProl5 <- run_LEC_prolif_DE(
  allLECcluster5,
  "5 dpf"
)

deMarkersLECProl7 <- run_LEC_prolif_DE(
  allLECcluster7,
  "7 dpf"
)

# volcano plots
p36 <- plot_volcano(
  deMarkersLECProl36,
  "36 hpf LEC G1 vs prol"
)

p48 <- plot_volcano(
  deMarkersLECProl48,
  "48 hpf LEC G1 vs prol"
)

p5 <- plot_volcano(
  deMarkersLECProl5,
  "5 dpf LEC G1 vs prol"
)

p7 <- plot_volcano(
  deMarkersLECProl7,
  "7 dpf LEC G1 vs prol"
)

p36
p48
p5
p7

# save plots
ggsave(
  file.path(outdir, "Volcano_LEC_Prolif_36h.pdf"),
  p36,
  width = 6,
  height = 5
)

ggsave(
  file.path(outdir, "Volcano_LEC_Prolif_48h.pdf"),
  p48,
  width = 6,
  height = 5
)

ggsave(
  file.path(outdir, "Volcano_LEC_Prolif_5d.pdf"),
  p5,
  width = 6,
  height = 5
)

ggsave(
  file.path(outdir, "Volcano_LEC_Prolif_7d.pdf"),
  p7,
  width = 6,
  height = 5
)


outdir <- "/Users/margl916/Documents/Main_LEC_Prolif_Volcano"
dir.create(outdir, showWarnings = FALSE)

# define main LEC clusters: 

adata_36h@active.ident <- as.factor(adata_36h$seurat_clusters)
adata_48h@active.ident <- as.factor(adata_48h$seurat_clusters)
adata_5d@active.ident  <- as.factor(adata_5d$seurat_clusters)
adata_7d@active.ident  <- as.factor(adata_7d$seurat_clusters)

mainLECcluster36 <- subset(adata_36h, idents = c("0"))
mainLECcluster48 <- subset(adata_48h, idents = c("2"))
mainLECcluster5  <- subset(adata_5d,  idents = c("0"))
mainLECcluster7  <- subset(adata_7d,  idents = c("1"))

run_LEC_prolif_DE <- function(seurat_obj, timepoint_label) {
  
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj@active.ident <- as.factor(seurat_obj$CellCycle)
  
  deMarkers <- FindMarkers(
    seurat_obj,
    assay = "RNA",
    ident.1 = "G1",
    ident.2 = "prol",
    test.use = "wilcox",
    min.pct = 0.5,
    logfc.threshold = 0
  )
  
  deMarkers$gene <- rownames(deMarkers)
  deMarkers$timepoint <- timepoint_label
  deMarkers$comparison <- "G1_vs_prol_main_LEC_cluster_0_1"
  
  deMarkers$diffexpressed <- "NO"
  
  deMarkers$diffexpressed[
    deMarkers$avg_log2FC > 0.5 &
      deMarkers$p_val < 0.1
  ] <- "DOWN"
  
  deMarkers$diffexpressed[
    deMarkers$avg_log2FC < -0.5 &
      deMarkers$p_val < 0.1
  ] <- "UP"
  
  deMarkers$delabel <- deMarkers$gene
  
  deMarkers <- deMarkers %>%
    dplyr::select(
      timepoint,
      comparison,
      gene,
      p_val,
      avg_log2FC,
      pct.1,
      pct.2,
      p_val_adj,
      diffexpressed,
      delabel
    )
  
  return(deMarkers)
}

plot_volcano <- function(deMarkers, title_text) {
  
  ggplot(
    data = deMarkers,
    aes(
      x = avg_log2FC,
      y = -log10(p_val),
      col = diffexpressed,
      label = delabel
    )
  ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      col = "gray",
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.1),
      col = "gray",
      linetype = "dashed"
    ) +
    geom_point(size = 2) +
    geom_text(size = 2.5, check_overlap = TRUE) +
    scale_color_manual(
      values = c(
        "DOWN" = "lightgrey",
        "NO" = "#737373",
        "UP" = "#01665e"
      ),
      labels = c(
        "DOWN" = "Downregulated",
        "NO" = "Not significant",
        "UP" = "Upregulated"
      )
    ) +
    ggtitle(title_text) +
    theme_classic() +
    theme(
      text = element_text(size = 10),
      axis.text = element_text(size = 10),
      legend.title = element_blank()
    )
}

# run DE
deMarkersMainLECProl36 <- run_LEC_prolif_DE(mainLECcluster36, "36 hpf")
deMarkersMainLECProl48 <- run_LEC_prolif_DE(mainLECcluster48, "48 hpf")
deMarkersMainLECProl5  <- run_LEC_prolif_DE(mainLECcluster5,  "5 dpf")
deMarkersMainLECProl7  <- run_LEC_prolif_DE(mainLECcluster7,  "7 dpf")

# volcano plots
p36 <- plot_volcano(deMarkersMainLECProl36, "36 hpf main LEC cluster 0/1 G1 vs prol")
p48 <- plot_volcano(deMarkersMainLECProl48, "48 hpf main LEC cluster 0/1 G1 vs prol")
p5  <- plot_volcano(deMarkersMainLECProl5,  "5 dpf main LEC cluster 0/1 G1 vs prol")
p7  <- plot_volcano(deMarkersMainLECProl7,  "7 dpf main LEC cluster 0/1 G1 vs prol")

p36
p48
p5
p7

# save plots
ggsave(file.path(outdir, "Volcano_Main_LEC_Cluster_0_1_Prolif_36h.pdf"), p36, width = 6, height = 5)
ggsave(file.path(outdir, "Volcano_Main_LEC_Cluster_0_1_Prolif_48h.pdf"), p48, width = 6, height = 5)
ggsave(file.path(outdir, "Volcano_Main_LEC_Cluster_0_1_Prolif_5d.pdf"),  p5,  width = 6, height = 5)
ggsave(file.path(outdir, "Volcano_Main_LEC_Cluster_0_1_Prolif_7d.pdf"),  p7,  width = 6, height = 5)

# combine source data
DE_all_main_LEC_prolif <- dplyr::bind_rows(
  deMarkersMainLECProl36,
  deMarkersMainLECProl48,
  deMarkersMainLECProl5,
  deMarkersMainLECProl7
)

# export Excel source data
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "36h")
openxlsx::writeData(wb, "36h", deMarkersMainLECProl36)

openxlsx::addWorksheet(wb, "48h")
openxlsx::writeData(wb, "48h", deMarkersMainLECProl48)

openxlsx::addWorksheet(wb, "5d")
openxlsx::writeData(wb, "5d", deMarkersMainLECProl5)

openxlsx::addWorksheet(wb, "7d")
openxlsx::writeData(wb, "7d", deMarkersMainLECProl7)

openxlsx::addWorksheet(wb, "all_timepoints")
openxlsx::writeData(wb, "all_timepoints", DE_all_main_LEC_prolif)

openxlsx::saveWorkbook(
  wb,
  file = file.path(
    outdir,
    "Main_LEC_Cluster_0_1_Prolif_Volcano_source_data_all_timepoints.xlsx"
  ),
  overwrite = TRUE
)

#Combining cluster then Prol vs G1 Four-way single cell analysis

#48 hpf:

allLECcluster <- subset(adata_48h, idents = c("0", "2", "4"))
allVECcluster <- subset(adata_48h, idents = c("1", "3"))

allLECcluster@active.ident=as.factor(allLECcluster$CellCycle)
allVECcluster@active.ident=as.factor(allVECcluster$CellCycle)
allLECprol48 <- subset(allLECcluster, ident = "prol")
allLECG148 <- subset(allLECcluster, ident = "G1")
allVECprol48 <- subset(allVECcluster, ident = "prol")
allVECG148 <- subset(allVECcluster, ident = "G1")

ccprolL = rep("48h-prol-LEC", ncol(allLECprol48))
ccG1L = rep("48h-G1-LEC", ncol(allLECG148))
ccprolV = rep("48h-prol-VEC", ncol(allVECprol48))
ccG1V = rep("48h-G1-VEC", ncol(allVECG148))

allLECprol48 = AddMetaData(allLECprol48, factor(ccprolL), col.name = "CellCycle")
allLECG148 = AddMetaData(allLECG148, factor(ccG1L), col.name = "CellCycle")
allVECprol48 = AddMetaData(allVECprol48, factor(ccprolV), col.name = "CellCycle")
allVECG148 = AddMetaData(allVECG148, factor(ccG1V), col.name = "CellCycle")

merge_allLECprol48_allLECG148 = merge(allLECprol48, allLECG148)
merge_allVECprol48_allVECG148 = merge(allVECprol48, allVECG148)
merge_allLEC_allVEC_48 = merge(merge_allLECprol48_allLECG148, merge_allVECprol48_allVECG148)

merge_allLEC_allVEC_48@active.ident=as.factor(merge_allLEC_allVEC_48$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_48, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_48@active.ident=as.factor(merge_allLEC_allVEC_48$"CellCycle")
DoHeatmap(merge_allLEC_allVEC_48, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

FWDEG_merge_allLEC_allVEC_48 <- FindAllMarkers(merge_allLEC_allVEC_48, assay = "RNA")
write.csv(FWDEG_merge_allLEC_allVEC_48, "/Users/margl916/Documents/FWDEG_merge_allLEC_allVEC_48.csv", row.names=TRUE)

#7 dpf:

allLECcluster7 <- subset(adata_7d, idents = c("1", "4"))
allVECcluster7 <- subset(adata_7d, idents = c("0", "2", "5", "6"))

allLECcluster7@active.ident=as.factor(allLECcluster7$CellCycle)
allVECcluster7@active.ident=as.factor(allVECcluster7$CellCycle)
allLECprol7 <- subset(allLECcluster7, ident = "prol")
allLECG17 <- subset(allLECcluster7, ident = "G1")
allVECprol7 <- subset(allVECcluster7, ident = "prol")
allVECG17 <- subset(allVECcluster7, ident = "G1")

ccprolL7 = rep("7d-prol-LEC", ncol(allLECprol7))
ccG1L7 = rep("7d-G1-LEC", ncol(allLECG17))
ccprolV7 = rep("7d-prol-VEC", ncol(allVECprol7))
ccG1V7 = rep("7d-G1-VEC", ncol(allVECG17))

allLECprol7 = AddMetaData(allLECprol7, factor(ccprolL7), col.name = "CellCycle")
allLECG17 = AddMetaData(allLECG17, factor(ccG1L7), col.name = "CellCycle")
allVECprol7 = AddMetaData(allVECprol7, factor(ccprolV7), col.name = "CellCycle")
allVECG17 = AddMetaData(allVECG17, factor(ccG1V7), col.name = "CellCycle")

merge_allLECprol7_allLECG17 = merge(allLECprol7, allLECG17)
merge_allVECprol7_allVECG17 = merge(allVECprol7, allVECG17)
merge_allLEC_allVEC_7 = merge(merge_allLECprol7_allLECG17, merge_allVECprol7_allVECG17)

merge_allLEC_allVEC_7@active.ident=as.factor(merge_allLEC_allVEC_7$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_7, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_7@active.ident=as.factor(merge_allLEC_allVEC_7$"CellCycle")
DoHeatmap(merge_allLEC_allVEC_7, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

FWDEG_merge_allLEC_allVEC_7 <- FindAllMarkers(merge_allLEC_allVEC_7, assay = "RNA")
write.csv(FWDEG_merge_allLEC_allVEC_7, "/Users/margl916/Documents/FWDEG_merge_allLEC_allVEC_7.csv", row.names=TRUE)

#DEG filtering criteria
#adjusted p-value < 1 and expression rate > 50% of proliferating LECs
#additional stringency criteria log2FC > +/- 0.5 and exclusion of mitochondrial and histone genes

DoHeatmap(
  merge_allLEC_allVEC_48,
  features = c(
    "fli1a", "lyve1b", "dab2", "mrc1a", "mafbb", "mafba",
    "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "cldn11a",
    "prox1a", "mrc1b", "tbx1",
    "cdh5", "kdrl", "esama",
    "pcna", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7",
    "cdk2", "ccnd1", "cdk7", "mad2l2", "ccne1", "ccne2",
    "cdc6", "top2a", "birc5a", "mad2l1", "mad2l1bp",
    "ccna2", "ccnb1", "ccnb3", "cdk1", "cdc20",
    "aurka", "aurkb", "ccnb2", "mki67", "cdk6"
  ),
  draw.lines = TRUE,
  lines.width = 6,
  slot = "data"
) +
  scale_fill_gradientn(
    colors = viridis::viridis(20),
    na.value = "white"
  ) &
  theme(
    text = element_text(face = "oblique", size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.title.y.right = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.line = element_line(size = 1)
  )


DoHeatmap(
  merge_allLEC_allVEC_7,
  features = c(
    "fli1a", "lyve1b", "dab2", "mrc1a", "mafbb", "mafba",
    "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "cldn11a",
    "prox1a", "mrc1b", "tbx1",
    "cdh5", "kdrl", "esama",
    "pcna", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7",
    "cdk2", "ccnd1", "cdk7", "mad2l2", "ccne1", "ccne2",
    "cdc6", "top2a", "birc5a", "mad2l1", "mad2l1bp",
    "ccna2", "ccnb1", "ccnb3", "cdk1", "cdc20",
    "aurka", "aurkb", "ccnb2", "mki67", "cdk6"
  ),
  draw.lines = TRUE,
  lines.width = 6,
  slot = "data"
) +
  scale_fill_gradientn(
    colors = viridis::viridis(20),
    na.value = "white"
  ) &
  theme(
    text = element_text(face = "oblique", size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.title.y.right = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.line = element_line(size = 1)
  )

library(DoMultiBarHeatmap)

smerge_allLEC_allVEC_48 <- ScaleData(merge_allLEC_allVEC_48)
DoMultiBarHeatmap(smerge_allLEC_allVEC_48, assay = "RNA",  features = c(
  "fli1a", "lyve1b", "dab2", "mrc1a", "mafbb", "mafba",
  "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "cldn11a",
  "prox1a", "mrc1b", "tbx1",
  "cdh5", "kdrl", "esama",
  "pcna", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7",
  "cdk2", "ccnd1", "cdk7", "mad2l2", "ccne1", "ccne2",
  "cdc6", "top2a", "birc5a", "mad2l1", "mad2l1bp",
  "ccna2", "ccnb1", "ccnb3", "cdk1", "cdc20",
  "aurka", "aurkb", "ccnb2", "mki67", "cdk6"
), group.by="CellCycle", additional.group.by = "Phase", draw.lines = TRUE, lines.width = 6, size =3, group.bar.height = 0.05, angle = 45) + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

smerge_allLEC_allVEC_7 <- ScaleData(merge_allLEC_allVEC_7)
DoMultiBarHeatmap(smerge_allLEC_allVEC_7, assay = "RNA",  features = c(
  "fli1a", "lyve1b", "dab2", "mrc1a", "mafbb", "mafba",
  "stab2", "stab1", "hapln3", "cdh6", "cldn11b", "cldn11a",
  "prox1a", "mrc1b", "tbx1",
  "cdh5", "kdrl", "esama",
  "pcna", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7",
  "cdk2", "ccnd1", "cdk7", "mad2l2", "ccne1", "ccne2",
  "cdc6", "top2a", "birc5a", "mad2l1", "mad2l1bp",
  "ccna2", "ccnb1", "ccnb3", "cdk1", "cdc20",
  "aurka", "aurkb", "ccnb2", "mki67", "cdk6"
), group.by="CellCycle", additional.group.by = "Phase", draw.lines = TRUE, lines.width = 6, size =3, group.bar.height = 0.05, angle = 45) + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 



DoHeatmap(merge_allLEC_allVEC_7, features = c("lyve1b", "stab2", "cdh6", "tbx1", "cdh5", "kdrl", "pcna", "mcm6", "top2a", "mki67", "cdipt", "lbr", "strip1", "syvn1", "tor1", "b3gnt7", "aph1b", "znf521", "fkbp1ab", "psma6b"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DotPlot(merge_allLEC_allVEC_7, features = c("psma6b", "fkbp1ab", "znf521", "aph1b", "b3gnt7", "tor1", "syvn1", "strip1", "lbr", "cdipt", "mki67", "top2a", "mcm6", "pcna"), dot.scale = 5, scale=TRUE, group.by=c("CellCycle")) + coord_flip() + theme(text = element_text(face = "oblique", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) + scale_color_gradientn(colours = viridis::viridis(20))

DoHeatmap(merge_allLEC_allVEC_48, features = c("lyve1b", "stab2", "cdh6", "tbx1", "cdh5", "kdrl", "pcna", "mcm6", "top2a", "mki67", "lbr", "cse1l", "pom121", "nup58", "ubr7", "asf1ba", "aaas", "suv39h1b", "armc1l", "arhgef39"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 
DotPlot(merge_allLEC_allVEC_48, features = c("arhgef39", "armc1l", "suv39h1b", "aaas", "asf1ba", "ubr7", "nup58", "pom121", "cse1l", "lbr", "mki67", "top2a", "mcm6", "pcna"), dot.scale = 5, scale=TRUE, group.by=c("CellCycle")) + coord_flip() + theme(text = element_text(face = "oblique", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) + scale_color_gradientn(colours = viridis::viridis(20))


#Combining cluster then Phases extra G1, G2M, S 

#48 hpf:

allLECcluster <- subset(adata_48h, idents = c("0", "2", "4"))
allVECcluster <- subset(adata_48h, idents = c("1", "3"))

allLECcluster@active.ident=as.factor(allLECcluster$Phase)
allVECcluster@active.ident=as.factor(allVECcluster$Phase)

allLECS48 <- subset(allLECcluster, ident = "S")
allLECG2M48 <- subset(allLECcluster, ident = "G2M")
allLECG148 <- subset(allLECcluster, ident = "G1")
allVECS48 <- subset(allVECcluster, ident = "S")
allVECG2M48 <- subset(allVECcluster, ident = "G2M")
allVECG148 <- subset(allVECcluster, ident = "G1")

ccSL = rep("48h-S-LEC", ncol(allLECS48))
ccG2ML = rep("48h-G2M-LEC", ncol(allLECG2M48))
ccG1L = rep("48h-G1-LEC", ncol(allLECG148))
ccSV = rep("48h-S-VEC", ncol(allVECS48))
ccG2MV = rep("48h-G2M-VEC", ncol(allVECG2M48))
ccG1V = rep("48h-G1-VEC", ncol(allVECG148))

allLECS48 = AddMetaData(allLECS48, factor(ccSL), col.name = "Phase")
allLECG2M48 = AddMetaData(allLECG2M48, factor(ccG2ML), col.name = "Phase")
allLECG148 = AddMetaData(allLECG148, factor(ccG1L), col.name = "Phase")
allVECS48 = AddMetaData(allVECS48, factor(ccSV), col.name = "Phase")
allVECG2M48 = AddMetaData(allVECG2M48, factor(ccG2MV), col.name = "Phase")
allVECG148 = AddMetaData(allVECG148, factor(ccG1V), col.name = "Phase")

merge_allLECS48_allLECG2M48 = merge(allLECS48, allLECG2M48)
merge_allLECS48_allLECG2M48_allLECG148 = merge(merge_allLECS48_allLECG2M48, allLECG148)
merge_allVECS48_allVECG2M48 = merge(allVECS48, allVECG2M48)
merge_allVECS48_allVECG2M48_allVECG148 = merge(merge_allVECS48_allVECG2M48, allVECG148)
merge_allLEC_allVEC_48_phase = merge(merge_allLECS48_allLECG2M48_allLECG148, merge_allVECS48_allVECG2M48_allVECG148)


merge_allLEC_allVEC_48_phase@active.ident=as.factor(merge_allLEC_allVEC_48_phase$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_48_phase, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_48_phase@active.ident=as.factor(merge_allLEC_allVEC_48_phase$Phase)
DoHeatmap(merge_allLEC_allVEC_48_phase, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

#7 dpf:

allLECcluster7 <- subset(adata_7d, idents = c("1", "4"))
allVECcluster7 <- subset(adata_7d, idents = c("0", "2", "5", "6"))

allLECcluster7@active.ident=as.factor(allLECcluster7$Phase)
allVECcluster7@active.ident=as.factor(allVECcluster7$Phase)

allLECS7 <- subset(allLECcluster7, ident = "S")
allLECG2M7 <- subset(allLECcluster7, ident = "G2M")
allLECG17 <- subset(allLECcluster7, ident = "G1")
allVECS7 <- subset(allVECcluster7, ident = "S")
allVECG2M7 <- subset(allVECcluster7, ident = "G2M")
allVECG17 <- subset(allVECcluster7, ident = "G1")

ccSL7 = rep("7d-S-LEC", ncol(allLECS7))
ccG2ML7 = rep("7d-G2M-LEC", ncol(allLECG2M7))
ccG1L7 = rep("7d-G1-LEC", ncol(allLECG17))
ccSV7 = rep("7d-S-VEC", ncol(allVECS7))
ccG2MV7 = rep("7d-G2M-VEC", ncol(allVECG2M7))
ccG1V7 = rep("7d-G1-VEC", ncol(allVECG17))

allLECS7 = AddMetaData(allLECS7, factor(ccSL7), col.name = "Phase")
allLECG2M7 = AddMetaData(allLECG2M7, factor(ccG2ML7), col.name = "Phase")
allLECG17 = AddMetaData(allLECG17, factor(ccG1L7), col.name = "Phase")
allVECS7 = AddMetaData(allVECS7, factor(ccSV7), col.name = "Phase")
allVECG2M7 = AddMetaData(allVECG2M7, factor(ccG2MV7), col.name = "Phase")
allVECG17 = AddMetaData(allVECG17, factor(ccG1V7), col.name = "Phase")

merge_allLECS7_allLECG2M7 = merge(allLECS7, allLECG2M7)
merge_allLECS7_allLECG2M7_allLECG17 = merge(merge_allLECS7_allLECG2M7, allLECG17)
merge_allVECS7_allVECG2M7 = merge(allVECS7, allVECG2M7)
merge_allVECS7_allVECG2M7_allVECG17 = merge(merge_allVECS7_allVECG2M7, allVECG17)
merge_allLEC_allVEC_7_phase = merge(merge_allLECS7_allLECG2M7_allLECG17, merge_allVECS7_allVECG2M7_allVECG17)

merge_allLEC_allVEC_7_phase@active.ident=as.factor(merge_allLEC_allVEC_7_phase$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_7_phase, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_7_phase@active.ident=as.factor(merge_allLEC_allVEC_7_phase$Phase)
DoHeatmap(merge_allLEC_allVEC_7_phase, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 



merge_allLEC_allVEC_7_phase@active.ident=as.factor(merge_allLEC_allVEC_7_phase$Phase)
DotPlot(merge_allLEC_allVEC_7_phase, features = c("psma6b", "fkbp1ab", "znf521", "aph1b", "b3gnt7", "tor1", "syvn1", "strip1", "lbr", "cdipt"), dot.scale = 5, scale=TRUE, group.by=c("Phase")) + coord_flip() + theme(text = element_text(face = "oblique", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) + scale_color_gradientn(colours = viridis::viridis(20))

merge_allLEC_allVEC_48_phase@active.ident=as.factor(merge_allLEC_allVEC_48_phase$Phase)
DotPlot(merge_allLEC_allVEC_48_phase, features = c("arhgef39", "armc1l", "suv39h1b", "aaas", "asf1ba", "ubr7", "nup58", "pom121", "cse1l", "lbr"), dot.scale = 5, scale=TRUE, group.by=c("Phase")) + coord_flip() + theme(text = element_text(face = "oblique", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) + scale_color_gradientn(colours = viridis::viridis(20))


#Combining cluster then Prol vs G1 all time points

#48 hpf:

allLECcluster <- subset(adata_48h, idents = c("0", "2", "4"))
allVECcluster <- subset(adata_48h, idents = c("1", "3"))

allLECcluster@active.ident=as.factor(allLECcluster$CellCycle)
allVECcluster@active.ident=as.factor(allVECcluster$CellCycle)
allLECprol48 <- subset(allLECcluster, ident = "prol")
allLECG148 <- subset(allLECcluster, ident = "G1")
allVECprol48 <- subset(allVECcluster, ident = "prol")
allVECG148 <- subset(allVECcluster, ident = "G1")

ccprolL = rep("48h-prol-LEC", ncol(allLECprol48))
ccG1L = rep("48h-G1-LEC", ncol(allLECG148))
ccprolV = rep("48h-prol-VEC", ncol(allVECprol48))
ccG1V = rep("48h-G1-VEC", ncol(allVECG148))

allLECprol48 = AddMetaData(allLECprol48, factor(ccprolL), col.name = "CellCycle")
allLECG148 = AddMetaData(allLECG148, factor(ccG1L), col.name = "CellCycle")
allVECprol48 = AddMetaData(allVECprol48, factor(ccprolV), col.name = "CellCycle")
allVECG148 = AddMetaData(allVECG148, factor(ccG1V), col.name = "CellCycle")

merge_allLECprol48_allLECG148 = merge(allLECprol48, allLECG148)
merge_allVECprol48_allVECG148 = merge(allVECprol48, allVECG148)
merge_allLEC_allVEC_48 = merge(merge_allLECprol48_allLECG148, merge_allVECprol48_allVECG148)


merge_allLEC_allVEC_48@active.ident=as.factor(merge_allLEC_allVEC_48$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_48, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_48@active.ident=as.factor(merge_allLEC_allVEC_48$"CellCycle")
DoHeatmap(merge_allLEC_allVEC_48, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 




#7 dpf:

allLECcluster7 <- subset(adata_7d, idents = c("1", "4"))
allVECcluster7 <- subset(adata_7d, idents = c("0", "2", "5", "6"))

allLECcluster7@active.ident=as.factor(allLECcluster7$CellCycle)
allVECcluster7@active.ident=as.factor(allVECcluster7$CellCycle)
allLECprol7 <- subset(allLECcluster7, ident = "prol")
allLECG17 <- subset(allLECcluster7, ident = "G1")
allVECprol7 <- subset(allVECcluster7, ident = "prol")
allVECG17 <- subset(allVECcluster7, ident = "G1")

ccprolL7 = rep("7d-prol-LEC", ncol(allLECprol7))
ccG1L7 = rep("7d-G1-LEC", ncol(allLECG17))
ccprolV7 = rep("7d-prol-VEC", ncol(allVECprol7))
ccG1V7 = rep("7d-G1-VEC", ncol(allVECG17))

allLECprol7 = AddMetaData(allLECprol7, factor(ccprolL7), col.name = "CellCycle")
allLECG17 = AddMetaData(allLECG17, factor(ccG1L7), col.name = "CellCycle")
allVECprol7 = AddMetaData(allVECprol7, factor(ccprolV7), col.name = "CellCycle")
allVECG17 = AddMetaData(allVECG17, factor(ccG1V7), col.name = "CellCycle")

merge_allLECprol7_allLECG17 = merge(allLECprol7, allLECG17)
merge_allVECprol7_allVECG17 = merge(allVECprol7, allVECG17)
merge_allLEC_allVEC_7 = merge(merge_allLECprol7_allLECG17, merge_allVECprol7_allVECG17)

merge_allLEC_allVEC_7@active.ident=as.factor(merge_allLEC_allVEC_7$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_7, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_7@active.ident=as.factor(merge_allLEC_allVEC_7$"CellCycle")
DoHeatmap(merge_allLEC_allVEC_7, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 




#36 hpf:

allLECcluster36 <- subset(adata_36h, idents = c("0", "2", "4"))
allVECcluster36 <- subset(adata_36h, idents = c("1", "5"))

allLECcluster36@active.ident=as.factor(allLECcluster36$CellCycle)
allVECcluster36@active.ident=as.factor(allVECcluster36$CellCycle)
allLECprol36 <- subset(allLECcluster36, ident = "prol")
allLECG136 <- subset(allLECcluster36, ident = "G1")
allVECprol36 <- subset(allVECcluster36, ident = "prol")
allVECG136 <- subset(allVECcluster36, ident = "G1")

ccprolL = rep("36h-prol-LEC", ncol(allLECprol36))
ccG1L = rep("36h-G1-LEC", ncol(allLECG136))
ccprolV = rep("36h-prol-VEC", ncol(allVECprol36))
ccG1V = rep("36h-G1-VEC", ncol(allVECG136))

allLECprol36 = AddMetaData(allLECprol36, factor(ccprolL), col.name = "CellCycle")
allLECG136 = AddMetaData(allLECG136, factor(ccG1L), col.name = "CellCycle")
allVECprol36 = AddMetaData(allVECprol36, factor(ccprolV), col.name = "CellCycle")
allVECG136 = AddMetaData(allVECG136, factor(ccG1V), col.name = "CellCycle")

merge_allLECprol36_allLECG136 = merge(allLECprol36, allLECG136)
merge_allVECprol36_allVECG136 = merge(allVECprol36, allVECG136)
merge_allLEC_allVEC_36 = merge(merge_allLECprol36_allLECG136, merge_allVECprol36_allVECG136)


merge_allLEC_allVEC_36@active.ident=as.factor(merge_allLEC_allVEC_36$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_36, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_36@active.ident=as.factor(merge_allLEC_allVEC_36$"CellCycle")
DoHeatmap(merge_allLEC_allVEC_36, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 




#5 dpf:

allLECcluster5 <- subset(adata_5d, idents = c("0", "7", "8"))
allVECcluster5 <- subset(adata_5d, idents = c("1", "2", "3", "4", "5"))

allLECcluster5@active.ident=as.factor(allLECcluster5$CellCycle)
allVECcluster5@active.ident=as.factor(allVECcluster5$CellCycle)
allLECprol5 <- subset(allLECcluster5, ident = "prol")
allLECG15 <- subset(allLECcluster5, ident = "G1")
allVECprol5 <- subset(allVECcluster5, ident = "prol")
allVECG15 <- subset(allVECcluster5, ident = "G1")

ccprolL = rep("5d-prol-LEC", ncol(allLECprol5))
ccG1L = rep("5d-G1-LEC", ncol(allLECG15))
ccprolV = rep("5d-prol-VEC", ncol(allVECprol5))
ccG1V = rep("5d-G1-VEC", ncol(allVECG15))

allLECprol5 = AddMetaData(allLECprol5, factor(ccprolL), col.name = "CellCycle")
allLECG15 = AddMetaData(allLECG15, factor(ccG1L), col.name = "CellCycle")
allVECprol5 = AddMetaData(allVECprol5, factor(ccprolV), col.name = "CellCycle")
allVECG15 = AddMetaData(allVECG15, factor(ccG1V), col.name = "CellCycle")

merge_allLECprol5_allLECG15 = merge(allLECprol5, allLECG15)
merge_allVECprol5_allVECG15 = merge(allVECprol5, allVECG15)
merge_allLEC_allVEC_5 = merge(merge_allLECprol5_allLECG15, merge_allVECprol5_allVECG15)


merge_allLEC_allVEC_5@active.ident=as.factor(merge_allLEC_allVEC_5$seurat_clusters)
DoHeatmap(merge_allLEC_allVEC_5, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

merge_allLEC_allVEC_5@active.ident=as.factor(merge_allLEC_allVEC_5$"CellCycle")
DoHeatmap(merge_allLEC_allVEC_5, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 




#combining all LEC over time early and late:

merge_allLECprol_allLECG1_earlytime = merge(merge_allLECprol36_allLECG136, merge_allLECprol48_allLECG148)

merge_allLECprol_allLECG1_earlytime@active.ident=as.factor(merge_allLECprol_allLECG1_earlytime$"CellCycle")
DoHeatmap(merge_allLECprol_allLECG1_earlytime, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

DotPlot(merge_allLECprol_allLECG1_earlytime, features = c("arhgef39", "armc1l", "suv39h1b", "aaas", "asf1ba", "ubr7", "nup58", "pom121", "cse1l", "lbr", "mki67", "top2a", "mcm6", "pcna"), dot.scale = 5, scale=TRUE, group.by=c("CellCycle")) + coord_flip() + theme(text = element_text(face = "oblique", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) + scale_color_gradientn(colours = viridis::viridis(20))




merge_allLECprol_allLECG1_latetime = merge(merge_allLECprol5_allLECG15, merge_allLECprol7_allLECG17)

merge_allLECprol_allLECG1_latetime@active.ident=as.factor(merge_allLECprol_allLECG1_latetime$"CellCycle")
DoHeatmap(merge_allLECprol_allLECG1_latetime, features = c("fli1a", "lyve1b", "prox1a", "tbx1", "mrc1a", "mafba", "mafbb", "cldn11b", "cldn11a", "gpr183a", "cdh6", "hapln3", "flt4", "cdh5", "kdrl", "pcna", "cdk2", "cdk4", "cdk6", "cdk7", "mad2l2", "mad2l1", "mcm2", "mcm3", "mcm4", "mcm5", "mcm6", "mcm7", "ccna2", "ccna1", "ccnb1", "ccnb2", "ccnd1", "ccnd3", "ccne1", "ccne2", "cdk1", "cdc6", "cdc20", "aurkb"), draw.lines = TRUE, lines.width = 6, slot = "data") + (scale_fill_gradientn(colors = viridis::viridis(20), na.value = "white")) & theme(text = element_text(face = "oblique", size = 10),  axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) 

DotPlot(merge_allLECprol_allLECG1_latetime, features = c("psma6b", "fkbp1ab", "znf521", "aph1b", "b3gnt7", "tor1", "syvn1", "strip1", "lbr", "cdipt", "mki67", "top2a", "mcm6", "pcna"), dot.scale = 5, scale=TRUE, group.by=c("CellCycle")) + coord_flip() + theme(text = element_text(face = "oblique", size = 10), axis.text.x=element_text(size=10), axis.title = element_text(size=10), axis.title.y.right = element_text(size=10), legend.text=element_text(size=10), legend.title=element_text(size=10), axis.line = element_line(size=1)) + scale_color_gradientn(colours = viridis::viridis(20))




