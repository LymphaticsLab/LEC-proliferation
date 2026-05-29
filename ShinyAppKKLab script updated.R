reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "ggplot2", "gridExtra", "glue", "readr", "RColorBrewer", "R.utils", "Seurat")

newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "hdf5r", "reticulate", "ggplot2", "gridExtra", "magrittr", "ggdendro")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

devtools::install_github("SGDDNB/ShinyCell")

knitr::opts_chunk$set(message=FALSE, warning=FALSE, result='hold',fig.width=10, fig.height = 8)
library(SeuratObject)
library(Seurat)
library(ShinyCell)
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

shiny::runApp('shinyAppMulti')








#this really needed??? Test, if i send the scConf or the other scripts files and folder files of shinyappmulti

scConf1 = createConfig(adata_36h)
scConf2 = createConfig(adata_48h)
scConf3 = createConfig(adata_5d)
scConf4 = createConfig(adata_7d)

makeShinyFiles(adata_36h, scConf1, gex.assay = "SCT", gex.slot = "data",
               gene.mapping = FALSE, shiny.prefix = "sc1",
               shiny.dir = "shinyAppMulti/",
               default.gene1 = "prox1a", default.gene2 = "cdh6",
               default.multigene = c("fli1a", "lyve1b", "stab2", "cldn11b", "cdh6", "prox1a", "cdh5", "kdrl"),
               default.dimred = c("UMAP_1", "UMAP_2"))

makeShinyFiles(adata_48h, scConf2, gex.assay = "SCT", gex.slot = "data",
               gene.mapping = FALSE, shiny.prefix = "sc2",
               shiny.dir = "shinyAppMulti/",
               default.gene1 = "prox1a", default.gene2 = "cdh6",
               default.multigene = c("fli1a", "lyve1b", "stab2", "cldn11b", "cdh6", "prox1a", "cdh5", "kdrl"),
               default.dimred = c("UMAP_1", "UMAP_2"))

makeShinyFiles(adata_5d, scConf3, gex.assay = "SCT", gex.slot = "data",
               gene.mapping = FALSE, shiny.prefix = "sc3",
               shiny.dir = "shinyAppMulti/",
               default.gene1 = "prox1a", default.gene2 = "cdh6",
               default.multigene = c("fli1a", "lyve1b", "stab2", "cldn11b", "cdh6", "prox1a", "cdh5", "kdrl"),
               default.dimred = c("UMAP_1", "UMAP_2"))

makeShinyFiles(adata_7d, scConf4, gex.assay = "SCT", gex.slot = "data",
               gene.mapping = FALSE, shiny.prefix = "sc4",
               shiny.dir = "shinyAppMulti/",
               default.gene1 = "prox1a", default.gene2 = "cdh6",
               default.multigene = c("fli1a", "lyve1b", "stab2", "cldn11b", "cdh6", "prox1a", "cdh5", "kdrl"),
               default.dimred = c("UMAP_1", "UMAP_2"))

makeShinyCodesMulti(
  shiny.title = "ShinyCell KKlab ZF scRNAseq dataset",  shiny.footnotes = "Test MG citation xx",
  shiny.prefix = c("sc1", "sc2", "sc3", "sc4"),
  shiny.headers = c("36 hpf", "48 hpf", "5 dpf", "7 dpf"), 
  shiny.dir = "shinyAppMulti/")


scConf1 = delMeta(scConf1, c("orig.ident", "PU", "plate_row", "Timepoint", "Transgenic.line", "Date", "Flowcell", "Prox1a", "RNA_snn_res.0.5", "integrated_snn_res.0.5", "integrated_snn_res.2", "Flt4"))
scConf1 = scConf1[-3,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-5,]
scConf1 = scConf1[-13,]

showLegend(scConf1)

scConf1 = reorderMeta(scConf1, scConf1$ID[c(13,8,12,11,9,10,14,15,1:7)])
showOrder(scConf1)
scConf1 = modDefault(scConf1, "seurat_clusters", "Part")
showOrder(scConf1)

# XXXX make all others

scConf1 = modColours(scConf1, meta.to.mod = "seurat_clusters", 
                     new.colours= c("0" = "#109630","1"="#762a83", "2" = "#a4da84", "3" = "#fcf2a7", "4" = "#f49262", "5" = "#c5a9ce"))
scConf1 = modLabels(scConf1, meta.to.mod = "seurat_clusters", 
                    new.labels = c("LEC", "VEC", "2nd sprouts", "heart EC", "prol LEC", "trunk VEC"))
showLegend(scConf1)


scConf2 = modColours(scConf2, meta.to.mod = "seurat_clusters", 
                     new.colours= c("0" = "#a4da84","1"="#762a83", "2" = "#109630", "3" = "#96c7dd", "4" = "#f49262"))
scConf2 = modLabels(scConf2, meta.to.mod = "seurat_clusters", 
                    new.labels = c("maturing 2nd sprouts", "VEC", "LEC", "facial VEC", "prol LEC"))
showLegend(scConf2)



scConf3 = modColours(scConf3, meta.to.mod = "seurat_clusters", 
                     new.colours= c("0" = "#109630", "1"="#762a83", "2" = "#96c7dd", "3" = "#54859b", "4" = "#344a82", "5" = "#c5a9ce", "6" = "#4c4c4c", "7" = "#f49262", "8" = "#a4da84"))
scConf3 = modLabels(scConf3, meta.to.mod = "seurat_clusters", 
                    new.labels = c("LEC", "VEC", "facial VEC", "metabolic VEC", "cerebral VEC", "sprouting VEC", "non-EC", "OLV LEC", "cerebral LEC"))
showLegend(scConf3)



scConf4 = modColours(scConf4, meta.to.mod = "seurat_clusters", 
                     new.colours= c("0" = "#762a83","1"="#109630", "2" = "#c5a9ce", "3" = "#4c4c4c", "4" = "#f49262", "5" = "#344a82", "6" = "#96c7dd"))
scConf4 = modLabels(scConf4, meta.to.mod = "seurat_clusters", 
                    new.labels = c("VEC", "LEC", "sprouting VEC", "non-EC", "OLV LEC", "cerebral LEC", "prol VEC"))
showLegend(scConf4)



scConf1 = modColours(scConf1, meta.to.mod = "Phase", new.colours= c("gray", "#5ab4ac", "#01665e"))
scConf1 = modColours(scConf1, meta.to.mod = "CellCycle", new.colours= c("gray", "#01665e"))
scConf1 = modColours(scConf1, meta.to.mod = "Part", new.colours= c("#5ab4ac", "#d8b365", "#737373"))
showLegend(scConf1)

# XXXX make all others



shiny::runApp('shinyAppMulti')



