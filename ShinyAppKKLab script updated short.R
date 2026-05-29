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

#datadir needs so be adjusted to your file location:

datadir = "/Users/margl916/Downloads/v3_SCT/"
adata_5d <- readRDS(file.path(datadir,"data_5days.Rds"))
adata_7d <- readRDS(file.path(datadir,"data_7days.Rds"))
adata_36h <- readRDS(file.path(datadir,"data_36h.Rds"))
adata_48h <- readRDS(file.path(datadir,"data_48h.Rds"))

#move/copy shinyAppMulti to Home directory

shiny::runApp('shinyAppMulti')



#next time, run only this:

shiny::runApp('shinyAppMulti')


