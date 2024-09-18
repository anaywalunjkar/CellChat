setwd ("C:/RStudiowork/prac/cellchat")
library(CellChat)
library(Seurat)
library(reticulate)
library(SeuratObject)
library(patchwork)
library(circlize)
library(tidyverse)
library(patchwork)
options(stringsAsFactors = FALSE)


Seuratobj <-readRDS("ependymal_cells.rds")
Seuratobj <- UpdateSeuratObject(Seuratobj)
data.input <- GetAssayData(Seuratobj, assay = "RNA", slot = "data") 
labels <- Idents(Seuratobj)
meta <- data.frame(group = labels, row.names = names(labels)) 
#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
#cellChat <- createCellChat(object = Seuratobj, group.by = "ident", assay = "RNA")
cellChat <- createCellChat(object = data.input, meta = meta, group.by = "group")

print(cellChat)

CellChatDB.use <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB.mouse, search = "Secreted Signaling")
cellChat@DB <- CellChatDB.use
cellChat <- subsetData(cellChat)


cellChat <- identifyOverExpressedGenes(cellChat, do.fast = FALSE)
cellChat <- identifyOverExpressedInteractions(cellChat)
cellChat <- computeCommunProb(cellChat, raw.use = FALSE)
dim(cellChat@data)

cellChat <- filterCommunication(cellChat, min.cells = 10)
cellChat <- computeCommunProbPathway(cellChat)
