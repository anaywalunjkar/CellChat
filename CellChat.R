library(CellChat)
library(Seurat)
library(reticulate)
setwd ("C:/RStudiowork/prac/cellchat")
cellchat <-readRDS("ependymal_cells.rds")
library(SeuratObject)
cts [1:10,1:10]
data.input <- GetAssayData(cellchat, assay = "RNA", slot = "data") 
class
data.input <- GetAssayData(cellchat, assay = "RNA", slot = "data")
cellchat <- UpdateSeuratObject(cellchat)
Assays(cellchat)
data.input <- GetAssayData(cellchat, assay = "RNA", slot = "data") 
labels <- Idents(cellchat)
meta <- data.frame(group = labels, row.names = names(labels)) 
cellchatobj <- createCellChat(object = data.input, meta = meta, group.by = "group")

#databases
CellChatDB.human <- CellChatDB.human
CellChatDB.mouse <- CellChatDB.mouse

showDatabaseCategory(CellChatDB.human)
showDatabaseCategory(CellChatDB.mouse)

CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
interaction <- CellChatDB.human$interaction
write.csv(interaction,  file = "interaction.csv")
CellChatDB.use <- CellChatDB.human
cellchatobj@DB <- CellChatDB.use
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchatobj <- subsetData(cellchat)
cellchatobj <- identifyOverExpressedGenes(cellchatobj)
cellchatobj <- subsetData(cellchatobj)
str(cellchatobj)
cellchatobj <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchatobj <- subsetData(cellchatobj)
str(cellchatobj)
CellChatDB.use <- CellChatDB.human
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchatobj@DB <- CellChatDB.use
cellchatobj <- subsetData(cellchatobj)
cellchatobj <- identifyOverExpressedGenes(cellchatobj)
class(cellchatobj@data.signaling)
if (is.numeric(cellchatobj@data.signaling)) {
  cellchatobj@data.signaling <- matrix(cellchatobj@data.signaling, nrow = 3)
}
validObject(cellchatobj)
cellchatobj <- subsetData(cellchat)
cellchatobj <- subsetData(cellchatobj)
class(cellchatobj@data.signaling)
str(cellchatobj@data.signaling)

if (is.numeric(cellchatobj@data.signaling)) {
  cellchatobj@data.signaling <- matrix(cellchatobj@data.signaling, n_columns = 3)
}
if (is.numeric(cellchatobj@data.signaling)) {
  cellchatobj@data.signaling <- matrix(cellchatobj@data.signaling, ncol = 3)
}
if (is.numeric(cellchatobj@data.signaling)) {
  
  n_columns <- 3  
  n_rows <- length(cellchatobj@data.signaling) / n_columns
  
  if (n_rows %% 1 != 0) {
    stop("The number of elements is not divisible by the specified number of columns.")
  }
  
  cellchatobj@data.signaling <- matrix(cellchatobj@data.signaling, ncol = n_columns)
}
# Example numeric vector for @data.signaling
cellchatobj@data.signaling <- 1:9  # Dummy data with 9 elements

# Reshape the vector into a matrix with 3 rows
cellchatobj@data.signaling <- matrix(cellchatobj@data.signaling, nrow = 3)

# Print the resulting matrix to verify
print(cellchatobj@data.signaling)

cellchatobj@data.signaling <- matrix(cellchatobj@data.signaling, ncol = 3, nrow = 3)

# Print the resulting matrix
print(cellchatobj@data.signaling)
cellchatobj <- subsetData(cellchatobj)
cellchatobj@data.signaling <- 1:9
CellChatDB.human <- CellChatDB.human
CellChatDB.use <- CellChatDB.human
cellchatobj@DB <- CellChatDB.use
cellchatobj <- subsetData(cellchatobj)
