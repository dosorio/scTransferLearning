library(Seurat)
library(ggplot2)
library(Matrix)
library(harmony)
library(symphony)
library(caret)
library(pbapply)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(ggforce)
library(patchwork)
library(reshape2)
library(statsExpressions)
library(ggdendro)
library(pvclust)

# quantile_normalisation <- function(df){
#   df_rank <- apply(df,2,rank,ties.method="min")
#   df_sorted <- data.frame(apply(df, 2, sort))
#   df_mean <- apply(df_sorted, 1, mean)
#
#   index_to_mean <- function(my_index, my_mean){
#     return(my_mean[my_index])
#   }
#
#   df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
#   rownames(df_final) <- rownames(df)
#   return(df_final)
# }
source('https://raw.githubusercontent.com/dosorio/utilities/master/data.frame2matrix.R')

ENSEMBL <- read.table('../Data/hsaENSEMBL-GENES.txt', sep = '\t', header = TRUE)
ENSEMBL <- ENSEMBL[ENSEMBL$Gene.name != '',]
GENEID <- ENSEMBL$Gene.name
names(GENEID) <- ENSEMBL$Gene.stable.ID

getCellLines <- function(X){
  cellLines <- colnames(X)
  cellLines <- unlist(lapply(strsplit(cellLines, '_'), function(X){X[1]}))
  return(cellLines)
}


X <- readRDS('../Data/RAW.UMI.counts.BC.cell.lines.rds')
X <- X[rownames(X) %in% ENSEMBL$Gene.stable.ID,]
rownames(X) <- as.vector(GENEID[rownames(X)])
X <- CreateSeuratObject(X)

CL <- pbsapply(unique(Idents(X)), function(CT){
  rowSums(X@assays$RNA@counts[,Idents(X) %in% CT])
})
colnames(CL) <- unique(Idents(X))

clBRCA <- buildReference(
  exp_ref = X@assays$RNA@counts,
  metadata_ref = data.frame(g = 1, cellLine = getCellLines(X@assays$RNA@counts)),
  vars = 'g',
  do_umap = TRUE,
  verbose = TRUE,
  d = 50, save_uwot_path = 'umapBRCA')


query <- readMM('../Data/count_matrix_sparse.mtx')
rownames(query) <- readLines('../Data/count_matrix_genes.tsv')
colnames(query) <- readLines('../Data/count_matrix_barcodes.tsv')
queryMetadata <- read.csv('../Data/metadata.csv', row.names = 1)

DC <- pbsapply(unique(queryMetadata$orig.ident), function(Donor){
  rowSums(query[,grepl(Donor, colnames(query))])
})
write.csv(DC, '../Results/allCellsPatientsPseudobulkProfiles.csv')

queryMetadata <- queryMetadata[grepl('Cancer',queryMetadata$celltype_minor),]
query <- query[,rownames(queryMetadata)]

DCC <- pbsapply(unique(queryMetadata$orig.ident), function(Donor){
  rowSums(query[,grepl(Donor, colnames(query))])
})
write.csv(DCC, '../Results/cancerCellsPatientPseudobulkProfiles.csv')

donorSubType <- queryMetadata$subtype
names(donorSubType) <- queryMetadata$orig.ident

donor <- 'CID44971'
donorData <- query[,(queryMetadata$orig.ident %in% donor)]
cNames <- colnames(donorData)

LOOCV <- pbsapply(colnames(donorData), function(C){
  donorData <- donorData[,!cNames %in% C]
  donorMD <- data.frame(donor = rep(donor, ncol(donorData)))
  qmap <- mapQuery(donorData, metadata_query = donorMD, ref_obj = clBRCA, do_umap = TRUE, do_normalize = TRUE)
  qmap <- knnPredict(qmap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  table(qmap$meta_data$cell_type_pred_knn)
})
write.csv(LOOCV, '../Results/ccLOOCV.csv')

geneList <- intersect(rownames(CL), rownames(DC))
COMBN <- data.frame(CL[geneList,], DC[geneList,])
COMBN <- as.matrix(COMBN)
COMBN <- log1p((t(t(COMBN)/colSums(COMBN)))*1e4)
spCor <- function(x){as.dist(cor(x, method = 'sp'))}
O <- pvclust(COMBN, method.dist = spCor, parallel = TRUE, nboot = 1000, method.hclust = 'complete')
png('../Figures/S3.png', width = 4000, height = 1000, res = 300)
par(mar=c(1,4,1,1))
plot(O, print.pv = 'bp', print.num = FALSE, main = '', sub = '', xlab = '')
dev.off()

geneList <- intersect(rownames(CL), rownames(DCC))
COMBN <- data.frame(CL[geneList,], DCC[geneList,])
COMBN <- as.matrix(COMBN)
COMBN <- log1p((t(t(COMBN)/colSums(COMBN)))*1e4)
spCor <- function(x){as.dist(cor(x, method = 'sp'))}
O <- pvclust(COMBN, method.dist = spCor, parallel = TRUE, nboot = 1000, method.hclust = 'complete')
png('../Figures/S4.png', width = 4000, height = 1000, res = 300)
par(mar=c(1,4,1,1))
plot(O, print.pv = 'bp', print.num = FALSE, main = '', sub = '', xlab = '')
dev.off()
