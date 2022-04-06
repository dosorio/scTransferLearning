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
library(reshape2)
library(statsExpressions)
library(BisqueRNA)
library(Biobase)

source('https://raw.githubusercontent.com/dosorio/utilities/master/data.frame2matrix.R')
source('https://raw.githubusercontent.com/dosorio/utilities/master/ggColors.R')

# Loading dictionary of gene symbols
ENSEMBL <- read.table('../Data/hsaENSEMBL-GENES.txt', sep = '\t', header = TRUE)
ENSEMBL <- ENSEMBL[ENSEMBL$Gene.name != '',]
GENEID <- ENSEMBL$Gene.name
names(GENEID) <- ENSEMBL$Gene.stable.ID

# Function to extract the cell lines ID from the data
getCellLines <- function(X){
  cellLines <- colnames(X)
  cellLines <- unlist(lapply(strsplit(cellLines, '_'), function(X){X[1]}))
  return(cellLines)
}

# Loading breast cancer cell atlas
X <- readRDS('../Data/RAW.UMI.counts.BC.cell.lines.rds')
X <- X[rownames(X) %in% ENSEMBL$Gene.stable.ID,]
rownames(X) <- as.vector(GENEID[rownames(X)])
X <- CreateSeuratObject(X)
X <- NormalizeData(X)

clMD <- data.frame(CL = unique(getCellLines(X)), Type = NA)
clMD$Type <- c('Her2+', 'TNBC-A', 'TNBC-A', 'TNBC-A', 'Lum',
               'TNBC-A', 'TNBC-A', 'Her2+', 'Lum', 'TNBC-A',
               'TNBC-B', 'Her2+', 'TNBC-B', 'Lum', 'TNBC-B',
               'Her2+', 'TNBC-B', 'Lum', 'Lum', 'TNBC-A',
               'Lum', 'Her2+', 'TNBC-A', 'Her2+', 'TNBC-A',
               'TNBC-A', 'NI', 'TNBC-A', 'Lum', 'Lum',
               'TNBC-B', 'Her2+')
clMD <- clMD[order(clMD$CL),]
cellSubtype <- clMD$Type
names(cellSubtype) <- clMD$CL
clColor <- ggColors(nrow(clMD))
names(clColor) <- clMD$CL

# Pseudobulk analysis
PB <- sapply(levels(Idents(X)), function(I){
  rowSums(X@assays$RNA@data[,Idents(X) %in% I])
})
write.table(PB, quote = FALSE, sep = '\t', file = '../Results/referenceCSx.tsv')
set.seed(0)
cibersortxTest <- sample(colnames(X), 1000)
cibersortxTest <- X@assays$RNA@data[,cibersortxTest]
csTrue <- table(getCellLines(cibersortxTest))/1000
write.table(csTrue, '../Results/trueProportionsCSx.tsv',quote = FALSE, sep = '\t')
cibersortxTest <- data.frame(T1=rowSums(cibersortxTest), T2=rowSums(cibersortxTest))
cibersortxTest <- cibersortxTest[rowSums(is.finite(as.matrix(cibersortxTest))) == 2,]
write.table(cibersortxTest, quote = FALSE, sep = '\t', file = '../Results/samplesCSx.tsv')

cibersortXoutput <- read.csv('../Results/CIBERSORTx_Job12_Results.csv')
cibersortXoutput <- cibersortXoutput[1,names(csTrue)]
PO <- t(rbind(as.vector(csTrue), c(cibersortXoutput)))
PO <- as.data.frame(PO)
PO <- data.frame(Expected=unlist(PO$V1), Predicted= unlist(PO$V2))
PO <- round(PO*100,2)
PO$CT <- rownames(PO)

png('../Figures/S2.png', width = 1800, height = 1800, res = 300)
ggplot(PO, aes(Expected, Predicted, label = CT)) +
  geom_smooth(method = 'lm', color = 'red') +
  geom_point() +
  theme_bw() +
  labs(title = 'CIBERSORTx') +
  theme(plot.title = element_text(face = 2)) +
  xlab('Expected Proportion (%)') +
  ylab('Predicted Proportion (%)') +
  geom_text_repel(min.segment.length = 0, segment.size =0.1) +
  labs(subtitle = statsExpressions::corr_test(PO, Expected, Predicted, type = 'nonp')$expression[[1]])
dev.off()

# Cross validation
steps <- seq(from = (ncol(X)-1000), to = 1000, by = -1000)
CM <- pbsapply(steps, function(nCells){
  set.seed(0)
  testData <- sample(seq_len(ncol(X)), nCells)
  trainData <- seq_len(ncol(X))[!seq_len(ncol(X)) %in% testData]

  clBRCA <- buildReference(
    exp_ref = X@assays$RNA@counts[,trainData],
    metadata_ref = data.frame(g = 1, cellLine = getCellLines(X@assays$RNA@counts[,trainData])),
    vars = 'g',
    do_umap = TRUE,
    verbose = TRUE,
    d = 50)

  testMap <- mapQuery(X@assays$RNA@counts[,testData],
                      metadata_query = data.frame(cl = getCellLines(X@assays$RNA@counts[,testData])),
                      ref_obj = clBRCA)
  testMap <- knnPredict(testMap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  CM <- confusionMatrix(data = as.factor(testMap$meta_data[,1]), reference = as.factor(testMap$meta_data[,2]), mode = 'everything')
  CM$overall
})

ACC <- data.frame(t(CM), steps)
ACC$steps <- (ncol(X)-ACC$steps)
ACC$Accuracy <- ACC$Accuracy * 100
ACC$AccuracyLower <- ACC$AccuracyLower * 100
ACC$AccuracyUpper <- ACC$AccuracyUpper * 100
write.csv(ACC, '../Results/F2B.csv')

P2 <- ggplot(ACC, aes(steps, Accuracy)) +
  stat_smooth(color = 'red') +
  geom_point() +
  geom_errorbar(aes(ymin = AccuracyLower, ymax = AccuracyUpper)) +
  xlab('Cells used for Training') +
  ylab('Symphony Accuracy (%)') +
  theme_bw()
print(P2)


# Creating the atlas for BRCA
clBRCA <- buildReference(
  exp_ref = X@assays$RNA@counts,
  metadata_ref = data.frame(g = 1, cellLine = getCellLines(X@assays$RNA@counts)),
  vars = 'g',
  do_umap = TRUE,
  verbose = TRUE,
  d = 50, save_uwot_path = 'umapBRCA')

umapBRCA <- data.frame(clBRCA$umap$embedding, cl=clBRCA$meta_data$cellLine)
save(clBRCA, umapBRCA, file = '../Results/refBRCA.RData')

labelPos <- split(umapBRCA, umapBRCA$cl)
labelPos <- lapply(labelPos, function(S){
  mDistance <- mahalanobis(S[,1:2], colMeans(S[,1:2]), cov(S[,1:2]))
  S <- S[!mDistance %in% boxplot.stats(mDistance)$out,]
  as.data.frame(c(apply(S[,1:2],2,median), CL = S[1,3]))
})
labelPos <- t(do.call(cbind.data.frame, labelPos))
labelPos <- as.data.frame(labelPos)
rownames(labelPos) <- NULL
labelPos$UMAP1 <- as.numeric(labelPos$UMAP1)
labelPos$UMAP2 <- as.numeric(labelPos$UMAP2)
labelPos$cl <- factor(labelPos$CL, levels = clMD$CL)
umapBRCA <- umapBRCA[order(umapBRCA$cl, decreasing = TRUE),]
umapBRCA$cl <- factor(umapBRCA$cl, levels = clMD$CL)
write.csv('umapBRCA', '../Results/F2A.csv')

P1 <- ggplot(umapBRCA, aes(UMAP1, UMAP2, color = cl)) +
  geom_point(cex = 0.01) +
  theme_bw() +
  theme(legend.position = 'None') +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  #scale_color_viridis_d() +
  geom_text_repel(aes(UMAP1, UMAP2, label = CL),
                  labelPos,
                  min.segment.length = 0,
                  nudge_y = 1.5,
                  bg.color = 'white')
print(P1)

# Testing MCF7
MCF7 <- read.csv('../Data/MCF7.csv')
MCF7 <- as.matrix(MCF7)

mcf7Map <- mapQuery(MCF7,
                    metadata_query = data.frame(rep('MCF7', ncol(MCF7))),
                    ref_obj = clBRCA, do_umap = TRUE)

mcf7Map <- knnPredict(mcf7Map, clBRCA, clBRCA$meta_data$cellLine, k = 5)
mcf7CM <- confusionMatrix(as.factor(mcf7Map$meta_data[,1]), as.factor(mcf7Map$meta_data[,2]))

mDistance <- mahalanobis(mcf7Map$umap, colMeans(mcf7Map$umap), cov(mcf7Map$umap))
plotData <- rbind(data.frame(clBRCA$umap$embedding, cl='Ref'), data.frame(mcf7Map$umap, cl='MCF7'))
plotData$ct <- c(clBRCA$meta_data$cellLine, ifelse(!mDistance %in% boxplot.stats(mDistance)$out, 'Q', 'O'))
write.csv(plotData, '../Results/F2C.csv')

P3 <- ggplot(plotData, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01,
             color = ifelse(plotData$ct %in% c('Q', 'O'), rgb(1,0,0,0.01), 'gray75'),
             alpha = 1) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_mark_ellipse(aes(filter = ct == 'Q', color = 'red'),expand = unit(2,'mm')) +
  annotate(x = -6, y = -7,
           geom = 'text', label = paste0(round(mcf7CM$overall[1]*100,1), '% MCF7'),
           color = 'red', family = 2, size = 3) +
  xlab('UMAP 1') +
  ylab('UMAP 2')
P3 <- P3 + labs(title = 'Wild-Type MCF7', subtitle = parse(text = 'italic(n)==14372~Cells')) +
  theme(plot.title = element_text(face = 2))
print(P3)

# Testing T47D
td47dData <- read.csv('../Data/GSM4285803_scRNA_RawCounts.csv.gz', row.names = 1)
td47dData <- t(td47dData)
td47dMetaData <- read.csv('../Data/GSM4285803_scRNA_metaInfo.csv.gz')
td47dData <- td47dData[,td47dMetaData$X[grepl('T47D KO', td47dMetaData$CellType)]]
td47dData <- as.matrix(td47dData)
td47Map <- mapQuery(exp_query = td47dData,
                    metadata_query = data.frame(rep('T47D', ncol(td47dData))),
                    ref_obj = clBRCA, do_umap = TRUE)
td47Map <- knnPredict(td47Map, clBRCA, clBRCA$meta_data$cellLine, k = 5)
td47CM <- confusionMatrix(as.factor(td47Map$meta_data[,1]), as.factor(td47Map$meta_data[,2]))
td47CM

plotData <- rbind(data.frame(clBRCA$umap$embedding, cl='Ref'), data.frame(td47Map$umap, cl='T47D'))
plotData$ct <- 'Q'
plotData$ct[seq_len(nrow(clBRCA$meta_data))] <- c(clBRCA$meta_data$cellLine)
write.csv(plotData, '../Results/F2D.csv')

P4 <- ggplot(plotData, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01, color = ifelse(plotData$ct == 'Q', 'red', 'gray75'), alpha = 1) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_mark_ellipse(aes(filter = ct == 'T47D', color = 'red'),expand = unit(2,'mm')) +
  annotate(x = median(td47Map$ umap[,1]), y = 6,
           geom = 'text', label = paste0(td47CM$overall[1]*100, '% T47D'),
           color = 'red', family = 2, size = 3) +
  xlab('UMAP 1') +
  ylab('UMAP 2')
P4 <- P4 + labs(title = 'T47D with CDH1 knockout',
          subtitle = parse(text = 'italic(n)==491~Cells')) +
  theme(plot.title = element_text(face = 2))
print(P4)

# Testing BT474
bt474Data <- read.csv('../Data/GSE150949_pooled_watermelon.count.matrix.csv.gz', row.names = 1)
bt474MetaData <- read.csv('../Data/GSE150949_pooled_watermelon.metadata.matrix.csv.gz')
bt474MetaData <- bt474MetaData[grepl('BT474', bt474MetaData$cell_line),]
bt474Data <- bt474Data[,gsub('-','.',bt474MetaData$cell)]
bt474Data <- as.matrix(bt474Data)
bt474Map <- mapQuery(exp_query = bt474Data,
                    metadata_query = data.frame(rep('BT474', ncol(bt474Data))),
                    ref_obj = clBRCA, do_umap = TRUE)
bt474Map <- knnPredict(bt474Map, clBRCA, clBRCA$meta_data$cellLine, k = 5)
bt474CM <- confusionMatrix(as.factor(bt474Map$meta_data[,1]), as.factor(bt474Map$meta_data[,2]))
bt474CM

mDistance <- mahalanobis(bt474Map$umap, colMeans(bt474Map$umap), cov(bt474Map$umap))
plotData <- rbind(data.frame(clBRCA$umap$embedding, cl='Ref'), data.frame(bt474Map$umap, cl='BT474'))
plotData$ct <- c(clBRCA$meta_data$cellLine, ifelse(!mDistance %in% boxplot.stats(mDistance)$out, 'Q', 'O'))
write.csv(plotData, '../Results/F2E.csv')

P5 <- ggplot(plotData, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01,
             color = ifelse(plotData$ct %in% c('Q', 'O'), 'red', 'gray75'),
             alpha = 1) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_mark_ellipse(aes(filter = ct == 'Q', color = 'red'),expand = unit(2,'mm')) +
    annotate(x = -9, y = -1,
             geom = 'text', label = paste0(round(bt474CM$overall[1]*100,1), '% BT474'),
             color = 'red', family = 2, size = 3) +
  xlab('UMAP 1') +
  ylab('UMAP 2')
P5 <- P5 + labs(title = 'BT474 + Lapatinib', subtitle = parse(text = 'italic(n)==131~Cells')) +
  theme(plot.title = element_text(face = 2))
print(P5)

# Figure 1
png('../Figures/F2.png', width = 3500, height = 1750, res = 300)
pLayout <- '
AABC
AADE'
P1 + P2 + P3 + P4 + P5 +  plot_layout(design = pLayout) + plot_annotation(tag_levels = 'A', theme = theme(plot.tag = element_text(face = 2)))
dev.off()


# Patients data
query <- readMM('../Data/count_matrix_sparse.mtx')
rownames(query) <- readLines('../Data/count_matrix_genes.tsv')
colnames(query) <- readLines('../Data/count_matrix_barcodes.tsv')
queryMetadata <- read.csv('../Data/metadata.csv', row.names = 1)
queryMetadata <- queryMetadata[grepl('Cancer',queryMetadata$celltype_minor),]
query <- query[,rownames(queryMetadata)]
donorSubType <- queryMetadata$subtype
names(donorSubType) <- queryMetadata$orig.ident

clProportion <- lapply(unique(queryMetadata$orig.ident), function(donor){
  donorData <- query[,(queryMetadata$orig.ident %in% donor)]
  donorMD <- data.frame(donor = rep(donor, ncol(donorData)))
  qmap <- mapQuery(donorData, metadata_query = donorMD, ref_obj = clBRCA, do_umap = TRUE, do_normalize = TRUE)
  qmap <- knnPredict(qmap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
  qc <- qmap$meta_data
  qc <- table(qc$cell_type_pred_knn)/nrow(qc)
  qc <- data.frame(donor = donor, qc)
  colnames(qc) <- c('donor', 'cellLine', 'proportion')
  qc
})

clProportion <- do.call(rbind.data.frame, clProportion)
clProportion$subtype <- donorSubType[clProportion$donor]
clProportion$cellsubtype <- cellSubtype[as.vector(clProportion$cellLine)]
clProportion$cellsubtype[grepl('TNBC', clProportion$cellsubtype)] <- 'TNBC'
clProportion$cellsubtype[grepl('Lum', clProportion$cellsubtype)] <- 'ER+'
clProportion$cellLine <- factor(clProportion$cellLine, levels = clMD$CL)
write.csv(clProportion, '../Results/F3C.csv')

clProportion <- read.csv('../Results/F3C.csv', row.names = 1)
O <- round(acast(data = clProportion, formula = donor~cellLine, value.var = 'proportion') * 100,2)
O <- O[,clMD$CL]
col_fun = colorRamp2(c(0, 50), c("gray99", "red"))

png('../Figures/S1.png', width = 4000, height = 2500, res = 300)
Heatmap(O, col = col_fun,
        column_split = clMD$Type,name = '%',
        row_split = donorSubType[rownames(O)],
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.1f", O[i, j]), x, y, gp = gpar(fontsize = 10))})
dev.off()

chm <- acast(clProportion, subtype ~ cellLine, value.var = 'proportion', fun.aggregate = mean)
chm <- chm/rowSums(chm)
chm <- as.matrix(chm)
write.csv(chm, '../Results/F4.csv')

chm <- read.csv('../Results/F4.csv', row.names = 1)
chm <- chm * 100
col_fun = colorRamp2(c(0, 36.5), c("gray99", "red"))
#chm <- scale(chm)
chm <- round(chm,1)
png('../Figures/F4.png', width = 3500, height = 750, res = 300)
Heatmap(chm,
        column_split = cellSubtype[colnames(chm)],
        col = col_fun, name = '%',
        show_row_dend = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.1f", chm[i, j]), x, y, gp = gpar(fontsize = 10))})
dev.off()

pInfo <- data.frame(donor = clProportion$donor, subtype = clProportion$subtype)
pInfo <- unique(pInfo)
pInfo <- pInfo[order(pInfo$subtype),]
clProportion$donor <- factor(clProportion$donor,levels = pInfo$donor)
clProportion$proportion <- clProportion$proportion * 100
clProportion <- clProportion[order(clProportion$cellLine),]
clProportion$cellLine <- factor(clProportion$cellLine, levels = clMD$CL)

P7 <- ggplot(clProportion, aes(x = donor, y = proportion , fill = cellLine)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.title = element_text(face = 2)) +
  xlab('Donor') +
  ylab('%') +
  labs(fill = 'Reference\nCell Line')
print(P7)

# Donor specific profile
donor <- 'CID44971'
donorData <- query[,(queryMetadata$orig.ident %in% donor)]
donorMD <- data.frame(donor = rep(donor, ncol(donorData)))
qmap <- mapQuery(donorData, metadata_query = donorMD, ref_obj = clBRCA, do_umap = TRUE, do_normalize = TRUE)
qmap <- knnPredict(qmap, clBRCA, clBRCA$meta_data$cellLine, k = 5)
A <- data.frame(clBRCA$umap$embedding, cl = NA)
B <- data.frame(qmap$umap, cl = qmap$meta_data$cell_type_pred_knn)
qmap <- rbind(A,B)
qmap$cl <- factor(qmap$cl, levels = clMD$CL)
lInfo <- round((table(qmap$cl)/sum(table(qmap$cl)))*100,2)
lInfo <- paste0(lInfo, '% ', names(lInfo))
lInfo <- gsub('0% ', '', lInfo)
names(lInfo) <- names(table(qmap$cl))
labelPos$pct <- lInfo[as.vector(labelPos$cl)]
labelPos$pct[!grepl('%',labelPos$pct)] <- NA
write.csv(qmap, '../Results/F3A.csv')

P6 <- ggplot(qmap, aes(UMAP1, UMAP2)) +
  geom_point(cex = 0.01, color = ifelse(is.na(qmap$cl), 'gray75', rgb(1,0,0,1))) +
  theme_bw() +
  theme(legend.position = 'None') +
  geom_text_repel(aes(UMAP1, UMAP2, label = pct),
                  labelPos,
                  min.segment.length = 0,
                  nudge_y = 1.7,
                  bg.color = 'white',
                  col = 'black', size = 4) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  labs(title = donor,
       subtitle = parse(text = paste0('italic(n)==', ncol(donorData), '~Cancer~Cells'))) +
  theme(plot.title = element_text(face = 2))
P6

source('../Code/S2-LOOCV.R')
LOOCV <- read.csv('../Results/ccLOOCV.csv', row.names = 1)
cvMean <- round((apply(LOOCV,1,mean)/ncol(donorData)) * 100,2)
cvSD <- round((apply(LOOCV,1,sd)/ncol(donorData)) * 100,2)
RMSE <- sqrt(mean((round(sort(table(qc$cell_type_pred_knn)/nrow(qc)*100, decreasing = TRUE),2)-sort(cvMean, decreasing = TRUE))^2))
CTest <- cor.test(table(qc$cell_type_pred_knn)/nrow(qc)*100, cvMean, method = 'sp', continuity = TRUE)
CTest$p.value
DF <- data.frame(CL = rownames(LOOCV), M = cvMean, LB = cvMean-cvSD, UB = cvMean+cvSD)
DF <- DF[order(DF$M),]
DF$CL <- factor(DF$CL, levels = DF$CL)
DF$CL2 <- factor(DF$CL, levels = clMD$CL)
DF <- DF[DF$M > 0,]
write.csv(DF, '../Results/F3B.csv')

P8 <- ggplot(DF, aes(M,CL)) +
  geom_bar(stat = 'identity', fill = clColor[DF$CL2]) +
  geom_errorbarh(mapping = aes(xmin = LB, xmax = UB), height = .25, size = 0.5) +
  theme_bw() +
  theme(legend.position = 'None') +
  xlab('%') +
  ylab('Cell Line')
print(P8)


png('../Figures/F3.png', width = 3500, height = 2000, res = 300)
pLayout <- '
AAACCC
AAACCC
BBBCCC
'
P6 + P8 + P7 + plot_layout(design = pLayout) + plot_annotation(tag_levels = 'A', theme = theme(plot.tag = element_text(face = 2)))
dev.off()
