library(ggplot2)
library(ggrepel)
library(circlize)
library(patchwork)
library(ggpubr)
library(fgsea)

cellLineType <- list()
cellLineType$ER <- c('EFM19', 'BT483', 'CAL148', 'UACC812', 'HCC1428', 'MFM223','MCF7',
                     'DU4475', 'CAMA1', 'HCC1500')
cellLineType$HER2 <- c('MDAMB175VII', 'MDAMB361', 'HCC202', 'UACC893', 'HCC2218', 'HCC1419',
                       'AU565', 'EFM192A', 'MDAMB415', 'MDAMB453', 'ZR7530', 'BT474')
cellLineType$TNBC <- c('JIMT1', 'BT20', 'HDQP1', 'CAL851', 'MDAMB231', 'HCC70',
                       'HCC1599', 'CAL51', 'HCC1937', 'CAL120', 'HCC38', 'MDAMB468', 'HCC1806',
                       'BT549', 'HCC2157', 'HCC1954', 'HCC1187', 'HCC1143', 'HCC1395', 'MDAMB436')

iData <- read.csv('../Data/Original screen_All tissues_fitted data.csv')
iData <- iData[iData$Tissue == 'Breast',]
iData$CELL_LINE_NAME <- toupper(gsub('-','',iData$CELL_LINE_NAME))

cancerSubtypes <- read.csv('../Results/F4.csv', row.names = 1)
BRCA <- colSums(cancerSubtypes)/sum(cancerSubtypes)

patientComposition <- read.csv('../Results/F3C.csv', row.names = 1)
CID44971 <- patientComposition$proportion[patientComposition$donor == 'CID44971']
names(CID44971) <- patientComposition$cellLine[patientComposition$donor == 'CID44971']

findDrugCombinations <- function(iData, tumorComposition, plotTitle = NA, lSize = 2.2, nLabels = 5, outFile){
  tumorComposition <- tumorComposition[tumorComposition != 0]
  fData <- iData[iData$CELL_LINE_NAME %in% names(tumorComposition),]
  fData$COMB <- paste0(fData$ANCHOR_CONC, ' ', fData$ANCHOR_NAME, ' + ', fData$LIBRARY_CONC, ' ', fData$LIBRARY_NAME)
  fData$COMB_COMP <- paste0(fData$ANCHOR_NAME, ' + ', fData$LIBRARY_NAME)
  fData <- split(fData, fData$CELL_LINE_NAME)
  for(i in seq_along(fData)){
    fData[[i]]$SYNERGY_DELTA_EMAX <- fData[[i]]$SYNERGY_DELTA_EMAX * tumorComposition[fData[[i]]$CELL_LINE_NAME]
  }
  fData <- do.call(rbind.data.frame, fData)
  cDrugs <- unlist(lapply(split(fData$SYNERGY_DELTA_EMAX, fData$COMB_COMP), mean))
  fDrugs <- unlist(lapply(split(fData$CELL_LINE_NAME, fData$COMB_COMP), function(X){sum(tumorComposition[unique(X)])}))
  rmseComb <- unlist(lapply(split(fData$SYNERGY_RMSE, fData$COMB_COMP), mean))
  rDrugs <- (cDrugs * fDrugs)/rmseComb
  B <- abs(sapply(1:1e3, function(X){sample(rDrugs, replace = TRUE)}))
  P <- sapply(abs(rDrugs), function(X){mean(X <= B)})
  colFun <- circlize::colorRamp2(c(min(rDrugs),0, max(rDrugs)), colors = c('gray90','gray80', 'red'))
  df <- data.frame(sActivity = cDrugs*fDrugs, P = P, id = names(cDrugs))
  df$id[!df$id %in% names(sort(rDrugs, decreasing = TRUE)[seq_len(nLabels)])] <- NA
  dPCH <- ifelse(((df$sActivity > 0) & (df$P < 0.05)),8,16)
  df$sActivity <- df$sActivity/max(abs(df$sActivity))
  write.csv(df, outFile)
  n = sum(df$sActivity > 0)
  oPlot <- ggplot(df, aes(sActivity, -log10(P), label = id)) +
    geom_point(color = colFun(rDrugs), cex = 0.5, pch = dPCH) +
    geom_text_repel(min.segment.length = 0,
                    size = lSize,
                    bg.color = 'white',
                    segment.size = 0.1,
                    #nudge_x = .15,
                    box.padding = 0.25,
                    seed = 123,
                    max.time = 1,
                    max.iter = Inf
                    #nudge_y = 0.1,
                    #segment.curvature = -0.15
                    #segment.ncp = 3,
                    #segment.angle = 5
                    ) +
    theme_bw() +
    xlab('Activity Score') +
    ylab(parse(text = '-log[10]~(P[Bootstrap])')) +
    labs(title = plotTitle, subtitle = parse(text = paste0('italic(n)==',n,'~Synergistic~Drug~Combinations'))) +
    xlim(c(0, max(df$sActivity))) +
    theme(plot.title = element_text(face = 2), plot.subtitle = element_text(size = 7))
  return(oPlot)
}

P1 <- findDrugCombinations(iData, CID44971, plotTitle = 'CID44971', lSize = 3, nLabels = 10, outFile = '../Results/F5A.csv')
P2 <- findDrugCombinations(iData, BRCA, plotTitle = 'BRCA', outFile = '../Results/F5B.csv')
P3 <- findDrugCombinations(iData, unlist(cancerSubtypes[1,]), plotTitle = 'ER+', outFile = '../Results/F5C.csv')
P4 <- findDrugCombinations(iData, unlist(cancerSubtypes[2,]), plotTitle = 'HER2+', outFile = '../Results/F5D.csv')
P5 <- findDrugCombinations(iData, unlist(cancerSubtypes[3,]), plotTitle = 'TNBC', outFile = '../Results/F5E.csv')

A <- lapply(c('../Results/F5C.csv', '../Results/F5D.csv', '../Results/F5E.csv'), function(X){
  X <- read.csv(X)
  #X <- X[X[,2] > 0,]
  A <- X[,2]
  names(A) <- X[,1]
  return(A)
})
cNames <- unique(unlist(lapply(A, names)))
A <- sapply(A, function(X){X[cNames]})
colnames(A) <- c('ER+', 'HER2+', 'TNBC')
A <- reshape2::melt(A)
colnames(A) <- c('Combination', 'Subtype', 'Activity')
A <-A[is.finite(A$Activity),]
write.csv(A, '../Results/F5F.csv')

colFun <- circlize::colorRamp2(c(0, 1), colors = c('gray90', 'red'))
A$Color <- colFun(A$Activity)
A <- A[A$Activity>0,]
P6 <- ggplot(A, aes(Activity, Subtype)) +
  geom_jitter(height = 0.25, cex = 0.1, color = A$Color) +
  geom_boxplot(outlier.colour = NA, fill = NA) +
  theme_bw() +
  xlab('Activity Score')
P6 <- P6 +
  stat_compare_means(comparisons = list(c('TNBC', 'HER2+'),c('TNBC','ER+')),
                     method = 'wilcox.test',
                     label = 'p.signif',
                     label.y = c(0.9, 0.95),
                     bracket.size = 0.2, tip.length = 0.005)
P6
pLayout <- '
AAAABBCC
AAAABBCC
AAAADDEE
AAAADDEE
FFFFFFFF'

png('../Figures/F5.png', width = 3000, height = 1800, res = 300)
P1 + P2 + P3 + P4 + P5 + P6 + plot_layout(design = pLayout) + plot_annotation(tag_levels = 'A')
dev.off()

iData$COMB_COMP <- paste0(iData$ANCHOR_NAME, ' + ', iData$LIBRARY_NAME)
nData <- iData[grepl('Navitoclax',iData$COMB_COMP),]
cLmedian <- sort(unlist(lapply(split(nData$SYNERGY_DELTA_EMAX,nData$CELL_LINE_NAME), median)), decreasing = TRUE)
nData$CELL_LINE_NAME <- factor(nData$CELL_LINE_NAME, levels = names(cLmedian))

set.seed(1)
E <- fgsea::fgseaMultilevel(cellLineType, cLmedian)
E

PA <- plotEnrichment(cellLineType$ER, cLmedian) +
  ylab('Enrichment Score') +
  xlab('Cell Line Rank') +
  theme_bw() +
  labs(title = 'ER+',
       subtitle = parse(text = paste0('NES == ',round(E$NES[1],2),'~~P-adj==',formatC(E$padj[1],digits = 2, format = 'g')))) +
  theme(plot.title = element_text(face = 2))
PB <- plotEnrichment(cellLineType$HER2, cLmedian) +
  ylab('Enrichment Score') +
  xlab('Cell Line Rank') +
  theme_bw() +
  labs(title = 'HER2+',
       subtitle = parse(text = paste0('NES == ',round(E$NES[2],2),'~~P-adj==',formatC(E$padj[2],digits = 2, format = 'g')))) +
  theme(plot.title = element_text(face = 2))
PC <- plotEnrichment(cellLineType$TNBC, cLmedian) +
  ylab('Enrichment Score') +
  xlab('Cell Line Rank') +
  theme_bw()+
  labs(title = 'TNBC',
       subtitle = parse(text = paste0('NES == ',round(E$NES[3],2),'~~P-adj==',formatC(E$padj[3],digits = 2, format = 'g')))) +
  theme(plot.title = element_text(face = 2))

PD <- ggplot(nData, aes(CELL_LINE_NAME,SYNERGY_DELTA_EMAX)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(width = 0.1, alpha = 0.05) +
  ylab(parse(text = 'Delta~Synergistic~E[max]')) +
  xlab('Cell Line') +
  labs(title = 'Navitoclax Combinations') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(face = 2))

pLayout <- '
AAAAAA
AAAAAA
AAAAAA
BBCCDD'

png('../Figures/S6.png', width = 3000, height = 2000, res = 300)
PD+PA+PB+PC + plot_layout(design = pLayout) + plot_annotation(tag_levels = 'A')
dev.off()

rER <- rep(1/length(cellLineType$ER), length(cellLineType$ER))
names(rER) <- cellLineType$ER
rHER2 <- rep(1/length(cellLineType$HER2), length(cellLineType$HER2))
names(rHER2) <- cellLineType$HER2
rTNBC <- rep(1/length(cellLineType$TNBC), length(cellLineType$TNBC))
names(rTNBC) <- cellLineType$TNBC

P7 <- findDrugCombinations(iData, rER, plotTitle = 'ER+', outFile = '../Results/SF5A.csv')
P8 <- findDrugCombinations(iData, rHER2, plotTitle = 'HER2+', outFile = '../Results/SF5B.csv')
P9 <- findDrugCombinations(iData, rTNBC, plotTitle = 'TNBC', outFile = '../Results/SF5C.csv')

png('../Figures/S5.png', width = 3000, height = 1000, res = 300)
P7 + P8 + P9
dev.off()
