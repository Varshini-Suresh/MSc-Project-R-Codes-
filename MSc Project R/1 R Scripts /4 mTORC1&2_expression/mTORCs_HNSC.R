#### mTORC1 AND mTORC2 EXPRESSION DATA - HNSC ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.hnsc <- log2(expr.hnsc+1) # log of exp data 
expr.genes.hnsc <- expr.hnsc[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.hnsc <-pheatmap(t(expr.genes.hnsc), show_colnames = FALSE,
                  cutree_cols = 4) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.hnsc$tree_col
assignments.hnsc <- cutree(clusteredSamples, k=4) # k = cutree_cols
groupAssignments.hnsc <- data.frame(Group=factor(assignments.hnsc))
p.hnsc <- pheatmap(t(expr.genes.hnsc), show_colnames = FALSE,
                   cutree_cols = 4, annotation = groupAssignments.hnsc)

### Merging group assignments 
df.expr.hnsc <- data.frame(expr.genes.hnsc)
df.expr.hnsc$Sample <- rownames(df.expr.hnsc)
groupAssignments.hnsc$SampleID <- rownames(groupAssignments.hnsc)
df.merged.hnsc <- merge(df.expr.hnsc, groupAssignments.hnsc,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.hnsc$Group[df.merged.hnsc$Group =="4"] <-"3"

# merging clinical data
clinical.hnsc$PatientID <- rownames(clinical.hnsc)
df.merged.hnsc$PatientID <- substr(df.merged.hnsc$Sample, start = 1, stop = 12)
df.merged.hnsc$PatientAge <- clinical.hnsc[match(df.merged.hnsc$PatientID, clinical.hnsc$PatientID), "yearstobirth"]
df.merged.hnsc$PatientAge <- as.numeric(df.merged.hnsc$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.hnsc$SubID <- substr(df.merged.hnsc$Sample, start = 1, stop = 19)
S1_corr_data_HNSC$SubID <- substr(S1_corr_data_HNSC$SampleID, start = 1, stop = 19)
df.merged.hnsc <- merge(df.merged.hnsc, S1_corr_data_HNSC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("lowDEPDC6.highPRR5L", "lowDEPDC6.lowPRR5L", "highDEPDC6.lowPRR5L")
df.merged.hnsc$Status <- values[df.merged.hnsc$Group]

## MP ##
library(ggpubr)
violin.hnsc.mp <- ggviolin(df.merged.hnsc, x = "Status", y = "MP", fill = "Status", palette = c("#f29b30", "#5dd4f5", "#95d108"), add = "boxplot", add.params = list(fill = "white"))
comparison.hnsc <- list(c("lowDEPDC6.highPRR5L", "lowDEPDC6.lowPRR5L"), c("lowDEPDC6.lowPRR5L", "highDEPDC6.lowPRR5L"), c("lowDEPDC6.highPRR5L", "highDEPDC6.lowPRR5L"))

p.mp <- violin.hnsc.mp + stat_compare_means(comparisons = comparison.hnsc) + 
  stat_compare_means(label.y = 0.45, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.hnsc.mp <- ggpar(p.mp, 
                     main = "HNSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.hnsc$MP ~ df.merged.hnsc$Status, df.merged.hnsc, mean)

## DS ##
violin.hnsc.ds <- ggviolin(df.merged.hnsc, x = "Status", y = "DS", fill = "Status", palette = c("#f27735", "#1155f5", "#2da10a"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.hnsc.ds + stat_compare_means(comparisons = comparison.hnsc) + 
  stat_compare_means(label.y = 1.6, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.hnsc.ds <- ggpar(p.ds, 
                     main = "HNSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.hnsc$DS ~ df.merged.hnsc$Status, df.merged.hnsc, mean)

save(df.merged.hnsc, file = "df.merged.hnsc.RData")
