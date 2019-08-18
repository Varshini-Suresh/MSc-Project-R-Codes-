#### mTORC1 AND mTORC2 EXPRESSION DATA - BLCA ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.blca <- log2(expr.blca+1) # log of exp data 
expr.genes.blca <- expr.blca[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.blca <-pheatmap(t(expr.genes.blca), show_colnames = FALSE,
                  cutree_cols = 1) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.blca$tree_col
assignments.blca <- cutree(clusteredSamples, k=5) # k = cutree_cols
groupAssignments.blca <- data.frame(Group=factor(assignments.blca))
p.blca <- pheatmap(t(expr.genes.blca), show_colnames = FALSE,
                   cutree_cols = 5, annotation = groupAssignments.blca)

### Merging group assignments 
df.expr.blca <- data.frame(expr.genes.blca)
df.expr.blca$Sample <- rownames(df.expr.blca)
groupAssignments.blca$SampleID <- rownames(groupAssignments.blca)
df.merged.blca <- merge(df.expr.blca, groupAssignments.blca,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.blca$Group[df.merged.blca$Group =="5"] <-"1"

# merging clinical data
clinical.blca$PatientID <- rownames(clinical.blca)
df.merged.blca$PatientID <- substr(df.merged.blca$Sample, start = 1, stop = 12)
df.merged.blca$PatientAge <- clinical.blca[match(df.merged.blca$PatientID, clinical.blca$PatientID), "yearstobirth"]
df.merged.blca$PatientAge <- as.numeric(df.merged.blca$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.blca$SubID <- substr(df.merged.blca$Sample, start = 1, stop = 19)
S1_corr_data_BLCA$SubID <- substr(S1_corr_data_BLCA$SampleID, start = 1, stop = 19)
df.merged.blca <- merge(df.merged.blca, S1_corr_data_BLCA,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("lowDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L", "highDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L")
df.merged.blca$Status <- values[df.merged.blca$Group]

## MP ##
library(ggpubr)
boxplot.mp <- ggboxplot(df.merged.blca, x="Status", y="MP", fill = "Status", palette = c("#5dd4f5", "#f291d0", "#95d108", "#f29b30"), shape = "Status")
comparison.blca <- list(c("highDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L"), c("lowDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L"), c("highDEPDC6.highPRR5L", "highDEPDC6.lowPRR5L"), 
                        c("highDEPDC6.highPRR5L", "lowDEPDC6.highPRR5L"), c("lowDEPDC6.lowPRR5L", "highDEPDC6.lowPRR5L"), c("lowDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L"))

p.mp <- boxplot.mp + stat_compare_means(comparisons = comparison.blca) + 
  stat_compare_means(label.y = 0.35, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.blca.mp <- ggpar(p.mp, 
                     main = "BLCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.blca$MP ~ df.merged.blca$Status, df.merged.blca, mean)

## DS ##
boxplot.ds <- ggboxplot(df.merged.blca, x="Status", y="DS", fill = "Status", palette = c("#5dd4f5", "#f291d0", "#95d108", "#f29b30"), shape = "Status")
p.ds <- boxplot.ds + stat_compare_means(comparisons = comparison.blca) + 
  stat_compare_means(label.y = 1.75, size = 6)
p.ds$layers[[2]]$aes_params$textsize <- 6

lab.blca.ds <- ggpar(p.ds, 
                     main = "BLCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.blca$DS ~ df.merged.blca$Status, df.merged.blca, mean)

save(df.merged.blca, file = "df.merged.blca.RData")
