#### mTORC1 AND mTORC2 EXPRESSION DATA - BRCA ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.brca <- log2(expr.brca+1) # log of exp data 
expr.genes.brca <- expr.brca[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.brca <-pheatmap(t(expr.genes.brca), show_colnames = FALSE,
                  cutree_cols = 2) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.brca$tree_col
assignments.brca <- cutree(clusteredSamples, k=2) # k = cutree_cols
groupAssignments.brca <- data.frame(Group=factor(assignments.brca))
p.brca <- pheatmap(t(expr.genes.brca), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments.brca)

### Merging group assignments 
df.expr.brca <- data.frame(expr.genes.brca)
df.expr.brca$Sample <- rownames(df.expr.brca)
groupAssignments.brca$SampleID <- rownames(groupAssignments.brca)
df.merged.brca <- merge(df.expr.brca, groupAssignments.brca,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# merging clinical data
clinical.brca$PatientID <- rownames(clinical.brca)
df.merged.brca$PatientID <- substr(df.merged.brca$Sample, start = 1, stop = 12)
df.merged.brca$PatientAge <- clinical.brca[match(df.merged.brca$PatientID, clinical.brca$PatientID), "yearstobirth"]
df.merged.brca$PatientAge <- as.numeric(df.merged.brca$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.brca$SubID <- substr(df.merged.brca$Sample, start = 1, stop = 19)
S1_corr_data_BRCA$SubID <- substr(S1_corr_data_BRCA$SampleID, start = 1, stop = 19)
df.merged.brca <- merge(df.merged.brca, S1_corr_data_BRCA,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("lowDEPDC6.lowPRR5L", "highDEPDC6.lowPRR5L")
df.merged.brca$Status <- values[df.merged.brca$Group]

## MP ##
library(ggpubr)
violin.brca.mp <- ggviolin(df.merged.brca, x = "Status", y = "MP", fill = "Status", palette = c("#5dd4f5", "#95d108"), add = "boxplot", add.params = list(fill = "white"))
comparison.brca <- list(c("lowDEPDC6.lowPRR5L", "highDEPDC6.lowPRR5L"))

p.mp <- violin.brca.mp + stat_compare_means(comparisons = comparison.brca) + 
  stat_compare_means(label.y = 0.55, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.brca.mp <- ggpar(p.mp, 
                     main = "BRCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

## DS ##
violin.brca.ds <- ggviolin(df.merged.brca, x = "Status", y = "DS", fill = "Status", palette = c("#1155f5", "#2da10a"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.brca.ds + stat_compare_means(comparisons = comparison.brca) + 
  stat_compare_means(label.y = 1.4, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.brca.ds <- ggpar(p.ds, 
                     main = "BRCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

save(df.merged.brca, file = "df.merged.brca.RData")
