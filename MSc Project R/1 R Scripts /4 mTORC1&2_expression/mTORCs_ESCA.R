#### mTORC1 AND mTORC2 EXPRESSION DATA - ESCA ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.esca <- log2(expr.esca+1) # log of exp data 
expr.genes.esca <- expr.esca[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.esca <-pheatmap(t(expr.genes.esca), show_colnames = FALSE,
                  cutree_cols = 3) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.esca$tree_col
assignments.esca <- cutree(clusteredSamples, k=3) # k = cutree_cols
groupAssignments.esca <- data.frame(Group=factor(assignments.esca))
p.esca <- pheatmap(t(expr.genes.esca), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments.esca)

### Merging group assignments 
df.expr.esca <- data.frame(expr.genes.esca)
df.expr.esca$Sample <- rownames(df.expr.esca)
groupAssignments.esca$SampleID <- rownames(groupAssignments.esca)
df.merged.esca <- merge(df.expr.esca, groupAssignments.esca,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# merging clinical data
clinical.esca$PatientID <- rownames(clinical.esca)
df.merged.esca$PatientID <- substr(df.merged.esca$Sample, start = 1, stop = 12)
df.merged.esca$PatientAge <- clinical.esca[match(df.merged.esca$PatientID, clinical.esca$PatientID), "yearstobirth"]
df.merged.esca$PatientAge <- as.numeric(df.merged.esca$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.esca$SubID <- substr(df.merged.esca$Sample, start = 1, stop = 19)
S1_corr_data_ESCA_filter$SubID <- substr(S1_corr_data_ESCA_filter$Sample_ID, start = 1, stop = 19)
df.merged.esca <- merge(df.merged.esca, S1_corr_data_ESCA_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L", "lowDEPDC6.lowPRR5L")
df.merged.esca$Status <- values[df.merged.esca$Group]

## MP ##
library(ggpubr)
violin.esca.mp <- ggviolin(df.merged.esca, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#f29b30", "#5dd4f5"), add = "boxplot", add.params = list(fill = "white"))
comparison.esca <- list(c("highDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L"), c("lowDEPDC6.highPRR5L", "lowDEPDC6.lowPRR5L"), c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L"))

p.mp <- violin.esca.mp + stat_compare_means(comparisons = comparison.esca) + 
        stat_compare_means(label.y = 0.35, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.esca.mp <- ggpar(p.mp, 
                 main = "ESCA",
                 font.main = c(16, "bold"),
                 xlab = "Gene expression status",
                 ylab = "% S1 (MP)",
                 font.x = c(14, "bold"), 
                 font.y = c(14, "bold"),
                 font.ytickslab = 14,
                 font.xtickslab = c(1, "white"),
                 legend = "none")
mean.ds <- aggregate(df.merged.esca$DS ~ df.merged.esca$Status, df.merged.esca, mean)

## DS ##
violin.esca.ds <- ggviolin(df.merged.esca, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#f27735", "#1155f5"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.esca.ds + stat_compare_means(comparisons = comparison.esca) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.esca.ds <- ggpar(p.ds, 
                     main = "ESCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

mean.ds <- aggregate(df.merged.esca$DS ~ df.merged.esca$Status, df.merged.esca, mean)

save(df.merged.esca, file = "df.merged.esca.RData")
