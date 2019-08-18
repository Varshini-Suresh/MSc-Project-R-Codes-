#### mTORC1 AND mTORC2 EXPRESSION DATA - OV ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.ov <- log2(expr.ov+1) # log of exp data 
expr.genes.ov <- expr.ov[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.ov <-pheatmap(t(expr.genes.ov), show_colnames = FALSE,
                  cutree_cols = 3) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.ov$tree_col
assignments.ov <- cutree(clusteredSamples, k=3) # k = cutree_cols
groupAssignments.ov <- data.frame(Group=factor(assignments.ov))
p.ov <- pheatmap(t(expr.genes.ov), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments.ov)

### Merging group assignments 
df.expr.ov <- data.frame(expr.genes.ov)
df.expr.ov$Sample <- rownames(df.expr.ov)
groupAssignments.ov$SampleID <- rownames(groupAssignments.ov)
df.merged.ov <- merge(df.expr.ov, groupAssignments.ov,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.ov$Group[df.merged.ov$Group =="2"] <-"1" # gp1 = 1,2 and gp2 = 3

# merging clinical data
clinical.ov$PatientID <- rownames(clinical.ov)
df.merged.ov$PatientID <- substr(df.merged.ov$Sample, start = 1, stop = 12)
df.merged.ov$PatientAge <- clinical.ov[match(df.merged.ov$PatientID, clinical.ov$PatientID), "yearstobirth"]
df.merged.ov$PatientAge <- as.numeric(df.merged.ov$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.ov$SubID <- substr(df.merged.ov$Sample, start = 1, stop = 19)
S1_corr_data_OV_filter$SubID <- substr(S1_corr_data_OV_filter$SampleID, start = 1, stop = 19)
df.merged.ov <- merge(df.merged.ov, S1_corr_data_OV_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L")
df.merged.ov$Status <- values[df.merged.ov$Group]

## MP ##
library(ggpubr)
violin.ov.mp <- ggviolin(df.merged.ov, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#5dd4f5"), add = "boxplot", add.params = list(fill = "white"))
comparison.ov <- list(c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L"))

p.mp <- violin.ov.mp + stat_compare_means(comparisons = comparison.ov) + 
  stat_compare_means(label.y = 0.20, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.ov.mp <- ggpar(p.mp, 
                     main = "OV",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.ov$MP ~ df.merged.ov$Status, df.merged.ov, mean)

## DS ##
violin.ov.ds <- ggviolin(df.merged.ov, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#1155f5"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.ov.ds + stat_compare_means(comparisons = comparison.ov) + 
  stat_compare_means(label.y = 1.0, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.ov.ds <- ggpar(p.ds, 
                     main = "OV",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.ov$DS ~ df.merged.ov$Status, df.merged.ov, mean)

save(df.merged.ov, file = "df.merged.ov.RData")
