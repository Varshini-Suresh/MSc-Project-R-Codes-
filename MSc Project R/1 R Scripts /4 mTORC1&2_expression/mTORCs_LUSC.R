#### mTORC1 AND mTORC2 EXPRESSION DATA - LUSC ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.lusc <- log2(expr.lusc+1) # log of exp data 
expr.genes.lusc <- expr.lusc[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.lusc <-pheatmap(t(expr.genes.lusc), show_colnames = FALSE,
                  cutree_cols = 4) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.lusc$tree_col
assignments.lusc <- cutree(clusteredSamples, k=4) # k = cutree_cols
groupAssignments.lusc <- data.frame(Group=factor(assignments.lusc))
p.lusc <- pheatmap(t(expr.genes.lusc), show_colnames = FALSE,
                   cutree_cols = 4, annotation = groupAssignments.lusc)

### Merging group assignments 
df.expr.lusc <- data.frame(expr.genes.lusc)
df.expr.lusc$Sample <- rownames(df.expr.lusc)
groupAssignments.lusc$SampleID <- rownames(groupAssignments.lusc)
df.merged.lusc <- merge(df.expr.lusc, groupAssignments.lusc,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.lusc$Group[df.merged.lusc$Group =="3"] <-"2" # gp1 = 1,4 and gp2 = 2,3

# merging clinical data
clinical.lusc$PatientID <- rownames(clinical.lusc)
df.merged.lusc$PatientID <- substr(df.merged.lusc$Sample, start = 1, stop = 12)
df.merged.lusc$PatientAge <- clinical.lusc[match(df.merged.lusc$PatientID, clinical.lusc$PatientID), "yearstobirth"]
df.merged.lusc$PatientAge <- as.numeric(df.merged.lusc$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.lusc$SubID <- substr(df.merged.lusc$Sample, start = 1, stop = 19)
S1_corr_data_LUSC$SubID <- substr(S1_corr_data_LUSC$SampleID, start = 1, stop = 19)
df.merged.lusc <- merge(df.merged.lusc, S1_corr_data_LUSC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("lowDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L")
df.merged.lusc$Status <- values[df.merged.lusc$Group]

## MP ##
library(ggpubr)
violin.lusc.mp <- ggviolin(df.merged.lusc, x = "Status", y = "MP", fill = "Status", palette = c("#5dd4f5", "#f29b30"), add = "boxplot", add.params = list(fill = "white"))
comparison.lusc <- list(c("lowDEPDC6.lowPRR5L", "lowDEPDC6.highPRR5L"))

p.mp <- violin.lusc.mp + stat_compare_means(comparisons = comparison.lusc) + 
  stat_compare_means(label.y = 0.32, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.lusc.mp <- ggpar(p.mp, 
                     main = "LUSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.lusc$MP ~ df.merged.lusc$Status, df.merged.lusc, mean)

## DS ##
violin.lusc.ds <- ggviolin(df.merged.lusc, x = "Status", y = "DS", fill = "Status", palette = c("#1155f5", "#f27735"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.lusc.ds + stat_compare_means(comparisons = comparison.lusc) + 
  stat_compare_means(label.y = 1., size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.lusc.ds <- ggpar(p.ds, 
                     main = "LUSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.lusc$DS ~ df.merged.lusc$Status, df.merged.lusc, mean)

save(df.merged.lusc, file = "df.merged.lusc.RData")
