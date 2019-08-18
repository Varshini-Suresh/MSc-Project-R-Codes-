#### mTORC1 AND mTORC2 EXPRESSION DATA - LUAD ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.luad <- log2(expr.luad+1) # log of exp data 
expr.genes.luad <- expr.luad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.luad <-pheatmap(t(expr.genes.luad), show_colnames = FALSE,
                  cutree_cols = 9) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.luad$tree_col
assignments.luad <- cutree(clusteredSamples, k=9) # k = cutree_cols
groupAssignments.luad <- data.frame(Group=factor(assignments.luad))
p.luad <- pheatmap(t(expr.genes.luad), show_colnames = FALSE,
                   cutree_cols = 9, annotation = groupAssignments.luad)

### Merging group assignments 
df.expr.luad <- data.frame(expr.genes.luad)
df.expr.luad$Sample <- rownames(df.expr.luad)
groupAssignments.luad$SampleID <- rownames(groupAssignments.luad)
df.merged.luad <- merge(df.expr.luad, groupAssignments.luad,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.luad.copy <- df.merged.luad
filter.luad <- subset(df.merged.luad, df.merged.luad$Group != "9")
filter.luad$Group[filter.luad$Group =="2"] <-"1" # gp1 = 1,2,4,5,8 and gp2 = 3,6,7
df.merged.luad <- filter.luad

# merging clinical data
clinical.luad$PatientID <- rownames(clinical.luad)
df.merged.luad$PatientID <- substr(df.merged.luad$Sample, start = 1, stop = 12)
df.merged.luad$PatientAge <- clinical.luad[match(df.merged.luad$PatientID, clinical.luad$PatientID), "yearstobirth"]
df.merged.luad$PatientAge <- as.numeric(df.merged.luad$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.luad$SubID <- substr(df.merged.luad$Sample, start = 1, stop = 19)
S1_corr_data_LUAD_filter$SubID <- substr(S1_corr_data_LUAD_filter$SampleID, start = 1, stop = 19)
df.merged.luad <- merge(df.merged.luad, S1_corr_data_LUAD_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L")
df.merged.luad$Status <- values[df.merged.luad$Group]

## MP ##
library(ggpubr)
violin.luad.mp <- ggviolin(df.merged.luad, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#5dd4f5"), add = "boxplot", add.params = list(fill = "white"))
comparison.luad <- list(c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L"))

p.mp <- violin.luad.mp + stat_compare_means(comparisons = comparison.luad) + 
  stat_compare_means(label.y = 0.32, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.luad.mp <- ggpar(p.mp, 
                     main = "LUAD",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.luad$MP ~ df.merged.luad$Status, df.merged.luad, mean)

## DS ##
violin.luad.ds <- ggviolin(df.merged.luad, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#1155f5"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.luad.ds + stat_compare_means(comparisons = comparison.luad) + 
  stat_compare_means(label.y = 1.2, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.luad.ds <- ggpar(p.ds, 
                     main = "LUAD",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.luad$DS ~ df.merged.luad$Status, df.merged.luad, mean)

save(df.merged.luad, file = "df.merged.luad.RData")
