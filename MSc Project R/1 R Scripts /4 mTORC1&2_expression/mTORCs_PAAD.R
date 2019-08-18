#### mTORC1 AND mTORC2 EXPRESSION DATA - PAAD ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.paad <- log2(expr.paad+1) # log of exp data 
expr.genes.paad <- expr.paad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.paad <-pheatmap(t(expr.genes.paad), show_colnames = FALSE,
                cutree_cols = 3) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.paad$tree_col
assignments.paad <- cutree(clusteredSamples, k=3) # k = cutree_cols
groupAssignments.paad <- data.frame(Group=factor(assignments.paad))
p.paad <- pheatmap(t(expr.genes.paad), show_colnames = FALSE,
                 cutree_cols = 3, annotation = groupAssignments.paad)

### Merging group assignments 
df.expr.paad <- data.frame(expr.genes.paad)
df.expr.paad$Sample <- rownames(df.expr.paad)
groupAssignments.paad$SampleID <- rownames(groupAssignments.paad)
df.merged.paad <- merge(df.expr.paad, groupAssignments.paad,
                      by.x = "Sample", by.y = "SampleID",
                      all.x = FALSE, all.y = FALSE)

df.merged.paad$Group[df.merged.paad$Group =="3"] <-"2" # gp2 = 2,3

# merging clinical data
clinical.paad$PatientID <- rownames(clinical.paad)
df.merged.paad$PatientID <- substr(df.merged.paad$Sample, start = 1, stop = 12)
df.merged.paad$PatientAge <- clinical.paad[match(df.merged.paad$PatientID, clinical.paad$PatientID), "yearstobirth"]
df.merged.paad$PatientAge <- as.numeric(df.merged.paad$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.paad$SubID <- substr(df.merged.paad$Sample, start = 1, stop = 19)
S1_corr_data_PAAD_filter$SubID <- substr(S1_corr_data_PAAD_filter$SampleID, start = 1, stop = 19)
df.merged.paad <- merge(df.merged.paad, S1_corr_data_PAAD_filter,
                      by.x = "SubID", by.y = "SubID", 
                      all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.highPRR5L", "highDEPDC6.lowPRR5L")
df.merged.paad$Status <- values[df.merged.paad$Group]

## MP ##
library(ggpubr)
violin.paad.mp <- ggviolin(df.merged.paad, x = "Status", y = "MP", fill = "Status", palette = c("#f291d0", "#95d108"), add = "boxplot", add.params = list(fill = "white"))
comparison.paad <- list(c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L"))

p.mp <- violin.paad.mp + stat_compare_means(comparisons = comparison.paad) + 
  stat_compare_means(label.y = 0.39, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.paad.mp <- ggpar(p.mp, 
                   main = "PAAD",
                   font.main = c(16, "bold"),
                   xlab = "Gene expression status",
                   ylab = "% S1 (MP)",
                   font.x = c(14, "bold"), 
                   font.y = c(14, "bold"),
                   font.ytickslab = 14,
                   font.xtickslab = c(1, "white"),
                   legend = "none")
mean.mp <- aggregate(df.merged.paad$MP ~ df.merged.paad$Status, df.merged.paad, mean)

## DS ##
violin.paad.ds <- ggviolin(df.merged.paad, x = "Status", y = "DS", fill = "Status", palette = c("#eb57c8", "#2da10a"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.paad.ds + stat_compare_means(comparisons = comparison.paad) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.paad.ds <- ggpar(p.ds, 
                   main = "PAAD",
                   font.main = c(16, "bold"),
                   xlab = "Gene expression status",
                   ylab = "% S1 (DS)",
                   font.x = c(14, "bold"), 
                   font.y = c(14, "bold"),
                   font.ytickslab = 14,
                   font.xtickslab = c(1, "white"),
                   legend = "none")
mean.ds <- aggregate(df.merged.paad$DS ~ df.merged.paad$Status, df.merged.paad, mean)

save(df.merged.paad, file = "df.merged.paad.RData")
