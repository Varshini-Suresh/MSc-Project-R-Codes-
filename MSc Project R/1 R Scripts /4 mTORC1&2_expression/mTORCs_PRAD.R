#### mTORC1 AND mTORC2 EXPRESSION DATA - PRAD ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.prad <- log2(expr.prad+1) # log of exp data 
expr.genes.prad <- expr.prad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.prad <-pheatmap(t(expr.genes.prad), show_colnames = FALSE,
                cutree_cols = 1) # to be specified after visualising first 


#################

# Extracting cluster assignments for each sample
clusteredSamples <- p.prad$tree_col
assignments.prad <- cutree(clusteredSamples, k=1) # k = cutree_cols
groupAssignments.prad <- data.frame(Group=factor(assignments.prad))
p.prad <- pheatmap(t(expr.genes.prad), show_colnames = FALSE,
                 cutree_cols = 1, annotation = groupAssignments.prad)

### Merging group assignments 
df.expr.prad <- data.frame(expr.genes.prad)
df.expr.prad$Sample <- rownames(df.expr.prad)
groupAssignments.prad$SampleID <- rownames(groupAssignments.prad)
df.merged.prad <- merge(df.expr.prad, groupAssignments.prad,
                      by.x = "Sample", by.y = "SampleID",
                      all.x = FALSE, all.y = FALSE)

df.merged.prad$Group[df.merged.prad$Group =="2"] <-"1" # gp1 = 1,2 and gp2 = 3

# merging clinical data
clinical.prad$PatientID <- rownames(clinical.prad)
df.merged.prad$PatientID <- substr(df.merged.prad$Sample, start = 1, stop = 12)
df.merged.prad$PatientAge <- clinical.prad[match(df.merged.prad$PatientID, clinical.prad$PatientID), "yearstobirth"]
df.merged.prad$PatientAge <- as.numeric(df.merged.prad$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.prad$SubID <- substr(df.merged.prad$Sample, start = 1, stop = 19)
S1_corr_data_PRAD_filter$SubID <- substr(S1_corr_data_PRAD_filter$SampleID, start = 1, stop = 19)
df.merged.prad <- merge(df.merged.prad, S1_corr_data_PRAD_filter,
                      by.x = "SubID", by.y = "SubID", 
                      all.x = FALSE, all.y = FALSE)

df.merged.prad$Status <- "highDEPDC6.lowPRR5L"
save(df.merged.prad, file = "df.merged.prad.RData")

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L")
df.merged.prad$Status <- values[df.merged.prad$Group]

## MP ##
library(ggpubr)
violin.prad.mp <- ggviolin(df.merged.prad, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#5dd4f5"), add = "boxplot", add.params = list(fill = "white"))
comparison.prad <- list(c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L"))

p.mp <- violin.prad.mp + stat_compare_means(comparisons = comparison.prad) + 
  stat_compare_means(label.y = 0.20, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.prad.mp <- ggpar(p.mp, 
                   main = "PRAD",
                   font.main = c(16, "bold"),
                   xlab = "Gene expression status",
                   ylab = "% S1 (MP)",
                   font.x = c(14, "bold"), 
                   font.y = c(14, "bold"),
                   font.ytickslab = 14,
                   font.xtickslab = c(1, "white"),
                   legend = "none")
mean.mp <- aggregate(df.merged.prad$MP ~ df.merged.prad$Status, df.merged.prad, mean)

## DS ##
violin.prad.ds <- ggviolin(df.merged.prad, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#1155f5"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.prad.ds + stat_compare_means(comparisons = comparison.prad) + 
  stat_compare_means(label.y = 1.0, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.prad.ds <- ggpar(p.ds, 
                   main = "PRAD",
                   font.main = c(16, "bold"),
                   xlab = "Gene expression status",
                   ylab = "% S1 (DS)",
                   font.x = c(14, "bold"), 
                   font.y = c(14, "bold"),
                   font.ytickslab = 14,
                   font.xtickslab = c(1, "white"),
                   legend = "none")
mean.ds <- aggregate(df.merged.prad$DS ~ df.merged.prad$Status, df.merged.prad, mean)
