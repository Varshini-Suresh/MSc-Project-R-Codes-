#### mTORC1 AND mTORC2 EXPRESSION DATA - STAD ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.stad <- log2(expr.stad+1) # log of exp data 
expr.genes.stad <- expr.stad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.stad <-pheatmap(t(expr.genes.stad), show_colnames = FALSE,
                  cutree_cols = 6) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.stad$tree_col
assignments.stad <- cutree(clusteredSamples, k=6) # k = cutree_cols
groupAssignments.stad <- data.frame(Group=factor(assignments.stad))
p.stad <- pheatmap(t(expr.genes.stad), show_colnames = FALSE,
                   cutree_cols = 6, annotation = groupAssignments.stad)

### Merging group assignments 
df.expr.stad <- data.frame(expr.genes.stad)
df.expr.stad$Sample <- rownames(df.expr.stad)
groupAssignments.stad$SampleID <- rownames(groupAssignments.stad)
df.merged.stad <- merge(df.expr.stad, groupAssignments.stad,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.stad$Group[df.merged.stad$Group =="2"] <-"1" #gp1 = 1,5,6; gp2 = 2,3 & gp3 = 4

# merging clinical data
clinical.stad$PatientID <- rownames(clinical.stad)
df.merged.stad$PatientID <- substr(df.merged.stad$Sample, start = 1, stop = 12)
df.merged.stad$PatientAge <- clinical.stad[match(df.merged.stad$PatientID, clinical.stad$PatientID), "yearstobirth"]
df.merged.stad$PatientAge <- as.numeric(df.merged.stad$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.stad$SubID <- substr(df.merged.stad$Sample, start = 1, stop = 19)
S1_corr_data_STAD$SubID <- substr(S1_corr_data_STAD$SampleID, start = 1, stop = 19)
df.merged.stad <- merge(df.merged.stad, S1_corr_data_STAD,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

filter.stad <- subset(df.merged.stad, df.merged.stad$Group != "3")
df.merged.stad.copy <- df.merged.stad
df.merged.stad <- filter.stad

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L")
df.merged.stad$Status <- values[df.merged.stad$Group]

## MP ##
library(ggpubr)
violin.stad.mp <- ggviolin(df.merged.stad, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#f291d0"), add = "boxplot", add.params = list(fill = "white"))
comparison.stad <- list(c("highDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L"))

p.mp <- violin.stad.mp + stat_compare_means(comparisons = comparison.stad) + 
  stat_compare_means(label.y = 0.5, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.stad.mp <- ggpar(p.mp, 
                     main = "STAD",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.stad$MP ~ df.merged.stad$Status, df.merged.stad, mean)

## DS ##
violin.stad.ds <- ggviolin(df.merged.stad, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#eb57c8"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.stad.ds + stat_compare_means(comparisons = comparison.stad) + 
  stat_compare_means(label.y = 1.35, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.stad.ds <- ggpar(p.ds, 
                     main = "STAD", 
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.stad$DS ~ df.merged.stad$Status, df.merged.stad, mean)

save (df.merged.stad, file = "df.merged.stad1.RData")
