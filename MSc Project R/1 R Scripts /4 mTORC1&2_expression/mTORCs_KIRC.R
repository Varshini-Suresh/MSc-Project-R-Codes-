#### mTORC1 AND mTORC2 EXPRESSION DATA - KIRC ####
# load expr. and clinical data files 
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L")
expr.kirc <- log2(expr.kirc+1) # log of exp data 
expr.genes.kirc <- expr.kirc[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p.kirc <-pheatmap(t(expr.genes.kirc), show_colnames = FALSE,
                  cutree_cols = 11) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples <- p.kirc$tree_col
assignments.kirc <- cutree(clusteredSamples, k=11) # k = cutree_cols
groupAssignments.kirc <- data.frame(Group=factor(assignments.kirc))
p.kirc <- pheatmap(t(expr.genes.kirc), show_colnames = FALSE,
                   cutree_cols = 11, annotation = groupAssignments.kirc)

### Merging group assignments 
df.expr.kirc <- data.frame(expr.genes.kirc)
df.expr.kirc$Sample <- rownames(df.expr.kirc)
groupAssignments.kirc$SampleID <- rownames(groupAssignments.kirc)
df.merged.kirc <- merge(df.expr.kirc, groupAssignments.kirc,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)


df.merged.kirc.copy <- df.merged.kirc
filter.kirc <- subset(df.merged.kirc, df.merged.kirc$Group != "8" & df.merged.kirc$Group != "6" &
                      df.merged.kirc$Group != "9" & df.merged.kirc$Group != "10" & df.merged.kirc$Group != "11")

filter.kirc$Group[filter.kirc$Group =="4"] <-"2"

# merging clinical data
clinical.kirc$PatientID <- rownames(clinical.kirc)
df.merged.kirc$PatientID <- substr(df.merged.kirc$Sample, start = 1, stop = 12)
df.merged.kirc$PatientAge <- clinical.kirc[match(df.merged.kirc$PatientID, clinical.kirc$PatientID), "yearstobirth"]
df.merged.kirc$PatientAge <- as.numeric(df.merged.kirc$PatientAge)

# merging ageing signature data (DS and MP)
df.merged.kirc$SubID <- substr(df.merged.kirc$Sample, start = 1, stop = 19)
S1_corr_data_KIRC$SubID <- substr(S1_corr_data_KIRC$SampleID, start = 1, stop = 19)
df.merged.kirc <- merge(df.merged.kirc, S1_corr_data_KIRC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L")
df.merged.kirc$Status <- values[df.merged.kirc$Group]

## MP ##
library(ggpubr)
violin.kirc.mp <- ggviolin(df.merged.kirc, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#f291d0"), add = "boxplot", add.params = list(fill = "white"))
comparison.kirc <- list(c("highDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L"))

p.mp <- violin.kirc.mp + stat_compare_means(comparisons = comparison.kirc) + 
  stat_compare_means(label.y = 0.35, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.kirc.mp <- ggpar(p.mp, 
                     main = "KIRC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.kirc$MP ~ df.merged.kirc$Status, df.merged.kirc, mean)

## DS ##
violin.kirc.ds <- ggviolin(df.merged.kirc, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#eb57c8"), add = "boxplot", add.params = list(fill = "white"))

p.ds <- violin.kirc.ds + stat_compare_means(comparisons = comparison.kirc) + 
  stat_compare_means(label.y = 1.2, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.kirc.ds <- ggpar(p.ds, 
                     main = "KIRC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.kirc$DS ~ df.merged.kirc$Status, df.merged.kirc, mean)

save(df.merged.kirc, file = "df.merged.kirc.RData")
