#### EXPRESSION DATA: mTOR PATHWAY GENES - STAD ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.stad <- log2(expr.stad+1) 
expr.genes.stad <- expr.stad[,genes]

library(pheatmap)
pdf("heatmap.stad.pdf")
p.stad <-pheatmap(t(expr.genes.stad), show_colnames = FALSE,
                  cutree_cols = 2)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.stad$tree_col
assignments.stad <- cutree(clusteredSamples, k=2) # k = cutree_cols
groupAssignments.stad <- data.frame(Group=factor(assignments.stad))
p.stad <- pheatmap(t(expr.genes.stad), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments.stad)

### Merging group assignments 
df.expr.stad <- data.frame(expr.genes.stad)
df.expr.stad$Sample <- rownames(df.expr.stad)
groupAssignments.stad$SampleID <- rownames(groupAssignments.stad)
df.merged.stad <- merge(df.expr.stad, groupAssignments.stad,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# Additional columns  
df.merged.stad$Cancer <- "STAD"
for (i in 1:450) df.merged.stad$Score[i] = sum(df.merged.stad[i,2:39])/39
values <- c("lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L")
df.merged.stad$Status <- values[df.merged.stad$Group]

# merging clinical data
clinical.stad$PatientID <- rownames(clinical.stad)
df.merged.stad$PatientID <- substr(df.merged.stad$Sample, start = 1, stop = 12)
df.merged.stad$PatientAge <- clinical.stad[match(df.merged.stad$PatientID, clinical.stad$PatientID), "yearstobirth"]
df.merged.stad$PatientAge <- as.numeric(df.merged.stad$PatientAge)

# merging % S1 data (DS and MP)
df.merged.stad$SubID <- substr(df.merged.stad$Sample, start = 1, stop = 19)
S1_corr_data_STAD$SubID <- substr(S1_corr_data_STAD$SampleID, start = 1, stop = 19)
df.merged.stad <- merge(df.merged.stad, S1_corr_data_STAD,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.stad, file = "merged.file.stad.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.stad <- list(c("lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.stad <- ggboxplot(df.merged.stad, x="Status", y="MP", fill = "Status", palette = c("#93ff61", "#fa7393"), shape = "Status")
p.mp <- boxplot.mp.stad + stat_compare_means(comparisons = comparison.stad) + 
  stat_compare_means(label.y = 0.52, size = 5)
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

### DS ###
boxplot.ds.stad <- ggboxplot(df.merged.stad, x="Status", y="DS", fill = "Status", palette = c("#40c303", "#ff4c77"), shape = "Status")
p.ds <- boxplot.ds.stad + stat_compare_means(comparisons = comparison.stad) + 
  stat_compare_means(label.y = 1.3, size = 5)
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

cor.test(df.merged.stad$MP, df.merged.stad$Score)
cor.test(df.merged.stad$DS, df.merged.stad$Score)
