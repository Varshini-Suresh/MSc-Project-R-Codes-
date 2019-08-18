#### EXPRESSION DATA: mTOR PATHWAY GENES - PAAD ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.paad <- log2(expr.paad+1) 
expr.genes.paad <- expr.paad[,genes]

library(pheatmap)
pdf("heatmap.paad.pdf")
p.paad <-pheatmap(t(expr.genes.paad), show_colnames = FALSE,
                cutree_cols = 2)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.paad$tree_col
assignments.paad <- cutree(clusteredSamples, k=2) # k = cutree_cols
groupAssignments.paad <- data.frame(Group=factor(assignments.paad))
p.paad <- pheatmap(t(expr.genes.paad), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments.paad)

### Merging group assignments 
df.expr.paad <- data.frame(expr.genes.paad)
df.expr.paad$Sample <- rownames(df.expr.paad)
groupAssignments.paad$SampleID <- rownames(groupAssignments.paad)
df.merged.paad <- merge(df.expr.paad, groupAssignments.paad,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# Additional columns  
df.merged.paad$Cancer <- "PAAD"
for (i in 1:183) df.merged.paad$Score[i] = sum(df.merged.paad[i,2:39])/39
values <- c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L")
df.merged.paad$Status <- values[df.merged.paad$Group]

# merging clinical data
clinical.paad$PatientID <- rownames(clinical.paad)
df.merged.paad$PatientID <- substr(df.merged.paad$Sample, start = 1, stop = 12)
df.merged.paad$PatientAge <- clinical.paad[match(df.merged.paad$PatientID, clinical.paad$PatientID), "yearstobirth"]
df.merged.paad$PatientAge <- as.numeric(df.merged.paad$PatientAge)

# merging % S1 data (DS and MP)
df.merged.paad$SubID <- substr(df.merged.paad$Sample, start = 1, stop = 19)
S1_corr_data_PAAD_filter$SubID <- substr(S1_corr_data_PAAD_filter$SampleID, start = 1, stop = 19)
df.merged.paad <- merge(df.merged.paad, S1_corr_data_PAAD_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.paad, file = "merged.file.paad.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.paad <- list(c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.paad <- ggboxplot(df.merged.paad, x="Status", y="MP", fill = "Status", palette = c("#ffc861", "#fa7393"), shape = "Status")
p.mp <- boxplot.mp.paad + stat_compare_means(comparisons = comparison.paad) + 
  stat_compare_means(label.y = 0.35, size = 5)
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

### DS ###
boxplot.ds.paad <- ggboxplot(df.merged.paad, x="Status", y="DS", fill = "Status", palette = c("#ffae15", "#ff4c77"), shape = "Status")
p.ds <- boxplot.ds.paad + stat_compare_means(comparisons = comparison.paad) + 
  stat_compare_means(label.y = 1.2, size = 5)
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

cor.test(df.merged.paad$MP, df.merged.paad$Score)
cor.test(df.merged.paad$DS, df.merged.paad$Score)
