#### EXPRESSION DATA: mTOR PATHWAY GENES - BRCA ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.brca <- log2(expr.brca+1) 
expr.genes.brca <- expr.brca[,genes]

library(pheatmap)
pdf("heatmap.brca.pdf")
p.brca <-pheatmap(t(expr.genes.brca), show_colnames = FALSE,
                  cutree_cols = 4)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.brca$tree_col
assignments.brca <- cutree(clusteredSamples, k=4) # k = cutree_cols
groupAssignments.brca <- data.frame(Group=factor(assignments.brca))
p.brca <- pheatmap(t(expr.genes.brca), show_colnames = FALSE,
                   cutree_cols = 4, annotation = groupAssignments.brca)

### Merging group assignments 
df.expr.brca <- data.frame(expr.genes.brca)
df.expr.brca$Sample <- rownames(df.expr.brca)
groupAssignments.brca$SampleID <- rownames(groupAssignments.brca)
df.merged.brca <- merge(df.expr.brca, groupAssignments.brca,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.brca$Group[df.merged.brca$Group =="4"] <-"3"

# Additional columns  
df.merged.brca$Cancer <- "BRCA"
for (i in 1:1212) df.merged.brca$Score[i] = sum(df.merged.brca[i,2:39])/39
values <- c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L","highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", 
            "lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L")
df.merged.brca$Status <- values[df.merged.brca$Group]

# merging clinical data
clinical.brca$PatientID <- rownames(clinical.brca)
df.merged.brca$PatientID <- substr(df.merged.brca$Sample, start = 1, stop = 12)
df.merged.brca$PatientAge <- clinical.brca[match(df.merged.brca$PatientID, clinical.brca$PatientID), "yearstobirth"]
df.merged.brca$PatientAge <- as.numeric(df.merged.brca$PatientAge)

# merging % S1 data (DS and MP)
df.merged.brca$SubID <- substr(df.merged.brca$Sample, start = 1, stop = 19)
S1_corr_data_BRCA$SubID <- substr(S1_corr_data_BRCA$SampleID, start = 1, stop = 19)
df.merged.brca <- merge(df.merged.brca, S1_corr_data_BRCA,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.brca, file = "merged.file.brca.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.brca <- list(c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"), 
                        c("highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L"), 
                        c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.brca <- ggboxplot(df.merged.brca, x="Status", y="MP", fill = "Status", palette = c("#fa7393", "#ffc861", "#93ff61"), shape = "Status")
p.mp <- boxplot.mp.brca + stat_compare_means(comparisons = comparison.brca) + 
  stat_compare_means(label.y = 0.69, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.brca.mp <- ggpar(p.mp, 
                     main = "BRCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.brca$MP ~ df.merged.brca$Status, df.merged.brca, mean)

### DS ###
boxplot.ds.brca <- ggboxplot(df.merged.brca, x="Status", y="DS", fill = "Status", palette = c("#ff4c77", "#ffae15", "#40c303"), shape = "Status")
p.ds <- boxplot.ds.brca + stat_compare_means(comparisons = comparison.brca) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.brca.ds <- ggpar(p.ds, 
                     main = "BRCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.brca$DS ~ df.merged.brca$Status, df.merged.brca, mean)

cor.test(df.merged.brca$MP, df.merged.brca$Score)
cor.test(df.merged.brca$DS, df.merged.brca$Score)
