#### EXPRESSION DATA: mTOR PATHWAY GENES - ESCA ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.esca <- log2(expr.esca+1) 
expr.genes.esca <- expr.esca[,genes]

library(pheatmap)
pdf("heatmap.esca.pdf")
p.esca <-pheatmap(t(expr.genes.esca), show_colnames = FALSE,
                  cutree_cols = 3)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.esca$tree_col
assignments.esca <- cutree(clusteredSamples, k=3) # k = cutree_cols
groupAssignments.esca <- data.frame(Group=factor(assignments.esca))
p.esca <- pheatmap(t(expr.genes.esca), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments.esca)

### Merging group assignments 
df.expr.esca <- data.frame(expr.genes.esca)
df.expr.esca$Sample <- rownames(df.expr.esca)
groupAssignments.esca$SampleID <- rownames(groupAssignments.esca)
df.merged.esca <- merge(df.expr.esca, groupAssignments.esca,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# Additional columns  
df.merged.esca$Cancer <- "ESCA"
for (i in 1:196) df.merged.esca$Score[i] = sum(df.merged.esca[i,2:39])/39
values <- c("highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", 
            "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L")
df.merged.esca$Status <- values[df.merged.esca$Group]

# merging clinical data
clinical.esca$PatientID <- rownames(clinical.esca)
df.merged.esca$PatientID <- substr(df.merged.esca$Sample, start = 1, stop = 12)
df.merged.esca$PatientAge <- clinical.esca[match(df.merged.esca$PatientID, clinical.esca$PatientID), "yearstobirth"]
df.merged.esca$PatientAge <- as.numeric(df.merged.esca$PatientAge)

# merging % S1 data (DS and MP)
df.merged.esca$SubID <- substr(df.merged.esca$Sample, start = 1, stop = 19)
S1_corr_data_ESCA_filter$SubID <- substr(S1_corr_data_ESCA_filter$Sample_ID, start = 1, stop = 19)
df.merged.esca <- merge(df.merged.esca, S1_corr_data_ESCA_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.esca, file = "merged.file.esca.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.esca <- list(c("highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L"),
                        c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L"), 
                        c("lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L", "highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.esca <- ggboxplot(df.merged.esca, x="Status", y="MP", fill = "Status", palette = c("#47bfff", "#fa7393", "#b789f0"), shape = "Status")
p.mp <- boxplot.mp.esca + stat_compare_means(comparisons = comparison.esca) + 
  stat_compare_means(label.y = 0.35, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.esca.mp <- ggpar(p.mp, 
                     main = "ESCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.esca$MP ~ df.merged.esca$Status, df.merged.esca, mean)

### DS ###
boxplot.ds.esca <- ggboxplot(df.merged.esca, x="Status", y="DS", fill = "Status", palette = c("#358fff", "#ff4c77", "#9900ff"), shape = "Status")
p.ds <- boxplot.ds.esca + stat_compare_means(comparisons = comparison.esca) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.esca.ds <- ggpar(p.ds, 
                     main = "ESCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

mean.ds <- aggregate(df.merged.esca$DS ~ df.merged.esca$Status, df.merged.esca, mean)

cor.test(df.merged.esca$MP, df.merged.esca$Score)
cor.test(df.merged.esca$DS, df.merged.esca$Score)

