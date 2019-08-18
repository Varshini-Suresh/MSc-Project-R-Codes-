#### EXPRESSION DATA: mTOR PATHWAY GENES - HNSC ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.hnsc <- log2(expr.hnsc+1) 
expr.genes.hnsc <- expr.hnsc[,genes]

library(pheatmap)
pdf("heatmap.hnsc.pdf")
p.hnsc <-pheatmap(t(expr.genes.hnsc), show_colnames = FALSE,
                  cutree_cols = 5)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.hnsc$tree_col
assignments.hnsc <- cutree(clusteredSamples, k=5) # k = cutree_cols
groupAssignments.hnsc <- data.frame(Group=factor(assignments.hnsc))
p.hnsc <- pheatmap(t(expr.genes.hnsc), show_colnames = FALSE,
                   cutree_cols = 5, annotation = groupAssignments.hnsc)

### Merging group assignments 
df.expr.hnsc <- data.frame(expr.genes.hnsc)
df.expr.hnsc$Sample <- rownames(df.expr.hnsc)
groupAssignments.hnsc$SampleID <- rownames(groupAssignments.hnsc)
df.merged.hnsc <- merge(df.expr.hnsc, groupAssignments.hnsc,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.hnsc$Group[df.merged.hnsc$Group == "5"] <- "2" #gp1 = 1,2,4 & gp2 = 3,5

# Additional columns  
df.merged.hnsc$Cancer <- "HNSC"
for (i in 1:566) df.merged.hnsc$Score[i] = sum(df.merged.hnsc[i,2:39])/39
values <- c("lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L")
df.merged.hnsc$Status <- values[df.merged.hnsc$Group]

# merging clinical data
clinical.hnsc$PatientID <- rownames(clinical.hnsc)
df.merged.hnsc$PatientID <- substr(df.merged.hnsc$Sample, start = 1, stop = 12)
df.merged.hnsc$PatientAge <- clinical.hnsc[match(df.merged.hnsc$PatientID, clinical.hnsc$PatientID), "yearstobirth"]
df.merged.hnsc$PatientAge <- as.numeric(df.merged.hnsc$PatientAge)

# merging % S1 data (DS and MP)
df.merged.hnsc$SubID <- substr(df.merged.hnsc$Sample, start = 1, stop = 19)
S1_corr_data_HNSC$SubID <- substr(S1_corr_data_HNSC$SampleID, start = 1, stop = 19)
df.merged.hnsc <- merge(df.merged.hnsc, S1_corr_data_HNSC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.hnsc, file = "merged.file.hnsc.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.hnsc <- list(c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.hnsc <- ggboxplot(df.merged.hnsc, x="Status", y="MP", fill = "Status", palette = c("#b789f0", "#fa7393"), shape = "Status")
p.mp <- boxplot.mp.hnsc + stat_compare_means(comparisons = comparison.hnsc) + 
  stat_compare_means(label.y = 0.38, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.hnsc.mp <- ggpar(p.mp, 
                     main = "HNSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.hnsc$MP ~ df.merged.hnsc$Status, df.merged.hnsc, mean)

### DS ###
boxplot.ds.hnsc <- ggboxplot(df.merged.hnsc, x="Status", y="DS", fill = "Status", palette = c("#9900ff", "#ff4c77"), shape = "Status")
p.ds <- boxplot.ds.hnsc + stat_compare_means(comparisons = comparison.hnsc) + 
  stat_compare_means(label.y = 1.2, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.hnsc.ds <- ggpar(p.ds, 
                     main = "HNSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

mean.ds <- aggregate(df.merged.hnsc$DS ~ df.merged.hnsc$Status, df.merged.hnsc, mean)

cor.test(df.merged.hnsc$MP, df.merged.hnsc$Score)
cor.test(df.merged.hnsc$DS, df.merged.hnsc$Score)

