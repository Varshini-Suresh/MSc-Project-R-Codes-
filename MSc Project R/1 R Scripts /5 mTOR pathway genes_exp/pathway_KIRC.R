#### EXPRESSION DATA: mTOR PATHWAY GENES - KIRC ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.kirc <- log2(expr.kirc+1) 
expr.genes.kirc <- expr.kirc[,genes]

library(pheatmap)
pdf("heatmap.kirc.pdf")
p.kirc <-pheatmap(t(expr.genes.kirc), show_colnames = FALSE,
                  cutree_cols = 5)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.kirc$tree_col
assignments.kirc <- cutree(clusteredSamples, k=5) # k = cutree_cols
groupAssignments.kirc <- data.frame(Group=factor(assignments.kirc))
p.kirc <- pheatmap(t(expr.genes.kirc), show_colnames = FALSE,
                   cutree_cols = 5, annotation = groupAssignments.kirc)

### Merging group assignments 
df.expr.kirc <- data.frame(expr.genes.kirc)
df.expr.kirc$Sample <- rownames(df.expr.kirc)
groupAssignments.kirc$SampleID <- rownames(groupAssignments.kirc)
df.merged.kirc <- merge(df.expr.kirc, groupAssignments.kirc,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.kirc$Group[df.merged.kirc$Group == "3"] <- "2" #gp2 = 2,3,4,5

# Additional columns  
df.merged.kirc$Cancer <- "KIRC"
for (i in 1:606) df.merged.kirc$Score[i] = sum(df.merged.kirc[i,2:39])/39
values <- c("highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L")
df.merged.kirc$Status <- values[df.merged.kirc$Group]

# merging clinical data
clinical.kirc$PatientID <- rownames(clinical.kirc)
df.merged.kirc$PatientID <- substr(df.merged.kirc$Sample, start = 1, stop = 12)
df.merged.kirc$PatientAge <- clinical.kirc[match(df.merged.kirc$PatientID, clinical.kirc$PatientID), "yearstobirth"]
df.merged.kirc$PatientAge <- as.numeric(df.merged.kirc$PatientAge)

# merging % S1 data (DS and MP)
df.merged.kirc$SubID <- substr(df.merged.kirc$Sample, start = 1, stop = 19)
S1_corr_data_KIRC$SubID <- substr(S1_corr_data_KIRC$SampleID, start = 1, stop = 19)
df.merged.kirc <- merge(df.merged.kirc, S1_corr_data_KIRC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.kirc, file = "merged.file.kirc.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.kirc <- list(c("highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.kirc <- ggboxplot(df.merged.kirc, x="Status", y="MP", fill = "Status", palette = c("#ffc861", "#fa7393"), shape = "Status")
p.mp <- boxplot.mp.kirc + stat_compare_means(comparisons = comparison.kirc) + 
  stat_compare_means(label.y = 0.29, size = 5)
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

### DS ###
boxplot.ds.kirc <- ggboxplot(df.merged.kirc, x="Status", y="DS", fill = "Status", palette = c("#ffae15", "#ff4c77"), shape = "Status")
p.ds <- boxplot.ds.kirc + stat_compare_means(comparisons = comparison.kirc) + 
  stat_compare_means(label.y = 1.1, size = 5)
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

cor.test(df.merged.kirc$MP, df.merged.kirc$Score)
cor.test(df.merged.kirc$DS, df.merged.kirc$Score)
