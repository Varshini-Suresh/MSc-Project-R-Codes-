#### EXPRESSION DATA: mTOR PATHWAY GENES - LUSC ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.lusc <- log2(expr.lusc+1) 
expr.genes.lusc <- expr.lusc[,genes]

library(pheatmap)
pdf("heatmap.lusc.pdf")
p.lusc <-pheatmap(t(expr.genes.lusc), show_colnames = FALSE,
                  cutree_cols = 7)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.lusc$tree_col
assignments.lusc <- cutree(clusteredSamples, k=7) # k = cutree_cols
groupAssignments.lusc <- data.frame(Group=factor(assignments.lusc))
p.lusc <- pheatmap(t(expr.genes.lusc), show_colnames = FALSE,
                   cutree_cols = 7, annotation = groupAssignments.lusc)

### Merging group assignments 
df.expr.lusc <- data.frame(expr.genes.lusc)
df.expr.lusc$Sample <- rownames(df.expr.lusc)
groupAssignments.lusc$SampleID <- rownames(groupAssignments.lusc)
df.merged.lusc <- merge(df.expr.lusc, groupAssignments.lusc,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.lusc$Group[df.merged.lusc$Group =="6"] <-"3" #gp 1 = 1,7; gp2 = 2,5; gp3 = 3,4,6

# Additional columns  
df.merged.lusc$Cancer <- "LUSC"
for (i in 1:552) df.merged.lusc$Score[i] = sum(df.merged.lusc[i,2:39])/39
values <- c("lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L",
            "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L")
df.merged.lusc$Status <- values[df.merged.lusc$Group]

# merging clinical data
clinical.lusc$PatientID <- rownames(clinical.lusc)
df.merged.lusc$PatientID <- substr(df.merged.lusc$Sample, start = 1, stop = 12)
df.merged.lusc$PatientAge <- clinical.lusc[match(df.merged.lusc$PatientID, clinical.lusc$PatientID), "yearstobirth"]
df.merged.lusc$PatientAge <- as.numeric(df.merged.lusc$PatientAge)

# merging % S1 data (DS and MP)
df.merged.lusc$SubID <- substr(df.merged.lusc$Sample, start = 1, stop = 19)
S1_corr_data_LUSC$SubID <- substr(S1_corr_data_LUSC$SampleID, start = 1, stop = 19)
df.merged.lusc <- merge(df.merged.lusc, S1_corr_data_LUSC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.lusc, file = "merged.file.lusc.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.lusc <- list(c("lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L"),
                        c("lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"), 
                        c("lowMAPK15.lowPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.lusc <- ggboxplot(df.merged.lusc, x="Status", y="MP", fill = "Status", palette = c("#ffc861", "#93ff61", "#b789f0"), shape = "Status")
p.mp <- boxplot.mp.lusc + stat_compare_means(comparisons = comparison.lusc) + 
  stat_compare_means(label.y = 0.34, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.lusc.mp <- ggpar(p.mp, 
                     main = "LUSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.lusc$MP ~ df.merged.lusc$Status, df.merged.lusc, mean)

### DS ###
boxplot.ds.lusc <- ggboxplot(df.merged.lusc, x="Status", y="DS", fill = "Status", palette = c("#ffae15", "#40c303", "#9900ff"), shape = "Status")
p.ds <- boxplot.ds.lusc + stat_compare_means(comparisons = comparison.lusc) + 
  stat_compare_means(label.y = 1.2, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.lusc.ds <- ggpar(p.ds, 
                     main = "LUSC",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.lusc$DS ~ df.merged.lusc$Status, df.merged.lusc, mean)

cor.test(df.merged.lusc$MP, df.merged.lusc$Score)
cor.test(df.merged.lusc$DS, df.merged.lusc$Score)
