#### EXPRESSION DATA: mTOR PATHWAY GENES - LUAD ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.luad <- log2(expr.luad+1) 
expr.genes.luad <- expr.luad[,genes]

library(pheatmap)
pdf("heatmap.luad.pdf")
p.luad <-pheatmap(t(expr.genes.luad), show_colnames = FALSE,
                  cutree_cols = 5)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.luad$tree_col
assignments.luad <- cutree(clusteredSamples, k=5) # k = cutree_cols
groupAssignments.luad <- data.frame(Group=factor(assignments.luad))
p.luad <- pheatmap(t(expr.genes.luad), show_colnames = FALSE,
                   cutree_cols = 5, annotation = groupAssignments.luad)

### Merging group assignments 
df.expr.luad <- data.frame(expr.genes.luad)
df.expr.luad$Sample <- rownames(df.expr.luad)
groupAssignments.luad$SampleID <- rownames(groupAssignments.luad)
df.merged.luad <- merge(df.expr.luad, groupAssignments.luad,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

df.merged.luad$Group[df.merged.luad$Group =="4"] <-"3" #gp2 = 2,3,5; gp4 = 3

# Additional columns  
df.merged.luad$Cancer <- "LUAD"
for (i in 1:576) df.merged.luad$Score[i] = sum(df.merged.luad[i,2:39])/39
values <- c("highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L",
            "highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L")
df.merged.luad$Status <- values[df.merged.luad$Group]

# merging clinical data
clinical.luad$PatientID <- rownames(clinical.luad)
df.merged.luad$PatientID <- substr(df.merged.luad$Sample, start = 1, stop = 12)
df.merged.luad$PatientAge <- clinical.luad[match(df.merged.luad$PatientID, clinical.luad$PatientID), "yearstobirth"]
df.merged.luad$PatientAge <- as.numeric(df.merged.luad$PatientAge)

# merging % S1 data (DS and MP)
df.merged.luad$SubID <- substr(df.merged.luad$Sample, start = 1, stop = 19)
S1_corr_data_LUAD_filter$SubID <- substr(S1_corr_data_LUAD_filter$SampleID, start = 1, stop = 19)
df.merged.luad <- merge(df.merged.luad, S1_corr_data_LUAD_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.luad, file = "merged.file.luad.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.luad <- list(c("highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L"),
                        c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L"),
                        c("highMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L")) 

### MP ###
library(ggpubr)
boxplot.mp.luad <- ggboxplot(df.merged.luad, x="Status", y="MP", fill = "Status", palette = c("#ffc861", "#fa7393", "#47bfff"), shape = "Status")
p.mp <- boxplot.mp.luad + stat_compare_means(comparisons = comparison.luad) + 
  stat_compare_means(label.y = 0.4, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.luad.mp <- ggpar(p.mp, 
                     main = "LUAD",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.luad$MP ~ df.merged.luad$Status, df.merged.luad, mean)

### DS ###
boxplot.ds.luad <- ggboxplot(df.merged.luad, x="Status", y="DS", fill = "Status", palette = c("#ffae15", "#ff4c77", "#358fff"), shape = "Status")
p.ds <- boxplot.ds.luad + stat_compare_means(comparisons = comparison.luad) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.luad.ds <- ggpar(p.ds, 
                     main = "LUAD",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.ds <- aggregate(df.merged.luad$DS ~ df.merged.luad$Status, df.merged.luad, mean)

cor.test(df.merged.luad$MP, df.merged.luad$Score)
cor.test(df.merged.luad$DS, df.merged.luad$Score)
