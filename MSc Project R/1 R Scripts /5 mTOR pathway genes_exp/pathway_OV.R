#### EXPRESSION DATA: mTOR PATHWAY GENES - OV ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.ov <- log2(expr.ov+1) 
expr.genes.ov <- expr.ov[,genes]

library(pheatmap)
pdf("heatmap.ov.pdf")
p.ov <-pheatmap(t(expr.genes.ov), show_colnames = FALSE,
                  cutree_cols = 1)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.ov$tree_col
assignments.ov <- cutree(clusteredSamples, k=1) # k = cutree_cols
groupAssignments.ov <- data.frame(Group=factor(assignments.ov))
p.ov <- pheatmap(t(expr.genes.ov), show_colnames = FALSE,
                   cutree_cols = 1, annotation = groupAssignments.ov)

### Merging group assignments 
df.expr.ov <- data.frame(expr.genes.ov)
df.expr.ov$Sample <- rownames(df.expr.ov)
groupAssignments.ov$SampleID <- rownames(groupAssignments.ov)
df.merged.ov <- merge(df.expr.ov, groupAssignments.ov,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# Additional columns  
df.merged.ov$Cancer <- "OV"
for (i in 1:307) df.merged.ov$Score[i] = sum(df.merged.ov[i,2:39])/39
df.merged.ov$Status <- "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"

# merging clinical data
clinical.ov$PatientID <- rownames(clinical.ov)
df.merged.ov$PatientID <- substr(df.merged.ov$Sample, start = 1, stop = 12)
df.merged.ov$PatientAge <- clinical.ov[match(df.merged.ov$PatientID, clinical.ov$PatientID), "yearstobirth"]
df.merged.ov$PatientAge <- as.numeric(df.merged.ov$PatientAge)

# merging % S1 data (DS and MP)
df.merged.ov$SubID <- substr(df.merged.ov$Sample, start = 1, stop = 19)
S1_corr_data_OV_filter$SubID <- substr(S1_corr_data_OV_filter$SampleID, start = 1, stop = 19)
df.merged.ov <- merge(df.merged.ov, S1_corr_data_OV_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.ov, file = "merged.file.ov.RData")
