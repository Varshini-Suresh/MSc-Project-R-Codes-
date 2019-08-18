#### EXPRESSION DATA: mTOR PATHWAY GENES - PRAD ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.prad <- log2(expr.prad+1) 
expr.genes.prad <- expr.prad[,genes]

library(pheatmap)
pdf("heatmap.prad.pdf")
p.prad <-pheatmap(t(expr.genes.prad), show_colnames = FALSE,
                cutree_cols = 1)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.prad$tree_col
assignments.prad <- cutree(clusteredSamples, k=1) # k = cutree_cols
groupAssignments.prad <- data.frame(Group=factor(assignments.prad))
p.prad <- pheatmap(t(expr.genes.prad), show_colnames = FALSE,
                 cutree_cols = 1, annotation = groupAssignments.prad)

### Merging group assignments 
df.expr.prad <- data.frame(expr.genes.prad)
df.expr.prad$Sample <- rownames(df.expr.prad)
groupAssignments.prad$SampleID <- rownames(groupAssignments.prad)
df.merged.prad <- merge(df.expr.prad, groupAssignments.prad,
                      by.x = "Sample", by.y = "SampleID",
                      all.x = FALSE, all.y = FALSE)

# Additional columns  
df.merged.prad$Cancer <- "PRAD"
for (i in 1:550) df.merged.prad$Score[i] = sum(df.merged.prad[i,2:39])/39
df.merged.prad$Status <- "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L"

# merging clinical data
clinical.prad$PatientID <- rownames(clinical.prad)
df.merged.prad$PatientID <- substr(df.merged.prad$Sample, start = 1, stop = 12)
df.merged.prad$PatientAge <- clinical.prad[match(df.merged.prad$PatientID, clinical.prad$PatientID), "yearstobirth"]
df.merged.prad$PatientAge <- as.numeric(df.merged.prad$PatientAge)

# merging % S1 data (DS and MP)
df.merged.prad$SubID <- substr(df.merged.prad$Sample, start = 1, stop = 19)
S1_corr_data_PRAD$SubID <- substr(S1_corr_data_PRAD$SampleID, start = 1, stop = 19)
df.merged.prad <- merge(df.merged.prad, S1_corr_data_PRAD,
                      by.x = "SubID", by.y = "SubID", 
                      all.x = FALSE, all.y = FALSE)

save(df.merged.prad, file = "merged.file.prad.RData")
