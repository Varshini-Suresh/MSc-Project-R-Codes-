#### EXPRESSION DATA: mTOR PATHWAY GENES - BLCA ####
library(TCGA2STAT)
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8", "AKT1S1", "RICTOR", "MAPKAP1", "PRR5L", "WNT9B", 
           "LRP5", "DVL1", "GNAQ", "PIK3R1", "IRS1", "PDPK1", "PTEN", "HRAS", "KRAS", "NRAS", "MAPK15", 
           "RPS6KA3", "GSK3B", "DDIT4", "PRKAA2", "PRKCA", "SGK1", "TSC1", "TSC2", "TBC1D7", "FKBP1A", 
           "GRB10", "EIF4G1", "KIAA0652", "RB1CC1", "ULK1", "RPS6KB1", "EIF4EBP1", "EIF4EBP2")

expr.blca <- log2(expr.blca+1) 
expr.genes.blca <- expr.blca[,genes]

library(pheatmap)
pdf("heatmap.blca.pdf")
p.blca <-pheatmap(t(expr.genes.blca), show_colnames = FALSE,
                  cutree_cols = 1, fontsize = 16)
dev.off()

# Extracting cluster assignments for each sample
clusteredSamples <- p.blca$tree_col
assignments.blca <- cutree(clusteredSamples, k=3) # k = cutree_cols
groupAssignments.blca <- data.frame(Group=factor(assignments.blca))
p.blca <- pheatmap(t(expr.genes.blca), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments.blca)

### Merging group assignments 
df.expr.blca <- data.frame(expr.genes.blca)
df.expr.blca$Sample <- rownames(df.expr.blca)
groupAssignments.blca$SampleID <- rownames(groupAssignments.blca)
df.merged.blca <- merge(df.expr.blca, groupAssignments.blca,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# Additional columns  
df.merged.blca$Cancer <- "BLCA"
for (i in 1:427) df.merged.blca$Score[i] = sum(df.merged.blca[i,2:39])/39
values <- c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L",
            "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L")
df.merged.blca$Status <- values[df.merged.blca$Group]

# merging clinical data
clinical.blca$PatientID <- rownames(clinical.blca)
df.merged.blca$PatientID <- substr(df.merged.blca$Sample, start = 1, stop = 12)
df.merged.blca$PatientAge <- clinical.blca[match(df.merged.blca$PatientID, clinical.blca$PatientID), "yearstobirth"]
df.merged.blca$PatientAge <- as.numeric(df.merged.blca$PatientAge)

# merging % S1 data (DS and MP)
df.merged.blca$SubID <- substr(df.merged.blca$Sample, start = 1, stop = 19)
S1_corr_data_BLCA$SubID <- substr(S1_corr_data_BLCA$SampleID, start = 1, stop = 19)
df.merged.blca <- merge(df.merged.blca, S1_corr_data_BLCA,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

save(df.merged.blca, file = "merged.file.blca.RData")

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.blca <- list(c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L"), 
                        c("highMAPK15.lowPRKAA2.highDEPDC6.lowPRR5L", "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L"), 
                        c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L"))

### MP ###
library(ggpubr)
boxplot.mp.blca <- ggboxplot(df.merged.blca, x="Status", y="MP", fill = "Status", palette = c("#fa7393", "#47bfff", "#b789f0"), shape = "Status")
p.mp <- boxplot.mp.blca + stat_compare_means(comparisons = comparison.blca) + 
  stat_compare_means(label.y = 0.3, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.blca.mp <- ggpar(p.mp, 
                     main = "BLCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.blca$MP ~ df.merged.blca$Status, df.merged.blca, mean)

### DS ###
boxplot.ds.blca <- ggboxplot(df.merged.blca, x="Status", y="DS", fill = "Status", palette = c("#ff4c77", "#358fff", "#9459de"), shape = "Status")
p.ds <- boxplot.ds.blca + stat_compare_means(comparisons = comparison.blca) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.blca.ds <- ggpar(p.ds, 
                     main = "BLCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

mean.ds <- aggregate(df.merged.blca$DS ~ df.merged.blca$Status, df.merged.blca, mean)

cor.test(df.merged.blca$MP, df.merged.blca$Score)
cor.test(df.merged.blca$DS, df.merged.blca$Score)

#### PATIENT OVERALL SURVIVAL (%) ####
# Addition of a few more columns from clinical data
df.merged.blca$VitalStatus <- clinical.blca[match(df.merged.blca$PatientID, clinical.blca$PatientID), "vitalstatus"]
  df.merged.blca$VitalStatus<- as.numeric(as.character(df.merged.blca$VitalStatus))
df.merged.blca$DaystoDeath <- clinical.blca[match(df.merged.blca$PatientID, clinical.blca$PatientID), "daystodeath"]
  df.merged.blca$DaystoDeath<- as.numeric(as.character(df.merged.blca$DaystoDeath))
df.merged.blca$Followup <- clinical.blca[match(df.merged.blca$PatientID, clinical.blca$PatientID), "daystolastfollowup"]
  df.merged.blca$Followup<- as.numeric(as.character(df.merged.blca$Followup))

## Calculating overall survival 
df.merged.blca$OS <- NA
for (i in 1:nrow(df.merged.blca)) {
  if (!is.na(df.merged.blca[i,]$DaystoDeath)) {
    df.merged.blca[i,]$OS <- df.merged.blca[i,]$DaystoDeath
  } else {
    if (!is.na(df.merged.blca[i,]$Followup)) {
      df.merged.blca[i,]$OS <- df.merged.blca[i,]$Followup
    }}}

### Plotting survival 
require("survival")
library(survminer)
fit.blca <- survfit(Surv(OS, VitalStatus) ~ Group, data = df.merged.blca)

ggsurvplot(fit.blca, 
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = "strata", 
           ggtheme = theme_bw())

save (df.merged.blca, file = "merged.file.blca.RData")
