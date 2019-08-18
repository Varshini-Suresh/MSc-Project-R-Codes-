#### EXPRESSION DATA FOR MTORC1 IN OV ####
library(TCGA2STAT)
rnaseq.ov <- getTCGA(disease="OV", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.ov <- t(rnaseq.ov$dat)
save(expr.ov, file="expr.ov.RData")
clinical.ov <- rnaseq.ov$clinical
clinical.ov <- data.frame(clinical.ov)
save(clinical.ov, file="clinical.ov.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.ov) to check whether they're present 
expr.ov <- log2(expr.ov+1) # log of exp data 
expr.genes.ov <- expr.ov[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_OV <-pheatmap(t(expr.genes.ov), show_colnames = FALSE,
                 cutree_cols = 3) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_OV <- p_OV$tree_col
assignments_OV <- cutree(clusteredSamples_OV, k=3) # k = cutree_cols
groupAssignments_OV <- data.frame(Group=factor(assignments_OV))
p_OV <- pheatmap(t(expr.genes.ov), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments_OV)

## Merge expr and group assignments 
df.expr.ov <- data.frame(expr.genes.ov)
df.expr.ov$Sample <- rownames(df.expr.ov)
groupAssignments_OV$SampleID <- rownames(groupAssignments_OV)

df.merged.OV <- merge(df.expr.ov, groupAssignments_OV,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.OV$SubID <- substr(df.merged.OV$Sample, start = 1, stop = 19)
S1_corr_data_OV_filter$SubID <- substr(S1_corr_data_OV_filter$SampleID, start = 1, stop = 19)
df.merged.OV <- merge(df.merged.OV, S1_corr_data_OV_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.ov$PatientID <- rownames(clinical.ov)

df.merged.OV$PatientID <- substr(df.merged.OV$Sample, start = 1, stop = 12)
df.merged.OV$PatientAge <- clinical.ov[match(df.merged.OV$PatientID, clinical.ov$PatientID), "yearstobirth"]
df.merged.OV$PatientAge <- as.numeric(df.merged.OV$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.OV$PatientAge, df.merged.OV$`% S1 (DS)`)
cor.test(df.merged.OV$PatientAge, df.merged.OV$`% S1 (MP)`)

save(df.merged.OV, file = "merged.mTORC1.OV.RData")

#### S1 in mTORC1 ACTIVITY ####
# creating a column containing gene expression status
values <- c("highDEPDC6", "mediumDEPDC6", "lowDEPDC6") # refer to heatmap for order
df.merged.OV$Status <- values[df.merged.OV$Group]
subset_OV <- subset(df.merged.OV, df.merged.OV$Group != "3")

## for % S1 MP ##
library(ggpubr)
boxplot_MP_OV <- ggboxplot(subset_OV, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_OV <- list(c("highDEPDC6", "mediumDEPDC6"))
p_MP <- boxplot_MP_OV + stat_compare_means(comparisons = comparison_OV) + 
  stat_compare_means(label.y = 0.21, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_OV <- ggpar(p_MP, 
                       main = "OV",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.OV$MP ~ df.merged.OV$Status, df.merged.OV, mean)

## for % S1 DS ##
boxplot_DS_OV <- ggboxplot(subset_OV, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_OV + stat_compare_means(comparisons = comparison_OV) + 
  stat_compare_means(label.y = 0.97, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_OV <- ggpar(p_DS, 
                       main = "OV",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")
mean_DS <- aggregate(df.merged.OV$DS ~ df.merged.OV$Status, df.merged.OV, mean)

save(df.merged.OV, file = "merged.mTORC1.OV.RData")
