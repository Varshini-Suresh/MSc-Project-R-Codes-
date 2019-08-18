#### EXPRESSION DATA FOR MTORC1 IN BRCA ####
library(TCGA2STAT)
rnaseq.brca <- getTCGA(disease="BRCA", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.brca <- t(rnaseq.brca$dat)
save(expr.brca, file="expr.brca.RData")
clinical.brca <- rnaseq.brca$clinical
clinical.brca <- data.frame(clinical.brca)
save(clinical.brca, file="clinical.brca.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.brca) to check whether they're present 
expr.brca <- log2(expr.brca+1) # log of exp data 
expr.genes.brca <- expr.brca[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_BRCA <-pheatmap(t(expr.genes.brca), show_colnames = FALSE,
                  cutree_cols = 3) # to be specified after visualising first

# Extracting cluster assignments for each sample
clusteredSamples_BRCA <- p_BRCA$tree_col
assignments_BRCA <- cutree(clusteredSamples_BRCA, k=3) # k = cutree_cols
groupAssignments_BRCA <- data.frame(Group=factor(assignments_BRCA))
p_BRCA <- pheatmap(t(expr.genes.brca), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments_BRCA)

### Merge expr and group assignments 
df.expr.BRCA <- data.frame(expr.genes.brca)
df.expr.BRCA$Sample <- rownames(df.expr.BRCA)
groupAssignments_BRCA$SampleID <- rownames(groupAssignments_BRCA)

df.merged.BRCA <- merge(df.expr.BRCA, groupAssignments_BRCA,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

# Merging % S1 DS and MP estimates
df.merged.BRCA$SubID <- substr(df.merged.BRCA$Sample, start = 1, stop = 19)
S1_corr_data_BRCA$SubID <- substr(S1_corr_data_BRCA$SampleID, start = 1, stop = 19)

df.merged.BRCA <- merge(df.merged.BRCA, S1_corr_data_BRCA,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.brca$PatientID <- rownames(clinical.brca)
df.merged.BRCA$PatientID <- substr(df.merged.BRCA$Sample, start = 1, stop = 12)
df.merged.BRCA$PatientAge <- clinical.brca[match(df.merged.BRCA$PatientID, clinical.brca$PatientID), "yearstobirth"]
df.merged.BRCA$PatientAge <- as.numeric(df.merged.BRCA$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.BRCA$PatientAge, df.merged.BRCA$`% S1 (MP)`)
cor.test(df.merged.BRCA$PatientAge, df.merged.BRCA$`% S1 (DS)`)

save(df.merged.BRCA, file = "merged.mTORC1.BRCA.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status

values <- c("lowDEPDC6", "highDEPDC6", "veryhighDEPDC6") # refer to heatmap for order
df.merged.BRCA$Status <- values[df.merged.BRCA$Group]
df.merged.BRCA.filter <- subset(df.merged.BRCA, df.merged.BRCA$Group != "3")

## for % S1 MP ##
library(ggpubr)
boxplot_MP_BRCA <- ggboxplot(df.merged.BRCA.filter, x="Status", y="MP", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
comparison_BRCA <- list(c("lowDEPDC6", "highDEPDC6"))
p_MP <- boxplot_MP_BRCA + stat_compare_means(comparisons = comparison_BRCA) + 
  stat_compare_means(label.y = 0.6, size = 6) 
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_BRCA <- ggpar(p_MP, 
                       main = "BRCA",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.BRCA.filter$MP ~ df.merged.BRCA.filter$Status, df.merged.BRCA.filter, mean)

## for % S1 DS ##
boxplot_DS_BRCA <- ggboxplot(df.merged.BRCA.filter, x="Status", y="DS", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
p_DS <- boxplot_DS_BRCA + stat_compare_means(comparisons = comparison_BRCA) + 
    stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_BRCA <- ggpar(p_DS, 
                       main = "BRCA",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.BRCA.filter$DS ~ df.merged.BRCA.filter$Status, df.merged.BRCA.filter, mean)

save(df.merged.BRCA, file = "merged.mTORC1.BRCA.RData")
save(df.merged.BRCA.filter, file = "merged.mTORC1.BRCA.filter.RData")
