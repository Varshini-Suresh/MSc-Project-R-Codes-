#### EXPRESSION DATA FOR MTORC1 IN BLCA ####
library(TCGA2STAT)
rnaseq.blca <- getTCGA(disease="BLCA", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
# extract relevant matrix:
expr.blca <- t(rnaseq.blca$dat)
save(expr.blca, file="expr.blca.RData")
clinical.blca <- rnaseq.blca$clinical
save(clinical.blca, file="clinical.blca.RData")

### Specify our list of genes of interest:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8")
# Use "MTOR" %in% colnames(expr.blca) to check whether 
expr.blca <- log2(expr.blca+1) # log of exp data 
expr.genes.blca <- expr.blca[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_BLCA<-pheatmap(t(expr.genes.blca), show_colnames = FALSE,
                 cutree_cols = 2) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples_BLCA <- p_BLCA$tree_col
assignments_BLCA <- cutree(clusteredSamples_BLCA, k=2) # k = cutree_cols
groupAssignments_BLCA <- data.frame(Group=factor(assignments_BLCA))
p_BLCA <- pheatmap(t(expr.genes.blca), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments_BLCA)

## Merge expr and group assignments 
df.expr.BLCA <- data.frame(expr.genes.blca)
df.expr.BLCA$Sample <- rownames(df.expr.BLCA)
groupAssignments_BLCA$SampleID <- rownames(groupAssignments_BLCA)

df.merged.BLCA <- merge(df.expr.BLCA, groupAssignments_BLCA,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

## Merging % S1 DS and MP estimates
df.merged.BLCA$SubID <- substr(df.merged.BLCA$Sample, start = 1, stop = 19)
# < LOAD S1.Corr.data here >
S1_corr_data_BLCA$SubID <- substr(S1_corr_data_BLCA$SampleID, start = 1, stop = 19)
df.merged.BLCA <- merge(df.merged.BLCA, S1_corr_data_BLCA,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.blca <- data.frame(clinical.blca)
clinical.blca$PatientID <- rownames(clinical.blca)

df.merged.BLCA$PatientID <- substr(df.merged.BLCA$Sample, start = 1, stop = 12)
df.merged.BLCA$PatientAge <- clinical.blca[match(df.merged.BLCA$PatientID, clinical.blca$PatientID), "yearstobirth"]
df.merged.BLCA$PatientAge <- as.numeric(df.merged.BLCA$PatientAge)


# testing correlation of Patient age against MP and DS
cor.test(df.merged.BLCA$PatientAge, df.merged.BLCA$`% S1 (DS)`)
cor.test(df.merged.BLCA$PatientAge, df.merged.BLCA$`% S1 (MP)`)

save(df.merged.BLCA, file = "merged.mTORC1.BLCA.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("highDEPDC6", "lowDEPDC6")
df.merged.BLCA$Status <- values[df.merged.BLCA$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_BLCA <- ggboxplot(df.merged.BLCA, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_BLCA <- list(c("highDEPDC6", "lowDEPDC6"))
p_MP <- boxplot_MP_BLCA + stat_compare_means(comparisons = comparison_BLCA) + 
  stat_compare_means(label.y = 0.25, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_BLCA <- ggpar(p_MP, 
                       main = "BLCA",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.BLCA$MP ~ df.merged.BLCA$Status, df.merged.BLCA, mean)

## for % S1 DS ##
boxplot_DS_BLCA <- ggboxplot(df.merged.BLCA, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_BLCA + stat_compare_means(comparisons = comparison_BLCA) + 
  stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_BLCA <- ggpar(p_DS, 
                       main = "BLCA",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.BLCA$DS ~ df.merged.BLCA$Status, df.merged.BLCA, mean)

save(df.merged.BLCA, file = "merged.mTORC1.BLCA.RData")
