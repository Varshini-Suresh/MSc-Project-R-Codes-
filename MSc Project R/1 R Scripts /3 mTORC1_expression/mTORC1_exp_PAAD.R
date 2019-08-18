#### EXPRESSION DATA FOR MTORC1 IN PAAD ####
library(TCGA2STAT)
rnaseq.paad <- getTCGA(disease="PAAD", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.paad <- t(rnaseq.paad$dat)
save(expr.paad, file="expr.paad.RData")
clinical.paad <- rnaseq.paad$clinical
clinical.paad <- data.frame(clinical.paad)
save(clinical.paad, file="clinical.paad.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.paad) to check whether they're present 
expr.paad <- log2(expr.paad+1) # log of exp data 
expr.genes.paad <- expr.paad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_PAAD <-pheatmap(t(expr.genes.paad), show_colnames = FALSE,
                 cutree_cols = 3) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_PAAD <- p_PAAD$tree_col
assignments_PAAD <- cutree(clusteredSamples_PAAD, k=3) # k = cutree_cols
groupAssignments_PAAD <- data.frame(Group=factor(assignments_PAAD))
p_PAAD <- pheatmap(t(expr.genes.paad), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments_PAAD)

## Merge expr and group assignments 
df.expr.PAAD <- data.frame(expr.genes.paad)
df.expr.PAAD$Sample <- rownames(df.expr.PAAD)
groupAssignments_PAAD$SampleID <- rownames(groupAssignments_PAAD)

df.merged.PAAD <- merge(df.expr.PAAD, groupAssignments_PAAD,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.PAAD$SubID <- substr(df.merged.PAAD$Sample, start = 1, stop = 19)
S1_corr_data_PAAD_filter$SubID <- substr(S1_corr_data_PAAD_filter$SampleID, start = 1, stop = 19)
df.merged.PAAD <- merge(df.merged.PAAD, S1_corr_data_PAAD_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.paad$PatientID <- rownames(clinical.paad)

df.merged.PAAD$PatientID <- substr(df.merged.PAAD$Sample, start = 1, stop = 12)
df.merged.PAAD$PatientAge <- clinical.paad[match(df.merged.PAAD$PatientID, clinical.paad$PatientID), "yearstobirth"]
df.merged.PAAD$PatientAge <- as.numeric(df.merged.PAAD$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.PAAD$PatientAge, df.merged.PAAD$`% S1 (DS)`)
cor.test(df.merged.PAAD$PatientAge, df.merged.PAAD$`% S1 (MP)`)

save(df.merged.PAAD, file = "merged.mTORC1.PAAD.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("highDEPDC6", "lowDEPDC6", "lowRPTOR") # refer to heatmap for order
df.merged.PAAD$Status <- values[df.merged.PAAD$Group]
df.merged.PAAD.filter <- subset(df.merged.PAAD, df.merged.PAAD$Group != "3")

## for % S1 MP ##
library(ggpubr)
boxplot_MP_PAAD <- ggboxplot(df.merged.PAAD.filter, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_PAAD <- list(c("highDEPDC6", "lowDEPDC6"))
p_MP <- boxplot_MP_PAAD + stat_compare_means(comparisons = comparison_PAAD) + 
  stat_compare_means(label.y = 0.4, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_PAAD <- ggpar(p_MP, 
                       main = "PAAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")


mean_MP <- aggregate(df.merged.PAAD.filter$MP ~ df.merged.PAAD.filter$Status, df.merged.PAAD.filter, mean)

## for % S1 DS ##
boxplot_DS_PAAD <- ggboxplot(df.merged.PAAD.filter, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_PAAD + stat_compare_means(comparisons = comparison_PAAD) + 
  stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_PAAD <- ggpar(p_DS, 
                       main = "PAAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.PAAD.filter$DS ~ df.merged.PAAD.filter$Status, df.merged.PAAD.filter, mean)

save(df.merged.PAAD, file = "merged.mTORC1.PAAD.RData")
save(df.merged.PAAD.filter, file = "merged.mTORC1.PAAD.filter.RData")

