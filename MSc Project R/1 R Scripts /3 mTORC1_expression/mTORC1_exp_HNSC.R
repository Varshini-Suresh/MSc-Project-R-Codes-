#### EXPRESSION DATA FOR MTORC1 IN HNSC ####
library(TCGA2STAT)
rnaseq.hnsc <- getTCGA(disease="HNSC", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
# extract relevant matrix:
expr.hnsc <- t(rnaseq.hnsc$dat)
save(expr.hnsc, file="expr.hnsc.RData")
clinical.hnsc <- rnaseq.hnsc$clinical
save(clinical.hnsc, file="clinical.hnsc.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.hnsc) to check whether they're present 
expr.hnsc <- log2(expr.hnsc+1) # log of exp data 
expr.genes.hnsc <- expr.hnsc[,genes] # to select genes from the matrix

### Visualizing heatmaps 
library(pheatmap)
p_HNSC<-pheatmap(t(expr.genes.hnsc), show_colnames = FALSE,
                 cutree_cols = 2) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_HNSC <- p_HNSC$tree_col
assignments_HNSC <- cutree(clusteredSamples_HNSC, k=2) # k = cutree_cols
groupAssignments_HNSC <- data.frame(Group=factor(assignments_HNSC))
p_HNSC <- pheatmap(t(expr.genes.hnsc), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments_HNSC)

## Merge expr and group assignments 
df.expr.HNSC <- data.frame(expr.genes.hnsc)
df.expr.HNSC$Sample <- rownames(df.expr.HNSC)
groupAssignments_HNSC$SampleID <- rownames(groupAssignments_HNSC)

df.merged.HNSC <- merge(df.expr.HNSC, groupAssignments_HNSC,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.HNSC$SubID <- substr(df.merged.HNSC$Sample, start = 1, stop = 19)
S1_corr_data_HNSC$SubID <- substr(S1_corr_data_HNSC$SampleID, start = 1, stop = 19)
df.merged.HNSC <- merge(df.merged.HNSC, S1_corr_data_HNSC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.hnsc <- data.frame(clinical.hnsc)
clinical.hnsc$PatientID <- rownames(clinical.hnsc)

df.merged.HNSC$PatientID <- substr(df.merged.HNSC$Sample, start = 1, stop = 12)
df.merged.HNSC$PatientAge <- clinical.hnsc[match(df.merged.HNSC$PatientID, clinical.hnsc$PatientID), "yearstobirth"]
df.merged.HNSC$PatientAge <- as.numeric(df.merged.HNSC$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.HNSC$PatientAge, df.merged.HNSC$`% S1 (DS)`)
cor.test(df.merged.HNSC$PatientAge, df.merged.HNSC$`% S1 (MP)`)

save(df.merged.HNSC, file = "merged.mTORC1.HNSC.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("highDEPDC6", "lowDEPDC6")
df.merged.HNSC$Status <- values[df.merged.HNSC$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_HNSC <- ggboxplot(df.merged.HNSC, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_HNSC <- list(c("highDEPDC6", "lowDEPDC6"))
p_MP <- boxplot_MP_HNSC + stat_compare_means(comparisons = comparison_HNSC) + 
  stat_compare_means(label.y = 0.4, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_HNSC <- ggpar(p_MP, 
                       main = "HNSC",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.HNSC$MP ~ df.merged.HNSC$Status, df.merged.HNSC, mean)

## for % S1 DS ##
boxplot_DS_HNSC <- ggboxplot(df.merged.HNSC, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_HNSC + stat_compare_means(comparisons = comparison_HNSC) + 
  stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_HNSC <- ggpar(p_DS, 
                       main = "HNSC",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "bottom", 
                       legend.title = c("Expression status"),
                       font.legend = 16)

mean_DS <- aggregate(df.merged.HNSC$DS ~ df.merged.HNSC$Status, df.merged.HNSC, mean)

save(df.merged.HNSC, file = "merged.mTORC1.HNSC.RData")
