#### EXPRESSION DATA FOR MTORC1 IN LUSC ####
library(TCGA2STAT)
rnaseq.lusc <- getTCGA(disease="LUSC", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.lusc <- t(rnaseq.lusc$dat)
save(expr.lusc, file="expr.lusc.RData")
clinical.lusc <- rnaseq.lusc$clinical
save(clinical.lusc, file="clinical.lusc.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.lusc) to check whether they're present 
expr.lusc <- log2(expr.lusc+1) # log of exp data 
expr.genes.lusc <- expr.lusc[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_LUSC<-pheatmap(t(expr.genes.lusc), show_colnames = FALSE,
                 cutree_cols = 2) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_LUSC <- p_LUSC$tree_col
assignments_LUSC <- cutree(clusteredSamples_LUSC, k=2) # k = cutree_cols
groupAssignments_LUSC <- data.frame(Group=factor(assignments_LUSC))
p_LUSC <- pheatmap(t(expr.genes.lusc), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments_LUSC)

## Merge expr and group assignments 
df.expr.lusc <- data.frame(expr.genes.lusc)
df.expr.lusc$Sample <- rownames(df.expr.lusc)
groupAssignments_LUSC$SampleID <- rownames(groupAssignments_LUSC)

df.merged.LUSC <- merge(df.expr.lusc, groupAssignments_LUSC,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE) 

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.LUSC$SubID <- substr(df.merged.LUSC$Sample, start = 1, stop = 19)
S1_corr_data_LUSC$SubID <- substr(S1_corr_data_LUSC$SampleID, start = 1, stop = 19)
df.merged.LUSC <- merge(df.merged.LUSC, S1_corr_data_LUSC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.lusc <- data.frame(clinical.lusc)
clinical.lusc$PatientID <- rownames(clinical.lusc)

df.merged.LUSC$PatientID <- substr(df.merged.LUSC$Sample, start = 1, stop = 12)
df.merged.LUSC$PatientAge <- clinical.lusc[match(df.merged.LUSC$PatientID, clinical.lusc$PatientID), "yearstobirth"]
df.merged.LUSC$PatientAge <- as.numeric(df.merged.LUSC$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.LUSC$PatientAge, df.merged.LUSC$`% S1 (DS)`)
cor.test(df.merged.LUSC$PatientAge, df.merged.LUSC$`% S1 (MP)`)

save(df.merged.LUSC, file = "merged.mTORC1.LUSC.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("mediumDEPDC6", "lowDEPDC6")
df.merged.LUSC$Status <- values[df.merged.LUSC$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_LUSC <- ggboxplot(df.merged.LUSC, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_LUSC <- list(c("mediumDEPDC6", "lowDEPDC6"))
p_MP <- boxplot_MP_LUSC + stat_compare_means(comparisons = comparison_LUSC) + 
  stat_compare_means(label.y = 0.3, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_LUSC <- ggpar(p_MP, 
                       main = "LUSC",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.LUSC$MP ~ df.merged.LUSC$Status, df.merged.LUSC, mean)

## for % S1 DS ##
boxplot_DS_LUSC <- ggboxplot(df.merged.LUSC, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_LUSC + stat_compare_means(comparisons = comparison_LUSC) + 
  stat_compare_means(label.y = 1.0, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_LUSC <- ggpar(p_DS, 
                       main = "LUSC",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.LUSC$DS ~ df.merged.LUSC$Status, df.merged.LUSC, mean)

save(df.merged.LUSC, file = "merged.mTORC1.LUSC.RData")
