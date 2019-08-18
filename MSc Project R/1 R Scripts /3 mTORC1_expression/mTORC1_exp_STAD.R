#### EXPRESSION DATA FOR MTORC1 IN STAD ####
library(TCGA2STAT)
rnaseq.stad <- getTCGA(disease="STAD", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.stad <- t(rnaseq.stad$dat)
save(expr.stad, file="expr.stad.RData")
clinical.stad <- rnaseq.stad$clinical
clinical.stad <- data.frame(clinical.stad)
save(clinical.stad, file="clinical.stad.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.stad) to check whether they're present 
expr.stad <- log2(expr.stad+1) # log of exp data 
expr.genes.stad <- expr.stad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_STAD <-pheatmap(t(expr.genes.stad), show_colnames = FALSE,
                  cutree_cols = 2) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_STAD <- p_STAD$tree_col
assignments_STAD <- cutree(clusteredSamples_STAD, k=2) # k = cutree_cols
groupAssignments_STAD <- data.frame(Group=factor(assignments_STAD))
p_STAD <- pheatmap(t(expr.genes.stad), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments_STAD)

## Merge expr and group assignments 
df.expr.STAD <- data.frame(expr.genes.stad)
df.expr.STAD$Sample <- rownames(df.expr.STAD)
groupAssignments_STAD$SampleID <- rownames(groupAssignments_STAD)

df.merged.STAD <- merge(df.expr.STAD, groupAssignments_STAD,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE) 

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.STAD$SubID <- substr(df.merged.STAD$Sample, start = 1, stop = 19)
S1_corr_data_STAD$SubID <- substr(S1_corr_data_STAD$SampleID, start = 1, stop = 19)
df.merged.STAD <- merge(df.merged.STAD, S1_corr_data_STAD,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.stad$PatientID <- rownames(clinical.stad)

df.merged.STAD$PatientID <- substr(df.merged.STAD$Sample, start = 1, stop = 12)
df.merged.STAD$PatientAge <- clinical.stad[match(df.merged.STAD$PatientID, clinical.stad$PatientID), "yearstobirth"]
df.merged.STAD$PatientAge <- as.numeric(df.merged.STAD$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.STAD$PatientAge, df.merged.STAD$`% S1 (DS)`)
cor.test(df.merged.STAD$PatientAge, df.merged.STAD$`% S1 (MP)`)

save(df.merged.STAD, file = "merged.mTORC1.STAD.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("highDEPDC6", "mediumDEPDC6")
df.merged.STAD$Status <- values[df.merged.STAD$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_STAD <- ggboxplot(df.merged.STAD, x="Status", y="MP", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
comparison_STAD <- list(c("highDEPDC6", "mediumDEPDC6"))
p_MP <- boxplot_MP_STAD + stat_compare_means(comparisons = comparison_STAD) + 
  stat_compare_means(label.y = 0.52, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_STAD <- ggpar(p_MP, 
                       main = "STAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")
mean_MP <- aggregate(df.merged.STAD$MP ~ df.merged.STAD$Status, df.merged.STAD, mean)

## for % S1 DS ##
boxplot_DS_STAD <- ggboxplot(df.merged.STAD, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_STAD + stat_compare_means(comparisons = comparison_STAD) + 
  stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_STAD <- ggpar(p_DS, 
                       main = "STAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.STAD$DS ~ df.merged.STAD$Status, df.merged.STAD, mean)

save(df.merged.STAD, file = "merged.mTORC1.STAD.RData")
