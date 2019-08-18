#### EXPRESSION DATA FOR MTORC1 IN KIRC ####
library(TCGA2STAT)
rnaseq.kirc <- getTCGA(disease="KIRC", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.kirc <- t(rnaseq.kirc$dat)
save(expr.kirc, file="expr.kirc.RData")
clinical.kirc <- rnaseq.kirc$clinical
save(clinical.kirc, file="clinical.kirc.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.kirc) to check whether they're present 
expr.kirc <- log2(expr.kirc+1) # log of exp data 
expr.genes.kirc <- expr.kirc[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_KIRC<-pheatmap(t(expr.genes.kirc), show_colnames = FALSE,
                 cutree_cols = 3) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_KIRC <- p_KIRC$tree_col
assignments_KIRC <- cutree(clusteredSamples_KIRC, k=3) # k = cutree_cols
groupAssignments_KIRC <- data.frame(Group=factor(assignments_KIRC))
p_KIRC <- pheatmap(t(expr.genes.kirc), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments_KIRC)

## Merge expr and group assignments 
df.expr.KIRC <- data.frame(expr.genes.kirc)
df.expr.KIRC$Sample <- rownames(df.expr.KIRC)
groupAssignments_KIRC$SampleID <- rownames(groupAssignments_KIRC)

df.merged.KIRC <- merge(df.expr.KIRC, groupAssignments_KIRC,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)


#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.KIRC$SubID <- substr(df.merged.KIRC$Sample, start = 1, stop = 19)
S1_corr_data_KIRC$SubID <- substr(S1_corr_data_KIRC$SampleID, start = 1, stop = 19)
df.merged.KIRC <- merge(df.merged.KIRC, S1_corr_data_KIRC,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.kirc <- data.frame(clinical.kirc)
clinical.kirc$PatientID <- rownames(clinical.kirc)

df.merged.KIRC$PatientID <- substr(df.merged.KIRC$Sample, start = 1, stop = 12)
df.merged.KIRC$PatientAge <- clinical.kirc[match(df.merged.KIRC$PatientID, clinical.kirc$PatientID), "yearstobirth"]
df.merged.KIRC$PatientAge <- as.numeric(df.merged.KIRC$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.KIRC$PatientAge, df.merged.KIRC$`% S1 (DS)`)
cor.test(df.merged.KIRC$PatientAge, df.merged.KIRC$`% S1 (MP)`)

save(df.merged.KIRC, file = "merged.mTORC1.KIRC.RData")

#### S1 in mTORC1 ACTIVITY ####
# creating a column containing gene expression status
values <- c("mediumDEPDC6", "highDEPDC6", "lowMLST8") # refer to heatmap for order
df.merged.KIRC$Status <- values[df.merged.KIRC$Group]
df.merged.KIRC.filter <- subset(df.merged.KIRC, df.merged.KIRC$Group != "3")

## for % S1 MP ##
library(ggpubr)
boxplot_MP_KIRC <- ggboxplot(df.merged.KIRC.filter, x="Status", y="MP", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
comparison_KIRC <- list(c("mediumDEPDC6", "highDEPDC6", "lowMLST8"))
p_MP <- boxplot_MP_KIRC + stat_compare_means(comparisons = comparison_KIRC) + 
  stat_compare_means(label.y = 0.3, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_KIRC <- ggpar(p_MP, 
                       main = "KIRC",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.KIRC.filter$MP ~ df.merged.KIRC.filter$Status, df.merged.KIRC.filter, mean)

## for % S1 DS ##
boxplot_DS_KIRC <- ggboxplot(df.merged.KIRC.filter, x="Status", y="DS", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
p_DS <- boxplot_DS_KIRC + stat_compare_means(comparisons = comparison_KIRC) + 
  stat_compare_means(label.y = 1.1, size= 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_KIRC <- ggpar(p_DS, 
                       main = "KIRC",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")
mean_DS <- aggregate(df.merged.KIRC.filter$DS ~ df.merged.KIRC.filter$Status, df.merged.KIRC.filter, mean)

save(df.merged.KIRC, file = "merged.mTORC1.KIRC.RData")
save(df.merged.KIRC.filter, file = "merged.mTORC1.KIRC.filter.RData")
