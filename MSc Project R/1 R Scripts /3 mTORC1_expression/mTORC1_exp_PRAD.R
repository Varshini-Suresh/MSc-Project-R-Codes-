#### EXPRESSION DATA FOR MTORC1 IN PRAD ####
library(TCGA2STAT)
rnaseq.prad <- getTCGA(disease="PRAD", data.type="RNASeq2", 
                     type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.prad <- t(rnaseq.prad$dat)
save(expr.prad, file="expr.prad.RData")
clinical.prad <- rnaseq.prad$clinical
clinical.prad <- data.frame(clinical.prad)
save(clinical.prad, file="clinical.prad.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.prad) to check whether they're present 
expr.prad <- log2(expr.prad+1) # log of exp data 
expr.genes.prad <- expr.prad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_PRAD <-pheatmap(t(expr.genes.prad), show_colnames = FALSE,
                cutree_cols = 3) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_PRAD <- p_PRAD$tree_col
assignments_PRAD <- cutree(clusteredSamples_PRAD, k=3) # k = cutree_cols
groupAssignments_PRAD <- data.frame(Group=factor(assignments_PRAD))
p_PRAD <- pheatmap(t(expr.genes.prad), show_colnames = FALSE,
                   cutree_cols = 3, annotation = groupAssignments_PRAD)

## Merge expr and group assignments 
df.expr.PRAD <- data.frame(expr.genes.prad)
df.expr.PRAD$Sample <- rownames(df.expr.PRAD)
groupAssignments_PRAD$SampleID <- rownames(groupAssignments_PRAD)

df.merged.PRAD <- merge(df.expr.PRAD, groupAssignments_PRAD,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)
df.merged.PRAD$Group[df.merged.PRAD$Group =="3"] <- "2"

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.PRAD$SubID <- substr(df.merged.PRAD$Sample, start = 1, stop = 19)
S1_corr_data_PRAD_filter$SubID <- substr(S1_corr_data_PRAD_filter$SampleID, start = 1, stop = 19)
df.merged.PRAD <- merge(df.merged.PRAD, S1_corr_data_PRAD_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.prad$PatientID <- rownames(clinical.prad)

df.merged.PRAD$PatientID <- substr(df.merged.PRAD$Sample, start = 1, stop = 12)
df.merged.PRAD$PatientAge <- clinical.prad[match(df.merged.PRAD$PatientID, clinical.prad$PatientID), "yearstobirth"]
df.merged.PRAD$PatientAge <- as.numeric(df.merged.PRAD$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.PRAD$PatientAge, df.merged.PRAD$`% S1 (DS)`)
cor.test(df.merged.PRAD$PatientAge, df.merged.PRAD$`% S1 (MP)`)

save(df.merged.PRAD, file = "merged.mTORC1.PRAD.RData")

#### S1 in mTORC1 ACTIVITY ####
# creating a column containing gene expression status
values <- c("lowDEPDC6", "highDEPDC6") # refer to heatmap for order
df.merged.PRAD$Status <- values[df.merged.PRAD$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_PRAD <- ggboxplot(df.merged.PRAD, x="Status", y="MP", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
comparison_PRAD <- list(c("lowDEPDC6", "highDEPDC6"))
p_MP <- boxplot_MP_PRAD + stat_compare_means(comparisons = comparison_PRAD) + 
  stat_compare_means(label.y = 0.5, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_PRAD <- ggpar(p_MP, 
                       main = "PRAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.PRAD$MP ~ df.merged.PRAD$Status, df.merged.PRAD, mean)

## for % S1 DS ##
boxplot_DS_PRAD <- ggboxplot(df.merged.PRAD, x="Status", y="DS", fill = "Status", palette = c("#aa78f0", "#52a379"), shape = "Status")
p_DS <- boxplot_DS_PRAD + stat_compare_means(comparisons = comparison_PRAD) + 
  stat_compare_means(label.y = 1.5, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_PRAD <- ggpar(p_DS, 
                       main = "PRAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")
mean_DS <- aggregate(df.merged.PRAD$DS ~ df.merged.PRAD$Status, df.merged.PRAD, mean)

save(df.merged.PRAD, file = "merged.mTORC1.PRAD.RData")
