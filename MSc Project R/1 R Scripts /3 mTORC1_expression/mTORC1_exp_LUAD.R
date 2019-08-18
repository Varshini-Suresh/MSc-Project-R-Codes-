#### EXPRESSION DATA FOR MTORC1 IN LUAD ####
library(TCGA2STAT)
rnaseq.luad <- getTCGA(disease="LUAD", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.luad <- t(rnaseq.luad$dat)
save(expr.luad, file="expr.luad.RData")
clinical.luad <- rnaseq.luad$clinical
clinical.luad <- data.frame(clinical.luad)
save(clinical.luad, file="clinical.luad.RData")

### Specify list of mTORC1 genes:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.luad) to check whether they're present 
expr.luad <- log2(expr.luad+1) # log of exp data 
expr.genes.luad <- expr.luad[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_LUAD<-pheatmap(t(expr.genes.luad), show_colnames = FALSE,
                 cutree_cols = 2) # to be specified after visualising first

## Extracting cluster assignments for each sample
clusteredSamples_LUAD <- p_LUAD$tree_col
assignments_LUAD <- cutree(clusteredSamples_LUAD, k=2) # k = cutree_cols
groupAssignments_LUAD <- data.frame(Group=factor(assignments_LUAD))
p_LUAD <- pheatmap(t(expr.genes.luad), show_colnames = FALSE,
                   cutree_cols = 2, annotation = groupAssignments_LUAD)

## Merge expr and group assignments 
df.expr.luad <- data.frame(expr.genes.luad)
df.expr.luad$Sample <- rownames(df.expr.luad)
groupAssignments_LUAD$SampleID <- rownames(groupAssignments_LUAD)

df.merged.LUAD <- merge(df.expr.luad, groupAssignments_LUAD,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

#### DS vs MP #### 
## Merging % S1 DS and MP estimates
df.merged.LUAD$SubID <- substr(df.merged.LUAD$Sample, start = 1, stop = 19)
S1_corr_data_LUAD_filter$SubID <- substr(S1_corr_data_LUAD_filter$SampleID, start = 1, stop = 19)
df.merged.LUAD <- merge(df.merged.LUAD, S1_corr_data_LUAD_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

# Merging patient age from clinical data
clinical.luad <- data.frame(clinical.luad)
clinical.luad$PatientID <- rownames(clinical.luad)

df.merged.LUAD$PatientID <- substr(df.merged.LUAD$Sample, start = 1, stop = 12)
df.merged.LUAD$PatientAge <- clinical.luad[match(df.merged.LUAD$PatientID, clinical.luad$PatientID), "yearstobirth"]
df.merged.LUAD$PatientAge <- as.numeric(df.merged.LUAD$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.LUAD$PatientAge, df.merged.LUAD$`% S1 (DS)`)
cor.test(df.merged.LUAD$PatientAge, df.merged.LUAD$`% S1 (MP)`)

save(df.merged.LUAD, file = "merged.mTORC1.LUAD.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("highDEPDC6", "lowDEPDC6")
df.merged.LUAD$Status <- values[df.merged.LUAD$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_LUAD <- ggboxplot(df.merged.LUAD, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_LUAD <- list(c("highDEPDC6", "lowDEPDC6"))
p_MP <- boxplot_MP_LUAD + stat_compare_means(comparisons = comparison_LUAD) + 
  stat_compare_means(label.y = 0.35, size = 6) 
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_LUAD <- ggpar(p_MP, 
                       main = "LUAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (MP)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_MP <- aggregate(df.merged.LUAD$MP ~ df.merged.LUAD$Status, df.merged.LUAD, mean)

## for % S1 DS ##
boxplot_DS_LUAD <- ggboxplot(df.merged.LUAD, x="Status", y="DS", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_LUAD + stat_compare_means(comparisons = comparison_LUAD) + 
  stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_LUAD <- ggpar(p_DS, 
                       main = "LUAD",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.LUAD$DS ~ df.merged.LUAD$Status, df.merged.LUAD, mean)

save(df.merged.LUAD, file = "merged.mTORC1.LUAD.RData")
