#### EXPRESSION DATA FOR MTORC1 IN ESCA ####
library(TCGA2STAT)
rnaseq.esca <- getTCGA(disease="ESCA", data.type="RNASeq2", 
                       type="RPKM", clinical = TRUE) # this takes a few mins
#extract relevant matrix:
expr.esca <- t(rnaseq.esca$dat)
save(expr.esca, file="expr.esca.RData")
clinical.esca <- rnaseq.esca$clinical
clinical.esca <- data.frame(clinical.esca)
save(clinical.esca, file="clinical.esca.RData")

### Specify our list of genes of interest:
genes <- c("MTOR", "RPTOR", "DEPDC6", "MLST8") # Use "MTOR" %in% colnames(expr.esca) to check whether 
expr.esca <- log2(expr.esca+1) # log of exp data 
expr.genes.esca <- expr.esca[,genes] # to select genes from the matrix

### Visualizing heatmaps:
library(pheatmap)
p_ESCA <-pheatmap(t(expr.genes.esca), show_colnames = FALSE,
            cutree_cols = 2, fontsize = 16) # to be specified after visualising first 

# Extracting cluster assignments for each sample
clusteredSamples_ESCA <- p_ESCA$tree_col
assignments_ESCA <- cutree(clusteredSamples_ESCA, k=2) # k = cutree_cols
groupAssignments_ESCA <- data.frame(Group=factor(assignments_ESCA))
p_ESCA <- pheatmap(t(expr.genes.esca), show_colnames = FALSE,
              cutree_cols = 2, annotation = groupAssignments_ESCA)

### Merge expr and group assignments 
df.expr.ESCA <- data.frame(expr.genes.esca)
df.expr.ESCA$Sample <- rownames(df.expr.ESCA)
groupAssignments_ESCA$SampleID <- rownames(groupAssignments_ESCA)

df.merged.ESCA <- merge(df.expr.ESCA, groupAssignments_ESCA,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)

#### DS vs MP #### 
# Merging % S1 DS and MP estimates  
df.merged.ESCA$SubID <- substr(df.merged.ESCA$Sample, start = 1, stop = 19)
S1_corr_data_ESCA_filter$SubID <- substr(S1_corr_data_ESCA_filter$Sample_ID, start = 1, stop = 19)

df.merged.ESCA <- merge(df.merged.ESCA, S1_corr_data_ESCA_filter,
                        by.x = "SubID", by.y = "SubID", 
                        all.x = FALSE, all.y = FALSE)

plot(df.merged.ESCA$Group, df.merged.ESCA$`rel % S1 MP`, xlab = "Group", ylab = "% S1 MP")

# Merging patient age from clinical data
clinical.esca$PatientID <- rownames(clinical.esca)
df.merged.ESCA$PatientID <- substr(df.merged.ESCA$Sample, start = 1, stop = 12)
df.merged.ESCA$PatientAge <- clinical.esca[match(df.merged.ESCA$PatientID, clinical.esca$PatientID), "yearstobirth"]
df.merged.ESCA$PatientAge <- as.numeric(df.merged.ESCA$PatientAge)

# testing correlation of Patient age against MP and DS
cor.test(df.merged.ESCA$PatientAge, df.merged.ESCA$`rel % S1 MP`)
cor.test(df.merged.ESCA$PatientAge, df.merged.ESCA$`% S1 (DS)`)

save(df.merged.ESCA, file = "merged.mTORC1.ESCA.RData")

#### S1 in mTORC1 ACTIVITY ####
# duplicated and replaced columns '% S1 (MP)' as 'MP' and '% S1 (DS)' as 'DS'  
# creating a column containing gene expression status
values <- c("highDEPDC6", "lowDEPDC6")
df.merged.ESCA$Status <- values[df.merged.ESCA$Group]

## for % S1 MP ##
library(ggpubr)
boxplot_MP_ESCA <- ggboxplot(df.merged.ESCA, x="Status", y="MP", fill = "Status", palette = c("#52a379", "#aa78f0"), shape = "Status")
comparison_ESCA <- list(c("highDEPDC6", "lowDEPDC6"))
p_MP <- boxplot_MP_ESCA + stat_compare_means(comparisons = comparison_ESCA) + 
              stat_compare_means(label.y = 0.3, size = 6) # decide label based on maximum y-axis
p_MP$layers[[2]]$aes_params$textsize <- 6

Label_MP_ESCA <- ggpar(p_MP, 
                     main = "ESCA",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(16, "bold"), 
                     font.y = c(16, "bold"),
                     font.ytickslab = 16,
                     font.xtickslab = c(1, "white"),
                     legend = "none")


mean_MP <- aggregate(df.merged.ESCA$MP ~ df.merged.ESCA$Status, df.merged.ESCA, mean)

## for % S1 DS ##
boxplot_DS_ESCA <- ggboxplot(df.merged.ESCA, x="Status", y="DS", fill = "Status", palette = c("#52a379","#aa78f0"), shape = "Status")
p_DS <- boxplot_DS_ESCA + stat_compare_means(comparisons = comparison_ESCA) + 
  stat_compare_means(label.y = 1.25, size = 6) # decide label based on maximum y-axis
p_DS$layers[[2]]$aes_params$textsize <- 6

Label_DS_ESCA <- ggpar(p_DS, 
                       main = "ESCA",
                       font.main = c(16, "bold"),
                       xlab = "Gene expression status",
                       ylab = "% S1 (DS)",
                       font.x = c(16, "bold"), 
                       font.y = c(16, "bold"),
                       font.ytickslab = 16,
                       font.xtickslab = c(1, "white"),
                       legend = "none")

mean_DS <- aggregate(df.merged.ESCA$DS ~ df.merged.ESCA$Status, df.merged.ESCA, mean)

save(df.merged.ESCA, file = "merged.mTORC1.ESCA.RData")
