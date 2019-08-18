#### mTOR PATHWAY GENES EXPRESSION DATA - ALL CANCERS MERGED ####
df.merged.all <- as.data.frame(df.merged.brca)
df.merged.all <- rbind(df.merged.brca, df.merged.esca, df.merged.hnsc, df.merged.kirc, df.merged.luad,
                       df.merged.lusc, df.merged.ov, df.merged.paad, df.merged.stad)

expr_subset <- data.frame(df.merged.all[2:40])
sample.as.rows <- expr_subset[,-1]
rownames(sample.as.rows) <- expr_subset[,1]

library(pheatmap)
pdf("heatmap.all.pdf")
p <-pheatmap(t(sample.as.rows), show_colnames = FALSE,
                  cutree_cols = 5)
dev.off()

## Assign groups ##
clusteredSamples <- p$tree_col
assignments <- cutree(clusteredSamples, k=5) # k = cutree_cols
groupAssignments <- data.frame(Group=factor(assignments))
p_label <- pheatmap(t(sample.as.rows), show_colnames = FALSE,
                    cutree_cols = 5, annotation = groupAssignments)

# merging group assignments
groupAssignments$SampleID <- rownames(groupAssignments)
df.merged.subset <- merge(expr_subset, groupAssignments,
                          by.x = "Sample", by.y = "SampleID",
                          all.x = FALSE, all.y = FALSE)

# modifications to original file
df.merged.all$Group <- NULL
df.merged.all$Status <- NULL
df.merged.all$Group <- df.merged.subset[match(df.merged.subset$Sample, df.merged.all$Sample), "Group"]

df.merged.all <- subset(df.merged.all, df.merged.all$Group != "4" & df.merged.all$Group != "5")

values <- c("lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", 
            "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L")
df.merged.all$Status <- values[df.merged.all$Group]

#### CORRELATION OF S1 IN MTOR ACTIVITY ####
comparison.all <- list(c("lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L", "lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L"), 
                        c("lowMAPK15.lowPRKAA2.lowDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"),
                       c("lowMAPK15.highPRKAA2.highDEPDC6.highPRR5L", "highMAPK15.highPRKAA2.highDEPDC6.highPRR5L"))
                        
### MP ###
library(ggpubr)
boxplot.mp.all <- ggboxplot(df.merged.all, x="Status", y="MP", fill = "Status", palette = c("#b789f0", "#fa7393", "#ffc861"), shape = "Status")
p.mp <- boxplot.mp.all + stat_compare_means(comparisons = comparison.all) + 
  stat_compare_means(label.y = 1.3, size = 5)
p.mp$layers[[2]]$aes_params$textsize <- 5

lab.all.mp <- ggpar(p.mp, 
                     main = "ALL MERGED",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
mean.mp <- aggregate(df.merged.all$MP ~ df.merged.all$Status, df.merged.all, mean)

### DS ###
boxplot.ds.all <- ggboxplot(df.merged.all, x="Status", y="DS", fill = "Status", palette = c("#9900ff", "#ff4c77", "#ffae15"), shape = "Status")
p.ds <- boxplot.ds.all + stat_compare_means(comparisons = comparison.all) + 
  stat_compare_means(label.y = 1.5, size = 5)
p.ds$layers[[2]]$aes_params$textsize <- 5

lab.all.ds <- ggpar(p.ds, 
                     main = "ALL MERGED",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(14, "bold"), 
                     font.y = c(14, "bold"),
                     font.ytickslab = 14,
                     font.xtickslab = c(1, "white"),
                     legend = "none")

mean.ds <- aggregate(df.merged.all$DS ~ df.merged.all$Status, df.merged.all, mean)

cor.test(df.merged.all$MP, df.merged.all$Score)
cor.test(df.merged.all$DS, df.merged.all$Score)

save(df.merged.all, file = "merged.file.all.RData")
