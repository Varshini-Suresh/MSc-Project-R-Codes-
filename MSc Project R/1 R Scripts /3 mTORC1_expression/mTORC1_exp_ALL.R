#### EXPRESSION ANALYSIS ON MERGED DATA #### 
subset1 <- data.frame(df.merged.all$Sample, df.merged.all$MTOR, df.merged.all$RPTOR, df.merged.all$DEPDC6, df.merged.all$MLST8)
subset2 <- subset1[,-1]
rownames(subset2) <- subset1[,1]
subset1 <- subset2
colnames(subset) <- c("MTOR", "RPTOR", "DEPDC6", "MLST8")

library(pheatmap)
p <-pheatmap(t(subset1), show_colnames = FALSE,
                 cutree_cols = 4)

## Assign groups ##
clusteredSamples <- p$tree_col
assignments <- cutree(clusteredSamples, k=4) # k = cutree_cols
groupAssignments <- data.frame(Group=factor(assignments))
p_label <- pheatmap(t(subset1), show_colnames = FALSE,
                   cutree_cols = 4, annotation = groupAssignments)

### Labelling with Cancer Type ###
CancerType <- data.frame(groupAssignments) # duplicating groupAssignments
# need to match Cancer column and then remove group 
data.table::setDT(CancerType, keep.rownames = TRUE) []
colnames(CancerType)[colnames(CancerType)=="rn"] <- "Sample_ID" 
CancerType$Cancer <- df.merged.all[match(CancerType$Sample_ID, df.merged.all$Sample), "Cancer"]
CancerType$Group <- NULL

# reconvert into rownames and labelling
row.names(CancerType) <- CancerType$Sample_ID
CancerType$Sample_ID <- NULL
p_cancer <- pheatmap(t(subset1), show_colnames = FALSE,
                    cutree_cols = 4, annotation = df.merged.all.filter[,c("Cancer","Group")])

### MERGING THE GROUP ASSIGNMENTS 
subset$Sample <- rownames(subset1)
groupAssignments$SampleID <- rownames(groupAssignments)
df.merged.subset <- merge(subset1, groupAssignments,
                        by.x = "Sample", by.y = "SampleID",
                        all.x = FALSE, all.y = FALSE)
# into the main merged file 
df.merged.all$Group <- NULL # remove old group
df.merged.all$Group <- df.merged.subset[match(df.merged.subset$Sample, df.merged.all$Sample), "Group"]

df.merged.all$Status <- NULL
values <- c("mediumDEPDC6", "lowDEPDC6", "highDEPDC6", "lowMLST8")
df.merged.all$Status <- values[df.merged.all$Group]

df.merged.all.filter <- subset(df.merged.all, df.merged.all$Group != 4)

## for % S1 MP ##
library(ggpubr)
boxplot_MP <- ggboxplot(df.merged.all.filter, x="Status", y="MP", color = "Status", palette = c("#7e33e8", "#552599", "#aa78f0"), add="jitter", shape = "Status")
comparison_MP <- list(c("mediumDEPDC6", "highDEPDC6"), c("highDEPDC6", "lowDEPDC6"), c("mediumDEPDC6", "lowDEPDC6"))
boxplot_MP + stat_compare_means(comparisons = comparison_MP) + 
  stat_compare_means(label.y = 0.75)
mean_MP <- aggregate(df.merged.all.filter$MP ~ df.merged.all.filter$Status, df.merged.all.filter, mean)

violin_MP <- ggviolin(df.merged.all.filter, x = "Status", y = "MP", fill = "Status", palette = c("#7e33e8", "#552599", "#aa78f0"), add = "boxplot", add.params = list(fill = "white"))
violin_MP + stat_compare_means(comparisons = comparison_MP) + 
  +     stat_compare_means(label.y = 0.75)

## for % S1 DS ##
library(ggpubr)
boxplot_DS <- ggboxplot(df.merged.all.filter, x="Status", y="DS", color = "Status", palette = c("#34ad38", "#427d44", "#0fd916"), add="jitter", shape = "Status")
comparison_DS <- list(c("mediumDEPDC6", "highDEPDC6"), c("highDEPDC6", "lowDEPDC6"), c("mediumDEPDC6", "lowDEPDC6"))
boxplot_DS + stat_compare_means(comparisons = comparison_DS) + 
  stat_compare_means(label.y = 1.75)
mean_DS <- aggregate(df.merged.all.filter$DS ~ df.merged.all.filter$Status, df.merged.all.filter, mean)

violin_DS <- ggviolin(df.merged.all.filter, x = "Status", y = "DS", fill = "Status", palette = c("#34ad38", "#427d44", "#0fd916"), add = "boxplot", add.params = list(fill = "white"))
violin_DS + stat_compare_means(comparisons = comparison_MP) + 
  +     stat_compare_means(label.y = 1.75)

#### HIGH vs LOW DEPDC6 ####
# created 2 files with subsets, added a column called "low/high" and then merged them
DEPDC6_high$DEPDC6Status <- c("High")
DEPDC6_low$DEPDC6Status <- c("Low")
DEPDC6_all <- rbind(DEPDC6_high, DEPDC6_low)

## MP 
library(ggpubr)
comparison_DEPDC6 <- list(c("High", "Low"))
violin_DEPDC6_MP <- ggviolin(DEPDC6_all, x = "DEPDC6Status", y = "MP", fill = "DEPDC6Status", palette = c("#ffb60a", "#ffed87"), add = "boxplot", add.params = list(fill = "white"))
depdc6.mp <- violin_DEPDC6_MP + stat_compare_means(comparisons = comparison_DEPDC6, size = 6) + 
  stat_compare_means(label.y = 0.62, size = 6)

depdc6.label.mp <- ggpar(depdc6.mp, 
                        xlab = "DEPDC6 expression status",
                        ylab = "% S1 (MP)",
                        font.x = c(16, "bold"), 
                        font.y = c(16, "bold"),
                        font.ytickslab = 16,
                        font.xtickslab = c(1, "white"),
                        legend = "bottom", 
                        legend.title = "DEPDC6 expression status",
                        font.legend = 14)

## MP 
library(ggpubr)
boxplot_DEPDC6_MP <- ggboxplot(DEPDC6_all, x="DEPDC6Status", y="MP", fill = "DEPDC6Status", palette = c("#552599", "#aa78f0"), shape = "DEPDC6Status")
comparison_DEPDC6 <- list(c("DEPDC6.High", "DEPDC6.Low"))
boxplot_DEPDC6_MP + stat_compare_means(comparisons = comparison_DEPDC6) + 
  stat_compare_means(label.y = 0.6)
mean_DEPDC6_MP <- aggregate(DEPDC6_all$MP ~ DEPDC6_all$DEPDC6Status, DEPDC6_all, mean)

## DS
boxplot_DEPDC6_DS <- ggboxplot(DEPDC6_all, x="DEPDC6Status", y="DS", fill = "DEPDC6Status", palette = c("#427d44", "#0fd916"), shape = "DEPDC6Status")
boxplot_DEPDC6_DS + stat_compare_means(comparisons = comparison_DEPDC6) + 
  stat_compare_means(label.y = 1.3)
mean_DEPDC6_DS <- aggregate(DEPDC6_all$DS ~ DEPDC6_all$DEPDC6Status, DEPDC6_all, mean)

violin_DEPDC6_DS <- ggviolin(DEPDC6_all, x = "DEPDC6Status", y = "DS", fill = "DEPDC6Status", palette = c("#427d44", "#0fd916"), add = "boxplot", add.params = list(fill = "white"))
violin_DEPDC6_DS + stat_compare_means(comparisons = comparison_DEPDC6) + 
  stat_compare_means(label.y = 1.3)
