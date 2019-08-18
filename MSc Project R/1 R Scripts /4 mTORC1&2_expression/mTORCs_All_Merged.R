#### mTORCs EXP DATA - ALL CANCERS MERGED ####
# Load all files
df.merged.all <- as.data.frame(df.merged.blca)
df.merged.all <- rbind(df.merged.blca, df.merged.brca, df.merged.esca, df.merged.hnsc, df.merged.luad,
                       df.merged.lusc, df.merged.ov, df.merged.paad, df.merged.paad, df.merged.prad, df.merged.stad)

expr_subset <- data.frame(df.merged.all[2:10])
sample.as.rows <- expr_subset[,-1]
rownames(sample.as.rows) <- expr_subset[,1]

library(pheatmap)
p <-pheatmap(t(sample.as.rows), show_colnames = FALSE,
             cutree_cols = 4)

## Assign groups ##
clusteredSamples <- p$tree_col
assignments <- cutree(clusteredSamples, k=4) # k = cutree_cols
groupAssignments <- data.frame(Group=factor(assignments))
p_label <- pheatmap(t(sample.as.rows), show_colnames = FALSE,
                    cutree_cols = 4, annotation = groupAssignments)

# merging group assignments
groupAssignments$SampleID <- rownames(groupAssignments)
row.names(expr_subset) <- expr_subset$Sample
df.merged.subset <- merge(expr_subset, groupAssignments,
                          by.x = "Sample", by.y = "SampleID",
                          all.x = FALSE, all.y = FALSE)

# modifications to original file
df.merged.all$Group <- NULL
df.merged.all$Status <- NULL
df.merged.all$Group <- df.merged.subset[match(df.merged.subset$SampleID, df.merged.all$Sample), "Group"]

#### S1 correlation with mTORCs ####
values <- c("highDEPDC6.lowPRR5L", "lowDEPDC6.lowPRR5L", "highDEPDC6.highPRR5L", "lowDEPDC6.highPRR5L")
df.merged.all$Status <- values[df.merged.all$Group]

## MP ##
library(ggpubr)
violin.all.mp <- ggviolin(df.merged.all, x = "Status", y = "MP", fill = "Status", palette = c("#95d108", "#f291d0", "#f29b30", "#5dd4f5"), add = "boxplot", add.params = list(fill = "white"))
comparison.all <- list(c(1,2), c(2,3), c(3,4), c(1,3), c(2,4), c(1,4))

p.mp <- violin.all.mp + stat_compare_means(comparisons = comparison.all) + 
  stat_compare_means(label.y = 0.87, size = 6)
p.mp$layers[[2]]$aes_params$textsize <- 6

lab.all.mp <- ggpar(p.mp, 
                     main = "ALL MERGED",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(16, "bold"), 
                     font.y = c(16, "bold"),
                     font.ytickslab = 16,
                     font.xtickslab = c(1, "white"),
                     legend = "bottom", 
                     legend.title = "Gene expression status", 
                     font.legend = 16)
mean.mp <- aggregate(df.merged.all$MP ~ df.merged.all$Status, df.merged.all, mean)

## DS ##
violin.all.ds <- ggviolin(df.merged.all, x = "Status", y = "DS", fill = "Status", palette = c("#2da10a", "#eb57c8", "#f27735", "#1155f5"), add = "boxplot", shape = "Status", add.params = list(fill = "white"))

p.ds <- violin.all.ds + stat_compare_means(comparisons = comparison.all) + 
  stat_compare_means(label.y = 1.85, size = 5)
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

save (df.merged.all, file = "df.merged.all.RData")

#### PRR5L SUBSET ####
quantile (df.merged.all$PRR5L)
PRR5L_low <- subset(df.merged.all, df.merged.all$PRR5L <= 6.612706)
PRR5L_low$PRR5Lstatus <- c("Low")
PRR5L_high <- subset(df.merged.all, df.merged.all$PRR5L >= 8.501653)
PRR5L_high$PRR5Lstatus <- c("High")
PRR5L_all <- rbind(PRR5L_high, PRR5L_low)

## MP 
library(ggpubr)
comparison_PRR5L <- list(c("High", "Low"))
violin_PRR5L_MP <- ggviolin(PRR5L_all, x = "PRR5Lstatus", y = "MP", fill = "PRR5Lstatus", palette = c("#9c1616", "#f26161"), add = "boxplot", add.params = list(fill = "white"))
prr5l.mp <- violin_PRR5L_MP + stat_compare_means(comparisons = comparison_PRR5L, size = 6) + 
  stat_compare_means(label.y = 0.5, size = 6)
prr5l.mp$layers[[2]]$aes_params$textsize <- 6

prr5l.label.mp <- ggpar(prr5l.mp, 
                     xlab = "PRR5L expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(16, "bold"), 
                     font.y = c(16, "bold"),
                     font.ytickslab = 16,
                     font.xtickslab = c(1, "white"),
                     legend = "bottom", 
                     legend.title = "PRR5L expression status",
                     font.legend = 14)
mean.prr5l.mp <- aggregate(PRR5L_all$MP ~ PRR5L_all$PRR5Lstatus, PRR5L_all, mean)

# DS
violin_PRR5L_DS <- ggviolin(PRR5L_all, x = "PRR5Lstatus", y = "DS", fill = "PRR5Lstatus", palette = c("#427d44", "#0fd916"), add = "boxplot", add.params = list(fill = "white"))
violin_PRR5L_DS + stat_compare_means(comparisons = comparison_PRR5L, size = 6) + 
  stat_compare_means(label.y = 1.3, size = 6)
mean.prr5l.ds <- aggregate(PRR5L_all$DS ~ PRR5L_all$PRR5Lstatus, PRR5L_all, mean)

#### mTORC SCORE ####
for (i in 1:4455) df.merged.all$Score[i] = sum(df.merged.all[i,3:10])/10
quantile (df.merged.all$Score)
Score_low <- subset(df.merged.all, df.merged.all$Score <= 7.568378)
Score_low$Scorestatus <- c("Low")
Score_medium <- subset(df.merged.all, df.merged.all$Score >= 7.568378 & df.merged.all$Score <= 7.838314)
Score_medium$Scorestatus <- c("Medium")
Score_high <- subset(df.merged.all, df.merged.all$Score >= 7.838314)
Score_high$Scorestatus <- c("High")
Score_all <- rbind(Score_high, Score_medium, Score_low)

## MP 
library(ggpubr)
comparison_Score <- list(c("High", "Medium"), c("Medium", "Low"), c("High", "Low"))
violin_Score_MP <- ggviolin(Score_all, x = "Scorestatus", y = "MP", fill = "Scorestatus", palette = c("#b00476", "#ff7af6", "#faaff5"), add = "boxplot", add.params = list(fill = "white"))
score.mp <- violin_Score_MP + stat_compare_means(comparisons = comparison_Score, size = 6) + 
  stat_compare_means(label.y = 0.85, size = 6)
score.mp$layers[[2]]$aes_params$textsize <- 6

lab.score.mp <- ggpar(score.mp, 
                     main = "mTORC Score",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (MP)",
                     font.x = c(16, "bold"), 
                     font.y = c(16, "bold"),
                     font.ytickslab = 16,
                     font.xtickslab = c(1, "white"),
                     legend = "bottom", 
                     legend.title = "Gene expression status", 
                     font.legend = 14)

mean.score.mp <- aggregate(Score_all$MP ~ Score_all$Scorestatus, Score_all, mean)

# DS
violin_Score_DS <- ggviolin(Score_all, x = "Scorestatus", y = "DS", fill = "Scorestatus", palette = c("red", "blue", "yellow"), add = "boxplot", add.params = list(fill = "white"))
violin_Score_DS + stat_compare_means(comparisons = comparison_Score) + 
  stat_compare_means(label.y = 1.3)
mean.score.ds <- aggregate(Score_all$DS ~ Score_all$Scorestatus, Score_all, mean)

