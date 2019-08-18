#### Plotting SCORE ####
quantile(df.merged.subset$Score)
Score_low <- subset(df.merged.subset, df.merged.subset$Score <= 9.481132)
Score_low$Scorestatus <- c("Low")
Score_medium <- subset(df.merged.subset, df.merged.subset$Score >= 9.481132 & df.merged.subset$Score <= 9.740038)
Score_medium$Scorestatus <- c("Medium")
Score_high <- subset(df.merged.subset, df.merged.subset$Score >= 9.740038)
Score_high$Scorestatus <- c("High")

Score_subset <- rbind(Score_high, Score_medium, Score_low)

## Violin plots 
comparison_Score <- list(c("High", "Medium"), c("Medium", "Low"), c("Low", "High"))
library(ggpubr)

## MP
violin_Score_MP <- ggviolin(Score_all, x = "Scorestatus", y = "MP", fill = "Scorestatus", palette = c("#677d37", "#8ebd24", "#d7f299"), add = "boxplot", add.params = list(fill = "white"))
p.score.mp <- violin_Score_MP + stat_compare_means(comparisons = comparison_Score, size = 6) + 
  stat_compare_means(label.y = 1.5, size = 6)
p.score.mp$layers[[2]]$aes_params$textsize <- 6
lab.score.mp <- ggpar(p.score.mp, 
                      main = "mTOR Score",
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

## DS
boxplot.score.ds <- ggboxplot(Score_all, x="Scorestatus", y="DS", fill = "Scorestatus", palette = c("light blue", "pink", "green"), shape = "Scorestatus")
p.score.ds <- boxplot.score.ds + stat_compare_means(comparisons = comparison_Score) + 
  stat_compare_means(label.y = 1.6, size = 5)

mean.score.ds <- aggregate(Score_all$DS ~ Score_all$Scorestatus, Score_all, mean)

#### HIGH VS LOW MAPK15 ####
# Subsets 
quantile(df.merged.all$MAPK15)
MAPK15_low <- subset(df.merged.all, df.merged.all$MAPK15 <= 3.627618)
MAPK15_low$MAPK15status <- c("Low")
MAPK15_high <- subset(df.merged.all, df.merged.all$MAPK15 >= 6.900668)
MAPK15_high$MAPK15status <- c("High")

MAPK15_all <- rbind(MAPK15_high, MAPK15_low)

# Plots 
library(ggpubr)
comparison_MAPK15 <- list(c("High", "Low"))
## MP 
violin.mp.MAPK15 <- ggviolin(MAPK15_all, x = "MAPK15status", y = "MP", fill = "MAPK15status", palette = c("#4287f5", "#8ec6fa"), add = "boxplot", add.params = list(fill = "white"))
MAPK15.mp <- violin.mp.MAPK15 + stat_compare_means(comparisons = comparison_MAPK15, size = 6) + 
  stat_compare_means(label.y = 1.2, size = 6)
MAPK15.mp$layers[[2]]$aes_params$textsize <- 6

MAPK15.label.mp <- ggpar(MAPK15.mp, 
                        xlab = "MAPK15 expression status",
                        ylab = "% S1 (MP)",
                        font.x = c(16, "bold"), 
                        font.y = c(16, "bold"),
                        font.ytickslab = 16,
                        font.xtickslab = c(1, "white"),
                        legend = "bottom", 
                        legend.title = "MAPK15 expression status",
                        font.legend = 14)
mean.MAPK15.mp <- aggregate(MAPK15_all$MP ~ MAPK15_all$MAPK15status, MAPK15_all, mean)

## DS
violin.ds.MAPK15 <- ggviolin(MAPK15_all, x = "MAPK15status", y = "DS", fill = "MAPK15status", palette = c("#4287f5", "#8ec6fa"), add = "boxplot", add.params = list(fill = "white"))
MAPK15.ds <- violin.ds.MAPK15 + stat_compare_means(comparisons = comparison_MAPK15, size = 6) + 
  stat_compare_means(label.y = 1.4, size = 6)
MAPK15.ds$layers[[2]]$aes_params$textsize <- 6
MAPK15.label.ds <- ggpar(MAPK15.ds, 
                         xlab = "MAPK15 expression status",
                         ylab = "% S1 (DS)",
                         font.x = c(16, "bold"), 
                         font.y = c(16, "bold"),
                         font.ytickslab = 16,
                         font.xtickslab = c(1, "white"),
                         legend = "bottom", 
                         legend.title = "MAPK15 expression status",
                         font.legend = 14)
mean.MAPK15.ds <- aggregate(MAPK15_all$DS ~ MAPK15_all$MAPK15status, MAPK15_all, mean)

#### HIGH VS LOW PRKAA2 ####
# Subsets 
quantile(df.merged.all$PRKAA2)
PRKAA2_low <- subset(df.merged.all, df.merged.all$PRKAA2 <= 5.514670)
PRKAA2_low$PRKAA2status <- c("Low")
PRKAA2_high <- subset(df.merged.all, df.merged.all$PRKAA2 >= 8.695613)
PRKAA2_high$PRKAA2status <- c("High")

PRKAA2_all <- rbind(PRKAA2_high, PRKAA2_low)

# Plots 
library(ggpubr)
comparison_PRKAA2 <- list(c("High", "Low"))
## MP 
violin.mp.PRKAA2 <- ggviolin(PRKAA2_all, x = "PRKAA2status", y = "MP", fill = "PRKAA2status", palette = c("#825b44", "#cc9767"), add = "boxplot", add.params = list(fill = "white"))
PRKAA2.mp <- violin.mp.PRKAA2 + stat_compare_means(comparisons = comparison_PRKAA2, size = 6) + 
  stat_compare_means(label.y = 1.2, size = 6)
PRKAA2.mp$layers[[2]]$aes_params$textsize <- 6

PRKAA2.label.mp <- ggpar(PRKAA2.mp, 
                         xlab = "PRKAA2 expression status",
                         ylab = "% S1 (MP)",
                         font.x = c(16, "bold"), 
                         font.y = c(16, "bold"),
                         font.ytickslab = 16,
                         font.xtickslab = c(1, "white"),
                         legend = "bottom", 
                         legend.title = "PRKAA2 expression status",
                         font.legend = 14)
mean.PRKAA2.mp <- aggregate(PRKAA2_all$MP ~ PRKAA2_all$PRKAA2status, PRKAA2_all, mean)

## DS
violin.ds.PRKAA2 <- ggviolin(PRKAA2_all, x = "PRKAA2status", y = "DS", fill = "PRKAA2status", palette = c("#825b44", "#cc9767"), add = "boxplot", add.params = list(fill = "white"))
PRKAA2.ds <- violin.ds.PRKAA2 + stat_compare_means(comparisons = comparison_PRKAA2, size = 6) + 
  stat_compare_means(label.y = 1.4, size = 6)
PRKAA2.mp$layers[[2]]$aes_params$textsize <- 6
PRKAA2.label.ds <- ggpar(PRKAA2.ds, 
                         xlab = "PRKAA2 expression status",
                         ylab = "% S1 (DS)",
                         font.x = c(16, "bold"), 
                         font.y = c(16, "bold"),
                         font.ytickslab = 16,
                         font.xtickslab = c(1, "white"),
                         legend = "bottom", 
                         legend.title = "PRKAA2 expression status",
                         font.legend = 14)

mean.PRKAA2.ds <- aggregate(MAPK15_all$DS ~ MAPK15_all$MAPK15status, MAPK15_all, mean)

