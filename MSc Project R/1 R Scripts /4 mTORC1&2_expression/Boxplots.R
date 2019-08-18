## BOXPLOTS
## MP ##
library(ggpubr)
boxplot.mp <- ggboxplot(df.merged.stad, x="Status", y="MP", fill = "Status", palette = c("#95d108", "#f291d0"), shape = "Status")
p.mp <- boxplot.mp + stat_compare_means(comparisons = comparison.stad) + 
  stat_compare_means(label.y = 0.5, size = 6)
p.mp$layers[[2]]$aes_params$textsize <- 6

lab.stad.mp <- ggpar(p.mp, 
                     main = "STAD",
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

## DS ##
boxplot.ds <- ggboxplot(df.merged.stad, x="Status", y="DS", fill = "Status", palette = c("#95d108", "#f291d0"), shape = "Status")
p.ds <- boxplot.ds + stat_compare_means(comparisons = comparison.stad) + 
  stat_compare_means(label.y = 1.35, size = 6)
p.ds$layers[[2]]$aes_params$textsize <- 6

lab.stad.ds <- ggpar(p.ds, 
                     main = "STAD",
                     font.main = c(16, "bold"),
                     xlab = "Gene expression status",
                     ylab = "% S1 (DS)",
                     font.x = c(16, "bold"), 
                     font.y = c(16, "bold"),
                     font.ytickslab = 16,
                     font.xtickslab = c(1, "white"),
                     legend = "none")
