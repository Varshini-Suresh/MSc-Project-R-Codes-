#### mTORC1 EXPRESSION DATA FOR ALL CANCERS ####
# load all df.merged files
# create 2 columns on each 
# then merge using rbind 

df.merged.BLCA$Cancer <- "BLCA" #new column with same data
for (i in 1:408) df.merged.BLCA$Score[i] <- (df.merged.BLCA$MTOR[i] + df.merged.BLCA$RPTOR[i] + df.merged.BLCA$DEPDC6[i] + df.merged.BLCA$MLST8[i])/4
save(df.merged.BLCA, file = "merged.mTORC1.BLCA.RData") # save into folder
df.merged.all <- as.data.frame(df.merged.BLCA)

df.merged.BRCA$Cancer <- "BRCA"
for (i in 1:980) df.merged.BRCA$Score[i] <- (df.merged.BRCA$MTOR[i] + df.merged.BRCA$RPTOR[i] + df.merged.BRCA$DEPDC6[i] + df.merged.BRCA$MLST8[i])/4
save(df.merged.BRCA, file = "merged.mTORC1.BRCA.RData")
df.merged.all <- rbind(df.merged.BLCA, df.merged.BRCA)
# then removed it from workspace 

df.merged.ESCA$Cancer <- "ESCA"
for (i in 1:179) df.merged.ESCA$Score[i] <- (df.merged.ESCA$MTOR[i] + df.merged.ESCA$RPTOR[i] + df.merged.ESCA$DEPDC6[i] + df.merged.ESCA$MLST8[i])/4
  # error with merging
reorder_ESCA <- df.merged.ESCA[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 11, 13, 14, 15)]
colnames(reorder_ESCA) [colnames(reorder_ESCA)=="Sample_ID"] <- "SampleID"
df.merged.all <- rbind(df.merged.all, reorder_ESCA)
save(reorder_ESCA, file = "merged.mTORC1.ESCA.RData")

df.merged.HNSC$Cancer <- "HNSC"
for (i in 1:499) df.merged.HNSC$Score[i] <- (df.merged.HNSC$MTOR[i] + df.merged.HNSC$RPTOR[i] + df.merged.HNSC$DEPDC6[i] + df.merged.HNSC$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.HNSC)
save(df.merged.HNSC, file = "merged.mTORC1.HNSC.RData")

df.merged.KIRC$Cancer <- "KIRC"
for (i in 1:316) df.merged.KIRC$Score[i] <- (df.merged.KIRC$MTOR[i] + df.merged.KIRC$RPTOR[i] + df.merged.KIRC$DEPDC6[i] + df.merged.KIRC$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.KIRC)
save(df.merged.KIRC, file = "merged.mTORC1.KIRC.RData")

df.merged.LUAD$Cancer <- "LUAD"
for (i in 1:501) df.merged.LUAD$Score[i] <- (df.merged.LUAD$MTOR[i] + df.merged.LUAD$RPTOR[i] + df.merged.LUAD$DEPDC6[i] + df.merged.LUAD$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.LUAD)
save(df.merged.LUAD, file = "merged.mTORC1.LUAD.RData")

df.merged.LUSC$Cancer <- "LUSC"
for (i in 1:489) df.merged.LUSC$Score[i] <- (df.merged.LUSC$MTOR[i] + df.merged.LUSC$RPTOR[i] + df.merged.LUSC$DEPDC6[i] + df.merged.LUSC$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.LUSC)
save(df.merged.LUSC, file = "merged.mTORC1.LUSC.RData")

df.merged.OV$Cancer <- "OV"
for (i in 1:204) df.merged.OV$Score[i] <- (df.merged.OV$MTOR[i] + df.merged.OV$RPTOR[i] + df.merged.OV$DEPDC6[i] + df.merged.OV$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.OV)
save(df.merged.OV, file = "merged.mTORC1.OV.RData")

df.merged.PAAD$Cancer <- "PAAD"
for (i in 1:165) df.merged.PAAD$Score[i] <- (df.merged.PAAD$MTOR[i] + df.merged.PAAD$RPTOR[i] + df.merged.PAAD$DEPDC6[i] + df.merged.PAAD$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.PAAD)
save(df.merged.PAAD, file = "merged.mTORC1.PAAD.RData")

df.merged.PRAD$Cancer <- "PRAD"
for (i in 1:490) df.merged.PRAD$Score[i] <- (df.merged.PRAD$MTOR[i] + df.merged.PRAD$RPTOR[i] + df.merged.PRAD$DEPDC6[i] + df.merged.PRAD$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.PRAD)
save(df.merged.PRAD, file = "merged.mTORC1.PRAD.RData")

df.merged.STAD$Cancer <- "STAD"
for (i in 1:410) df.merged.STAD$Score[i] <- (df.merged.STAD$MTOR[i] + df.merged.STAD$RPTOR[i] + df.merged.STAD$DEPDC6[i] + df.merged.STAD$MLST8[i])/4
df.merged.all <- rbind(df.merged.all, df.merged.STAD)
save(df.merged.STAD, file = "merged.mTORC1.STAD.RData")


df.merged.stad$Cancer <- "STAD"
save(df.merged.stad, file = "df.merged.stad.RData")
df.merged.all <- rbind(df.merged.all, df.merged.stad)
rm(df.merged.stad)

