#### PATIENT stadERALL SURVIVAL (%) ####
# open the df.merged with all genes and 

# Addition of a few more columns from clinical data
clinical.stad$PatientID <- rownames(clinical.stad)
df.merged.stad$VitalStatus <- clinical.stad[match(df.merged.stad$PatientID, clinical.stad$PatientID), "vitalstatus"]
df.merged.stad$VitalStatus<- as.numeric(as.character(df.merged.stad$VitalStatus))
df.merged.stad$DaystoDeath <- clinical.stad[match(df.merged.stad$PatientID, clinical.stad$PatientID), "daystodeath"]
df.merged.stad$DaystoDeath<- as.numeric(as.character(df.merged.stad$DaystoDeath))
df.merged.stad$Followup <- clinical.stad[match(df.merged.stad$PatientID, clinical.stad$PatientID), "daystolastfollowup"]
df.merged.stad$Followup<- as.numeric(as.character(df.merged.stad$Followup))

## Calculating staderall survival 
df.merged.stad$OS <- NA
for (i in 1:nrow(df.merged.stad)) {
  if (!is.na(df.merged.stad[i,]$DaystoDeath)) {
    df.merged.stad[i,]$OS <- df.merged.stad[i,]$DaystoDeath
  } else {
    if (!is.na(df.merged.stad[i,]$Followup)) {
      df.merged.stad[i,]$OS <- df.merged.stad[i,]$Followup
    }}}

colnames(df.merged.stad)
df.merged.stad$OS [8:23]
rm(clinical.stad)
save(df.merged.stad, file = "merged.file.stad.RData")

df.merged.all <- rbind(df.merged.brca, df.merged.esca)

#########################################################

### Plotting survival 
require("survival")
library(survminer)
fit.mTOR.score <- survfit(Surv(OS, VitalStatus) ~ Scorestatus, data = Score_subset)

surv.mTOR <- ggsurvplot(fit.mTOR.score,
                        pval = TRUE,
                        risk.table = TRUE, risk.table.col = "strata", 
                        ggtheme = theme_bw(), 
                        ylab = "Overall Survival (%)", font.y = 16,
                        xlab = "Time (days)", font.x = 16, 
                        font.ticks.lab = 14)

