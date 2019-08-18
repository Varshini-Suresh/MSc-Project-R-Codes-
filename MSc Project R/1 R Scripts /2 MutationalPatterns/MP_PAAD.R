#### EXTRACTING MUT SIG FOR PAAD ####
library(MutationalPatterns)
snvs_toconvert_PAAD <- snvs_PAAD[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_PAAD) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_PAAD$end <- snvs_toconvert_PAAD$start
snvs_toconvert_PAAD$mutation <- apply(snvs_toconvert_PAAD[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_PAAD <- makeGRangesListFromDataFrame(snvs_toconvert_PAAD, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_PAAD <- mut_matrix(vcf_list = vcfs_PAAD, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_PAAD <- mut_mat_PAAD + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_PAAD <- nmf(mut_mat_PAAD, rank=2:10, method="brunet", nrun=100, seed=123456) # took about 5 hours to run 
plot(estimate_PAAD) # optimal rank chosen as 5
nmf_res_PAAD <-extract_signatures(mut_mat_PAAD, rank = 5 , nrun = 100) 

# compare with COSMIC data 
cos_sim_samples_signatures_PAAD = cos_sim_matrix(nmf_res_PAAD$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_PAAD[1,]) # from 1 to X

# plotting the signature contribution graph 
colnames(nmf_res_PAAD$signatures) <- c("Sig 4", "Sig 30", "Sig 3", "Sig 1", "Sig 14") 
rownames(nmf_res_PAAD$contribution) <- c("Sig 4", "Sig 30", "Sig 3", "Sig 1", "Sig 14")
plot_96_profile(nmf_res_PAAD$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_PAAD <- as.data.frame(t(nmf_res_PAAD$contribution))
for (i in 1:177) S1_rel_contr_PAAD$total[i] = sum(S1_rel_contr_PAAD[i,])
for (i in 1:177) S1_rel_contr_PAAD$S1cont[i] = S1_rel_contr_PAAD[i,4]/(S1_rel_contr_PAAD[i,6])  
S1_rel_contr_PAAD$Sample <- rownames(S1_rel_contr_PAAD)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_PAAD <- data.frame(matrix("", ncol=1, nrow = 177)) #nrow is [i]
S1_corr_data_PAAD$SampleID <- S1_rel_contr_PAAD$Sample
S1_corr_data_PAAD$matrix.....ncol...1..nrow...177. <- NULL

S1_corr_data_PAAD$'% S1 (DS)' <- sigs.PAAD.dup[match(S1_corr_data_PAAD$SampleID, sigs.PAAD.dup$Sample), "Signature.1"]
S1_corr_data_PAAD$'% S1 (MP)' <- S1_rel_contr_PAAD[match(S1_corr_data_PAAD$SampleID, S1_rel_contr_PAAD$Sample), "S1cont"]

cor.test(S1_corr_data_PAAD$'% S1 (MP)', S1_corr_data_PAAD$'% S1 (DS)')
S1_corr_data_PAAD_filter <- subset(S1_corr_data_PAAD, S1_corr_data_PAAD$`% S1 (MP)` <0.3)
save(S1_corr_data_PAAD_filter, file="S1.corr.PAAD_filter.RData")

# scatterplot
library(ggplot2)
scat_PAAD <- ggplot(S1_corr_data_PAAD_filter, aes(x=S1_corr_data_PAAD_filter$`% S1 (MP)`, y=S1_corr_data_PAAD_filter$`% S1 (DS)`)) + geom_point(shape=1)
scat_PAAD + labs(x= "% S1 MP", y="% S1 DS")
