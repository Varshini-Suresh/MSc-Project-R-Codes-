#### EXTRACTING MUT SIG FOR STAD ####
library(MutationalPatterns)
snvs_toconvert_STAD <- snvs_STAD[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_STAD) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_STAD$end <- snvs_toconvert_STAD$start
snvs_toconvert_STAD$mutation <- apply(snvs_toconvert_STAD[,c("seqnames","start","REF","ALT")],
                                    1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                          ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_STAD <- makeGRangesListFromDataFrame(snvs_toconvert_STAD, split.field = "sample",
                                        keep.extra.columns = TRUE,
                                        names.field = "mutation")

mut_mat_STAD <- mut_matrix(vcf_list = vcfs_STAD, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_STAD <- mut_mat_STAD + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_STAD <- nmf(mut_mat_STAD, rank=2:10, method="brunet", nrun=100, seed=123456) # took about 5 hours to run 
plot(estimate_STAD) # optimal rank chosen as 8
nmf_res_STAD <-extract_signatures(mut_mat_STAD, rank = 8 , nrun = 100) 

# compare with COSMIC data 
cos_sim_samples_signatures_STAD = cos_sim_matrix(nmf_res_STAD$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_STAD[1,]) # from 1 to 5

# plotting the signature contribution graph 
colnames(nmf_res_STAD$signatures) <- c("Sig 17", "Sig 15", "Sig 20", "Sig 6", "Sig 3", "Sig 10", "Sig 1", "Sig 21") 
rownames(nmf_res_STAD$contribution) <- c("Sig 17", "Sig 15", "Sig 20", "Sig 6", "Sig 3", "Sig 10", "Sig 1", "Sig 21")
plot_96_profile(nmf_res_STAD$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_STAD <- as.data.frame(t(nmf_res_STAD$contribution))
for (i in 1:437) S1_rel_contr_STAD$total[i] = sum(S1_rel_contr_STAD[i,])
for (i in 1:437) S1_rel_contr_STAD$S1cont[i] = S1_rel_contr_STAD[i,7]/(S1_rel_contr_STAD[i,9])
S1_rel_contr_STAD$Sample <- rownames(S1_rel_contr_STAD)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_STAD <- data.frame(matrix("", ncol=1, nrow = 437)) #nrow is [i]
S1_corr_data_STAD$SampleID <- S1_rel_contr_STAD$Sample
S1_corr_data_STAD$matrix.....ncol...1..nrow...437. <- NULL

S1_corr_data_STAD$'% S1 (DS)' <- sigs.STAD.dup[match(S1_corr_data_STAD$SampleID, sigs.STAD.dup$Sample), "Signature.1"]
S1_corr_data_STAD$'% S1 (MP)' <- S1_rel_contr_STAD[match(S1_corr_data_STAD$SampleID, S1_rel_contr_STAD$Sample), "S1cont"]

cor.test(S1_corr_data_STAD$'% S1 (MP)', S1_corr_data_STAD$'% S1 (DS)')

save(S1_corr_data_STAD, file="S1.corr.STAD.RData")

# scatterplot
library(ggplot2)
scat_STAD <- ggplot(S1_corr_data_STAD, aes(x=S1_corr_data_STAD$`% S1 (MP)`, y=S1_corr_data_STAD$`% S1 (DS)`)) + geom_point(shape=1)
scat_STAD + labs(x= "% S1 MP", y="% S1 DS")




