#### EXTRACTING MUT SIG FOR OV ####
library(MutationalPatterns)
snvs_toconvert_OV <- snvs_OV[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_OV) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_OV$end <- snvs_toconvert_OV$start
snvs_toconvert_OV$mutation <- apply(snvs_toconvert_OV[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_OV <- makeGRangesListFromDataFrame(snvs_toconvert_OV, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_OV <- mut_matrix(vcf_list = vcfs_OV, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_OV <- mut_mat_OV + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_OV <- nmf(mut_mat_OV, rank=2:10, method="brunet", nrun=100, seed=123456) # took about 2.5 hours to run 
plot(estimate_OV) # optimal rank chosen as 5 
nmf_res_OV <-extract_signatures(mut_mat_OV, rank = 5 , nrun = 100) # took 20 min to complete

# compare with COSMIC data, 
cos_sim_samples_signatures_OV = cos_sim_matrix(nmf_res_OV$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_OV[1,]) # from 1 to 5

# plotting the signature contribution graph 
colnames(nmf_res_OV$signatures) <- c("Sig 8", "Sig 5", "Sig 1", "Sig 30", "Sig 3")
rownames(nmf_res_OV$contribution) <- c("Sig 8", "Sig 5", "Sig 1", "Sig 30", "Sig 3")
plot_96_profile(nmf_res_OV$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_OV <- as.data.frame(t(nmf_res_OV$contribution))
for (i in 1:436) S1_rel_contr_OV$total[i] = sum(S1_rel_contr_OV[i,])
for (i in 1:436) S1_rel_contr_OV$S1cont[i] = S1_rel_contr_OV[i,3]/(S1_rel_contr_OV[i,6])
S1_rel_contr_OV$Sample <- rownames(S1_rel_contr_OV)

#### CORRELATION WITH DS AND MP ####
# correlation table 
S1_corr_data_OV <- data.frame(matrix("", ncol=1, nrow = 436)) #nrow is [i]
S1_corr_data_OV$SampleID <- S1_rel_contr_OV$Sample
S1_corr_data_OV$matrix.....ncol...1..nrow...436. <- NULL

S1_corr_data_OV$'% S1 (DS)' <- sigs.OV.dup[match(S1_corr_data_OV$SampleID, sigs.OV.dup$Sample), "Signature.1"]
S1_corr_data_OV$'% S1 (MP)' <- S1_rel_contr_OV[match(S1_corr_data_OV$SampleID, S1_rel_contr_OV$Sample), "S1cont"]

# had to filter the outlier (above 0.25)
cor.test(S1_corr_data_OV_filter$'% S1 (MP)', S1_corr_data_OV_filter$'% S1 (DS)')

save(S1_corr_data_OV_, file="S1.corr.OV.RData")

# scatterplot
library(ggplot2)
scat_OV <- ggplot(S1_corr_data_OV_filter, aes(x=S1_corr_data_OV_filter$`% S1 (MP)`, y=S1_corr_data_OV_filter$`% S1 (DS)`)) + geom_point(shape=1)
scat_OV + labs(x= "% S1 MP", y="% S1 DS")



