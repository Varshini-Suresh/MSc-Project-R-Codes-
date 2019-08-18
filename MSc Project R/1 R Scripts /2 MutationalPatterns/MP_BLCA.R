#### EXTRACTING MUT SIG FOR BLCA ####
library(MutationalPatterns)
snvs_toconvert_BLCA <- snvs_BLCA[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_BLCA) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_BLCA$end <- snvs_toconvert_BLCA$start
snvs_toconvert_BLCA$mutation <- apply(snvs_toconvert_BLCA[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_BLCA <- makeGRangesListFromDataFrame(snvs_toconvert_BLCA, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_BLCA <- mut_matrix(vcf_list = vcfs_BLCA, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_BLCA <- mut_mat_BLCA + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_BLCA <- nmf(mut_mat_BLCA, rank=2:10, method="brunet", nrun=100, seed=123456)
plot(estimate_BLCA) # optimal rank chosen as 6
nmf_res_BLCA <-extract_signatures(mut_mat_BLCA, rank = 6, nrun = 100)

# compare with COSMIC data 
cos_sim_samples_signatures_BLCA = cos_sim_matrix(nmf_res_BLCA$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_BLCA[1,]) # from 1 to X

# plotting the signature contribution graph 
colnames(nmf_res_BLCA$signatures) <- c("Sig 5", "Sig 10", "Sig 2", "Sig 4", "Sig 13", "Sig 1") 
rownames(nmf_res_BLCA$contribution) <- c("Sig 5", "Sig 10", "Sig 2", "Sig 4", "Sig 13", "Sig 1")
plot_96_profile(nmf_res_BLCA$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_BLCA <- as.data.frame(t(nmf_res_BLCA$contribution))
for (i in 1:412) S1_rel_contr_BLCA$total[i] = sum(S1_rel_contr_BLCA[i,])
for (i in 1:412) S1_rel_contr_BLCA$S1cont[i] = S1_rel_contr_BLCA[i,6]/(S1_rel_contr_BLCA[i,7])
S1_rel_contr_BLCA$Sample <- rownames(S1_rel_contr_BLCA)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_BLCA <- data.frame(matrix("", ncol=1, nrow = 412)) #nrow is [i]
S1_corr_data_BLCA$SampleID <- S1_rel_contr_BLCA$Sample
S1_corr_data_BLCA$matrix.....ncol...1..nrow...412. <- NULL

S1_corr_data_BLCA$'% S1 (DS)' <- sigs.BLCA.dup[match(S1_corr_data_BLCA$SampleID, sigs.BLCA.dup$Sample), "Signature.1"]
S1_corr_data_BLCA$'% S1 (MP)' <- S1_rel_contr_BLCA[match(S1_corr_data_BLCA$SampleID, S1_rel_contr_BLCA$Sample), "S1cont"]

cor.test(S1_corr_data_BLCA$'% S1 (MP)', S1_corr_data_BLCA$'% S1 (DS)')
save(S1_corr_data_BLCA, file="S1.corr.BLCA.RData")

# scatterplot
library(ggplot2)
scat_BLCA <- ggplot(S1_corr_data_BLCA, aes(x=S1_corr_data_BLCA$`% S1 (MP)`, y=S1_corr_data_BLCA$`% S1 (DS)`)) + geom_point(shape=1)
scat_BLCA + labs(x= "% S1 MP", y="% S1 DS")



