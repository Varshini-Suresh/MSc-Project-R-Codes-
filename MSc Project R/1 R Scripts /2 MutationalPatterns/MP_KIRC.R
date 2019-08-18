#### EXTRACTING MUT SIG FOR KIRC ####
library(MutationalPatterns)
snvs_toconvert_KIRC <- snvs_KIRC[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_KIRC) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_KIRC$end <- snvs_toconvert_KIRC$start
snvs_toconvert_KIRC$mutation <- apply(snvs_toconvert_KIRC[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_KIRC <- makeGRangesListFromDataFrame(snvs_toconvert_KIRC, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_KIRC <- mut_matrix(vcf_list = vcfs_KIRC, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_KIRC <- mut_mat_KIRC + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_KIRC <- nmf(mut_mat_KIRC, rank=2:10, method="brunet", nrun=100, seed=123456)
plot(estimate_KIRC) # optimal rank chosen as 4
nmf_res_KIRC <-extract_signatures(mut_mat_KIRC, rank = 4, nrun = 100)

# compare with COSMIC data 
cos_sim_samples_signatures_KIRC = cos_sim_matrix(nmf_res_KIRC$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_KIRC[1,]) # from 1 to X

# plotting the signature contribution graph 
colnames(nmf_res_KIRC$signatures) <- c("Sig 1", "Sig 3", "Sig 20", "Sig 5") 
rownames(nmf_res_KIRC$contribution) <- c("Sig 1", "Sig 3", "Sig 20", "Sig 5")
plot_96_profile(nmf_res_KIRC$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_KIRC <- as.data.frame(t(nmf_res_KIRC$contribution))
for (i in 1:336) S1_rel_contr_KIRC$total[i] = sum(S1_rel_contr_KIRC[i,])
for (i in 1:336) S1_rel_contr_KIRC$S1cont[i] = S1_rel_contr_KIRC[i,1]/(S1_rel_contr_KIRC[i,5]) 
S1_rel_contr_KIRC$Sample <- rownames(S1_rel_contr_KIRC)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_KIRC <- data.frame(matrix("", ncol=1, nrow = 336)) #nrow is [i]
S1_corr_data_KIRC$SampleID <- S1_rel_contr_KIRC$Sample
S1_corr_data_KIRC$matrix.....ncol...1..nrow...336. <- NULL

S1_corr_data_KIRC$'% S1 (DS)' <- sigs.KIRC.dup[match(S1_corr_data_KIRC$SampleID, sigs.KIRC.dup$Sample), "Signature.1"]
S1_corr_data_KIRC$'% S1 (MP)' <- S1_rel_contr_KIRC[match(S1_corr_data_KIRC$SampleID, S1_rel_contr_KIRC$Sample), "S1cont"]

cor.test(S1_corr_data_KIRC_filter$'% S1 (MP)', S1_corr_data_KIRC_filter$'% S1 (DS)')
save(S1_corr_data_KIRC, file="S1.corr.KIRC.RData")

# scatterplot
library(ggplot2)
scat_KIRC <- ggplot(S1_corr_data_KIRC_filter, aes(x=S1_corr_data_KIRC_filter$`% S1 (MP)`, y=S1_corr_data_KIRC_filter$`% S1 (DS)`)) + geom_point(shape=1)
scat_KIRC + labs(x= "% S1 MP", y="% S1 DS")

# You might need this step: 
S1_corr_data_KIRC_filter <- subset(S1_corr_data_KIRC, S1_corr_data_KIRC$`% S1 (DS)` <0.50)


