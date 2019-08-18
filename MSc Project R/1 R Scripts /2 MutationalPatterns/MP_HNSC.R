#### EXTRACTING MUT SIG FOR HSNC ####
library(MutationalPatterns)
snvs_toconvert_HNSC <- snvs_HNSC[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_HNSC) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_HNSC$end <- snvs_toconvert_HNSC$start
snvs_toconvert_HNSC$mutation <- apply(snvs_toconvert_HNSC[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_HNSC <- makeGRangesListFromDataFrame(snvs_toconvert_HNSC, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_HNSC <- mut_matrix(vcf_list = vcfs_HNSC, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_HNSC <- mut_mat_HNSC + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_HNSC <- nmf(mut_mat_HNSC, rank=2:10, method="brunet", nrun=100, seed=123456)
plot(estimate_HNSC) # optimal rank chosen as 5 
nmf_res_HNSC_5 <-extract_signatures(mut_mat_HNSC, rank = 5 , nrun = 100)

# compare with COSMIC data
cos_sim_samples_signatures_HNSC_5 = cos_sim_matrix(nmf_res_HNSC_5$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_HNSC_5[1,]) # from 1 to 5

# plotting the signature contribution graph 
colnames(nmf_res_HNSC_5$signatures) <- c("Sig 4", "Sig 6", "Sig 13", "Sig 1", "Sig 7") 
rownames(nmf_res_HNSC_5$contribution) <- c("Sig 4", "Sig 6", "Sig 13", "Sig 1", "Sig 7")
plot_96_profile(nmf_res_HNSC_5$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_HNSC <- as.data.frame(t(nmf_res_HNSC_5$contribution))
for (i in 1:507) S1_rel_contr_HNSC$total[i] = sum(S1_rel_contr_HNSC[i,])
for (i in 1:507) S1_rel_contr_HNSC$S1cont[i] = S1_rel_contr_HNSC[i,4]/(S1_rel_contr_HNSC[i,6])
S1_rel_contr_HNSC$Sample <- rownames(S1_rel_contr_HNSC)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_HNSC <- data.frame(matrix("", ncol=1, nrow = 507)) #nrow is [i]
S1_corr_data_HNSC$SampleID <- S1_rel_contr_HNSC$Sample
S1_corr_data_HNSC$matrix.....ncol...1..nrow...507. <- NULL

S1_corr_data_HNSC$'% S1 (DS)' <- sigs.HNSC.dup[match(S1_corr_data_HNSC$SampleID, sigs.HNSC.dup$Sample), "Signature.1"]
S1_corr_data_HNSC$'% S1 (MP)' <- S1_rel_contr_HNSC[match(S1_corr_data_HNSC$SampleID, S1_rel_contr_HNSC$Sample), "S1cont"]

cor.test(S1_corr_data_HNSC$'% S1 (MP)', S1_corr_data_HNSC$'% S1 (DS)')

save(S1_corr_data_HNSC, file="S1.corr.HNSC.RData")

# scatterplot
library(ggplot2)
scat_HNSC <- ggplot(S1_corr_data_HNSC, aes(x=S1_corr_data_HNSC$`% S1 (MP)`, y=S1_corr_data_HNSC$`% S1 (DS)`)) + geom_point(shape=1)
scat_HNSC + labs(x= "% S1 MP", y="% S1 DS")

