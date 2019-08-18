#### EXTRACTING MUT SIG FOR LUSC ####
library(MutationalPatterns)
snvs_toconvert_LUSC <- snvs_LUSC[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_LUSC) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_LUSC$end <- snvs_toconvert_LUSC$start
snvs_toconvert_LUSC$mutation <- apply(snvs_toconvert_LUSC[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_LUSC <- makeGRangesListFromDataFrame(snvs_toconvert_LUSC, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_LUSC <- mut_matrix(vcf_list = vcfs_LUSC, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_LUSC <- mut_mat_LUSC + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_LUSC <- nmf(mut_mat_LUSC, rank=2:10, method="brunet", nrun=100, seed=123456) 
plot(estimate_LUSC) # optimal rank chosen as 8 
nmf_res_LUSC <-extract_signatures(mut_mat_LUSC, rank = 8 , nrun = 100) 

# compare with COSMIC data
# LOAD COSMIC.RDATA 
cos_sim_samples_signatures_LUSC = cos_sim_matrix(nmf_res_LUSC$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_LUSC[1,]) # from 1 to 8

# plotting the signature contribution graph 
colnames(nmf_res_LUSC$signatures) <- c("Sig 12", "Sig 1", "Sig 4", "Sig 15", "Sig 18", "Sig 3", "Sig 7", "Sig 13")
rownames(nmf_res_LUSC$contribution) <- c("Sig 12", "Sig 1", "Sig 4", "Sig 15", "Sig 18", "Sig 3", "Sig 7", "Sig 13")
plot_96_profile(nmf_res_LUSC$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_LUSC <- as.data.frame(t(nmf_res_LUSC$contribution))
for (i in 1:492) S1_rel_contr_LUSC$total[i] = sum(S1_rel_contr_LUSC[i,])
for (i in 1:492) S1_rel_contr_LUSC$S1cont[i] = S1_rel_contr_LUSC[i,2]/(S1_rel_contr_LUSC[i,9])
S1_rel_contr_LUSC$Sample <- rownames(S1_rel_contr_LUSC)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_LUSC <- data.frame(matrix("", ncol=1, nrow = 492)) #nrow is [i]
S1_corr_data_LUSC$SampleID <- S1_rel_contr_LUSC$Sample
S1_corr_data_LUSC$matrix.....ncol...1..nrow...492. <- NULL

S1_corr_data_LUSC$'% S1 (DS)' <- sigs.LUSC.dup[match(S1_corr_data_LUSC$SampleID, sigs.LUSC.dup$Sample), "Signature.1"]
S1_corr_data_LUSC$'% S1 (MP)' <- S1_rel_contr_LUSC[match(S1_corr_data_LUSC$SampleID, S1_rel_contr_LUSC$Sample), "S1cont"]

cor.test(S1_corr_data_LUSC$'% S1 (MP)', S1_corr_data_LUSC$'% S1 (DS)')

save(S1_corr_data_LUSC, file="S1.corr.LUSC.RData")

# scatterplot
library(ggplot2)
scat_LUSC <- ggplot(S1_corr_data_LUSC, aes(x=S1_corr_data_LUSC$`% S1 (MP)`, y=S1_corr_data_LUSC$`% S1 (DS)`)) + geom_point(shape=1)
scat_LUSC + labs(x= "% S1 MP", y="% S1 DS")




