#### EXTRACTING MUT SIG FOR ESCA ####
library(MutationalPatterns)
library(GenomicRanges)

### Convering file format ### 
snvs_toconvert_ESCA <- snvs[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                               "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_ESCA) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_ESCA$end <- snvs_toconvert_ESCA$start
snvs_toconvert_ESCA$mutation <- apply(snvs_toconvert_ESCA[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_ESCA <- makeGRangesListFromDataFrame(snvs_toconvert_ESCA, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_for_ESCA <- mut_matrix(vcf_list = vcfs_ESCA, ref_genome = ref_genome) 
mut_mat_for_ESCA <- mut_mat_for_ESCA + 0.0001

### NMF analysis ###
library("NMF")
estimate_ESCA <- nmf(mut_mat_for_ESCA, rank=2:10, method="brunet", nrun=100, seed=123456)
plot(estimate_ESCA) # optimal rank chosen is 5 
nmf_res_ESCA <-extract_signatures(mut_mat_for_ESCA, rank = 5, nrun = 100) 

# compare with COSMIC data
cos_sim_samples_signatures_ESCA = cos_sim_matrix(nmf_res_ESCA$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_ESCA[1,]) # this was done to [1,] to [5,] and noted.

# plotting the signature contribution graph 
colnames(nmf_res_ESCA$signatures) <- c("Sig 13", "Sig 1", "Sig 16", "Sig 17", "Sig 18")
rownames(nmf_res_ESCA$contribution) <- c("Sig 13", "Sig 1", "Sig 16", "Sig 17", "Sig 18")
plot_96_profile(nmf_res_ESCA$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_ESCA <- as.data.frame(t(nmf_res_ESCA$contribution))
for (i in 1:184) S1_rel_contr_ESCA$total[i] = sum(S1_rel_contr_ESCA[i,])
for (i in 1:184) S1_rel_contr_ESCA$S1cont[i] = S1_rel_contr_ESCA[i,2]/(S1_rel_contr_ESCA[i,6])
S1_rel_contr_ESCA$Sample <- rownames(S1_rel_contr_ESCA)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_ESCA <- data.frame(matrix("", ncol=1, nrow = 184)) #nrow is [i]
S1_corr_data_ESCA$SampleID <- S1_rel_contr_ESCA$Sample
S1_corr_data_ESCA$matrix.....ncol...1..nrow...184. <- NULL

S1_corr_data_ESCA$'% S1 (DS)' <- sigs.ESCA[match(S1_corr_data_ESCA$SampleID, sigs.ESCA$Sample), "Signature.1"]
S1_corr_data_ESCA$'% S1 (MP)' <- S1_rel_contr_ESCA[match(S1_corr_data_ESCA$SampleID, S1_rel_contr_ESCA$Sample), "S1cont"]

cor.test(S1_corr_data_ESCA$'% S1 (MP)', S1_corr_data_ESCA$'% S1 (DS)')

save(S1_corr_data_ESCA, file="S1.corr.ESCA.RData")

# scatterplot
library(ggplot2)
scat_ESCA <- ggplot(S1_corr_data_ESCA, aes(x=S1_corr_data_ESCA$`% S1 (MP)`, y=S1_corr_data_ESCA$`% S1 (DS)`)) + geom_point(shape=1)
scat_ESCA + labs(x= "% S1 MP", y="% S1 DS")




