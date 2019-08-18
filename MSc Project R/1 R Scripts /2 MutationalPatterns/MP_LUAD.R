#### EXTRACTING MUT SIG FOR LUAD ####
library(MutationalPatterns)
snvs_toconvert_LUAD <- snvs_LUAD[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_LUAD) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_LUAD$end <- snvs_toconvert_LUAD$start
snvs_toconvert_LUAD$mutation <- apply(snvs_toconvert_LUAD[,c("seqnames","start","REF","ALT")],
                                    1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                          ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_LUAD <- makeGRangesListFromDataFrame(snvs_toconvert_LUAD, split.field = "sample",
                                        keep.extra.columns = TRUE,
                                        names.field = "mutation")

mut_mat_LUAD <- mut_matrix(vcf_list = vcfs_LUAD, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_LUAD <- mut_mat_LUAD + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_LUAD <- nmf(mut_mat_LUAD, rank=2:10, method="brunet", nrun=100, seed=123456) 
plot(estimate_LUAD) # Rank chosen as 5 
nmf_res_LUAD <-extract_signatures(mut_mat_LUAD, rank = 5 , nrun = 100) 

# Compare with COSMIC data
cos_sim_samples_signatures_LUAD = cos_sim_matrix(nmf_res_LUAD$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_LUAD[1,]) # from 1 to XX

# plotting the signature contribution graph 
colnames(nmf_res_LUAD$signatures) <- c("Sig 1", "Sig 4", "Sig 17", "Sig 5", "Sig 2")
rownames(nmf_res_LUAD$contribution) <- c("Sig 1", "Sig 4", "Sig 17", "Sig 5", "Sig 2")
plot_96_profile(nmf_res_LUAD$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_LUAD <- as.data.frame(t(nmf_res_LUAD$contribution))
for (i in 1:567) S1_rel_contr_LUAD$total[i] = sum(S1_rel_contr_LUAD[i,])
for (i in 1:567) S1_rel_contr_LUAD$S1cont[i] = S1_rel_contr_LUAD[i,1]/(S1_rel_contr_LUAD[i,6])
S1_rel_contr_LUAD$Sample <- rownames(S1_rel_contr_LUAD)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_LUAD <- data.frame(matrix("", ncol=1, nrow = 567)) #nrow is [i]
S1_corr_data_LUAD$SampleID <- S1_rel_contr_LUAD$Sample
S1_corr_data_LUAD$matrix.....ncol...1..nrow...567. <- NULL

S1_corr_data_LUAD$'% S1 (DS)' <- sigs.LUAD.dup[match(S1_corr_data_LUAD$SampleID, sigs.LUAD.dup$Sample), "Signature.1"]
S1_corr_data_LUAD$'% S1 (MP)' <- S1_rel_contr_LUAD[match(S1_corr_data_LUAD$SampleID, S1_rel_contr_LUAD$Sample), "S1cont"]

# had to filter:
S1_corr_data_LUAD_filter <- subset(S1_corr_data_LUAD, S1_corr_data_LUAD$`% S1 (MP)` <0.3)
cor.test(S1_corr_data_LUAD_filter$'% S1 (MP)', S1_corr_data_LUAD_filter$'% S1 (DS)')

save(S1_corr_data_LUAD_filter, file="S1.corr.LUAD_filter.RData")

# scatterplot
library(ggplot2)
scat_LUAD <- ggplot(S1_corr_data_LUAD_filter, aes(x=S1_corr_data_LUAD_filter$`% S1 (MP)`, y=S1_corr_data_LUAD_filter$`% S1 (DS)`)) + geom_point(shape=1)
scat_LUAD + labs(x= "% S1 MP", y="% S1 DS")
