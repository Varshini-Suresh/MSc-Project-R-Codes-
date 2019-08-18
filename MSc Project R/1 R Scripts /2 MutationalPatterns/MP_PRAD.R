#### EXTRACTING MUT SIG FOR PRAD ####
library(MutationalPatterns)
snvs_toconvert_PRAD <- snvs_PRAD[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                                    "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_PRAD) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_PRAD$end <- snvs_toconvert_PRAD$start
snvs_toconvert_PRAD$mutation <- apply(snvs_toconvert_PRAD[,c("seqnames","start","REF","ALT")],
                                      1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                            ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_PRAD <- makeGRangesListFromDataFrame(snvs_toconvert_PRAD, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_PRAD <- mut_matrix(vcf_list = vcfs_PRAD, ref_genome = ref_genome) # this step takes about 5 min
mut_mat_PRAD <- mut_mat_PRAD + 0.0001 # adding pseudogene count 

### Estimating optimal number of signatures - NMF analysis ###
library ("NMF")
estimate_PRAD <- nmf(mut_mat_PRAD, rank=2:10, method="brunet", nrun=100, seed=123456) # took about 5 hours to run 
plot(estimate_PRAD)
nmf_res_PRAD <-extract_signatures(mut_mat_PRAD, rank = 3 , nrun = 100)

# CHOOSE A RANK!! woah looks like 3 will do! 

# compare with COSMIC data 
cos_sim_samples_signatures_PRAD = cos_sim_matrix(nmf_res_PRAD$signatures, cancer_signatures)
which.max(cos_sim_samples_signatures_PRAD[1,]) # from 1 to 5

# plotting the signature contribution graph 
colnames(nmf_res_PRAD$signatures) <- c("Sig 5", "Sig 6", "Sig 1")
rownames(nmf_res_PRAD$contribution) <- c("Sig 5", "Sig 6", "Sig 1")
plot_96_profile(nmf_res_PRAD$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_PRAD <- as.data.frame(t(nmf_res_PRAD$contribution))
for (i in 1:495) S1_rel_contr_PRAD$total[i] = sum(S1_rel_contr_PRAD[i,])
for (i in 1:495) S1_rel_contr_PRAD$S1cont[i] = S1_rel_contr_PRAD[i,3]/(S1_rel_contr_PRAD[i,4])
S1_rel_contr_PRAD$Sample <- rownames(S1_rel_contr_PRAD)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_PRAD <- data.frame(matrix("", ncol=1, nrow = 495)) #nrow is [i]
S1_corr_data_PRAD$SampleID <- S1_rel_contr_PRAD$Sample
S1_corr_data_PRAD$matrix.....ncol...1..nrow...495. <- NULL

S1_corr_data_PRAD$'% S1 (DS)' <- sigs.PRAD.dup[match(S1_corr_data_PRAD$SampleID, sigs.PRAD.dup$Sample), "Signature.1"]
S1_corr_data_PRAD$'% S1 (MP)' <- S1_rel_contr_PRAD[match(S1_corr_data_PRAD$SampleID, S1_rel_contr_PRAD$Sample), "S1cont"]

# HAD TO FILTER: 
S1_corr_data_PRAD_filter <- subset(S1_corr_data_PRAD, S1_corr_data_PRAD$`% S1 (MP)`<0.4)
cor.test(S1_corr_data_PRAD_filter$'% S1 (MP)', S1_corr_data_PRAD_filter$'% S1 (DS)')

save(S1_corr_data_PRAD, file="S1.corr.PRAD.RData")

# scatterplot
library(ggplot2)
scat_PRAD <- ggplot(S1_corr_data_PRAD_filter, aes(x=S1_corr_data_PRAD_filter$`% S1 (MP)`, y=S1_corr_data_PRAD_filter$`% S1 (DS)`)) + geom_point(shape=1)
scat_PRAD + labs(x= "% S1 MP", y="% S1 DS")
