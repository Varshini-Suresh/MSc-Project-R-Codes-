library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library(GenomicRanges)

#### EXTRACTING MUT SIG FOR BRCA ####
## Converting file for BRCA ##
snvs_toconvert_BRCA <- snvs_BRCA[,c("Tumor_Sample_Barcode","Chromosome","Start_Position",
                               "Reference_Allele","Tumor_Seq_Allele2")]

colnames(snvs_toconvert_BRCA) <- c("sample","seqnames","start","REF","ALT")
snvs_toconvert_BRCA$end <- snvs_toconvert_BRCA$start
snvs_toconvert_BRCA$mutation <- apply(snvs_toconvert_BRCA[,c("seqnames","start","REF","ALT")],
                                 1, function(x) paste0(strsplit(x[1],"chr")[[1]][2],
                                                       ":",trimws(x[2]),"_",x[3],"/",x[4]))

vcfs_BRCA <- makeGRangesListFromDataFrame(snvs_toconvert_BRCA, split.field = "sample",
                                          keep.extra.columns = TRUE,
                                          names.field = "mutation")

mut_mat_BRCA <- mut_matrix(vcf_list = vcfs_BRCA, ref_genome = ref_genome)

### NMF Rank for BRCA ###
mut_mat_BRCA <- mut_mat_BRCA + 0.0001
library ("NMF")
estimate_BRCA <- nmf(mut_mat_BRCA, rank=2:10, method="brunet", nrun=100, seed=123456)
plot(estimate_BRCA) # 7 signatures chosen 
nmf_res_BRCA <-extract_signatures(mut_mat_BRCA, rank = 7 , nrun = 100)

# compare with COSMIC data
cosmic_corr_BRCA <- cos_sim_matrix(nmf_res_BRCA$signatures, cancer_signatures)
which.max(cosmic_corr_BRCA[1,]) # from 1 to 7

# plotting the signature contribution graph 
colnames(nmf_res_BRCA$signatures) <- c("Sig 10", "Sig 13", "Sig 2", "Sig 19", "Sig 1", "Sig 3", "Sig 6")
rownames(nmf_res_BRCA$contribution) <- c("Sig 10", "Sig 13", "Sig 2", "Sig 19", "Sig 1", "Sig 3", "Sig 6")
plot_96_profile(nmf_res_BRCA$signatures, condensed = TRUE)

# Calculate relative S1 contribution
S1_rel_contr_BRCA <- as.data.frame(t(nmf_res_BRCA$contribution))
for (i in 1:985) S1_rel_contr_BRCA$total[i] = sum(S1_rel_contr_BRCA[i,])
for (i in 1:985) S1_rel_contr_BRCA$S1cont[i] = S1_rel_contr_BRCA[i,5]/(S1_rel_contr_BRCA[i,8])  
S1_rel_contr_BRCA$Sample <- rownames(S1_rel_contr_BRCA)

#### CORRELATION OF DS AND MP ####
# correlation table 
S1_corr_data_BRCA <- data.frame(matrix("", ncol=1, nrow = 985)) #nrow is [i]
S1_corr_data_BRCA$SampleID <- S1_rel_contr_BRCA$Sample
S1_corr_data_BRCA$matrix.....ncol...1..nrow...985. <- NULL

S1_corr_data_BRCA$'% S1 (DS)' <- sigs.BRCA.dup[match(S1_corr_data_BRCA$SampleID, sigs.BRCA.dup$Sample), "Signature.1"]
S1_corr_data_BRCA$'% S1 (MP)' <- S1_rel_contr_BRCA[match(S1_corr_data_BRCA$SampleID, S1_rel_contr_BRCA$Sample), "S1cont"]

cor.test(S1_corr_data_BRCA$'% S1 (MP)', S1_corr_data_BRCA$'% S1 (DS)')
save(S1_corr_data_BRCA, file="S1.corr.BRCA.RData")

# scatterplot
library(ggplot2)
scat_BRCA <- ggplot(S1_corr_data_BRCA, aes(x=S1_corr_data_BRCA$`% S1 (MP)`, y=S1_corr_data_BRCA$`% S1 (DS)`)) + geom_point(shape=1)
scat_BRCA + labs(x= "% S1 MP", y="% S1 DS")
