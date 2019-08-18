#### DECONSTRUCT SIGS ON KIRC ####
library(TCGAbiolinks)
mutations_KIRC <- GDCquery_Maf("KIRC", pipelines = "mutect2")
snvs_KIRC <- data.frame(mutations_KIRC[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_KIRC <- snvs_KIRC[which((snvs_KIRC$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_KIRC$Reference_Allele %in%
                                  c("A","C","G","T"))),] # 26693 filtered to 21833 obs

## Running DS on TCGA dataset ##
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

library(deconstructSigs)
sigs.input_KIRC <- mut.to.sigs.input(mut.ref = snvs_KIRC, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_KIRC <- NULL
i <- 0
for (sample in rownames(sigs.input_KIRC)) {
  i <- i+1
  print(i)
  sigs_1_KIRC = whichSignatures(tumor.ref = sigs.input_KIRC, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_KIRC <- rbind(sigs_KIRC,sigs_1_KIRC$weights)
}

sigs.KIRC <- data.frame(sigs_KIRC)
save(sigs.KIRC, file="sigs.KIRC.RData")

# duplicating file to have a column on rowname:
sigs.KIRC.dup <- data.frame(sigs_KIRC)
sigs.KIRC.dup$Sample <- rownames(sigs.KIRC.dup)

