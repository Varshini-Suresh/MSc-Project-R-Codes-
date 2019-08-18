#### DECONSTRUCT SIGS ON BLCA ####
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Retrieving relevant columns from TCGA ##
library(TCGAbiolinks)
mutations_BLCA <- GDCquery_Maf("BLCA", pipelines = "mutect2")
snvs_BLCA <- data.frame(mutations_BLCA[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_BLCA <- snvs_BLCA[which((snvs_BLCA$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_BLCA$Reference_Allele %in%
                                  c("A","C","G","T"))),] # filtered to 130255 obs

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_BLCA <- mut.to.sigs.input(mut.ref = snvs_BLCA, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_BLCA <- NULL
i <- 0
for (sample in rownames(sigs.input_BLCA)) {
  i <- i+1
  print(i)
  sigs_1_BLCA = whichSignatures(tumor.ref = sigs.input_BLCA, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_BLCA <- rbind(sigs_BLCA,sigs_1_BLCA$weights)
}

sigs.BLCA <- data.frame(sigs_BLCA)
save(sigs.BLCA, file="sigs.BLCA.RData")

# duplicating file to have a column on rowname:
sigs.BLCA.dup <- data.frame(sigs_BLCA)
sigs.BLCA.dup$Sample <- rownames(sigs.BLCA.dup)




