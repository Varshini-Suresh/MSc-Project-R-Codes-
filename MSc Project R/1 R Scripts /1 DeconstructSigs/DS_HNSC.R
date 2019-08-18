#### DECONSTRUCT SIGS ON HNSC ####
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Retrieving relevant columns from TCGA ##
library(TCGAbiolinks)
mutations_HNSC <- GDCquery_Maf("HNSC", pipelines = "mutect2")
snvs_HNSC <- data.frame(mutations_HNSC[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_HNSC <- snvs_HNSC[which((snvs_HNSC$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_HNSC$Reference_Allele %in%
                                  c("A","C","G","T"))),] # filtered from 102309 to 96932 obs

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_HNSC <- mut.to.sigs.input(mut.ref = snvs_HNSC, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_HNSC <- NULL
i <- 0
for (sample in rownames(sigs.input_HNSC)) {
  i <- i+1
  print(i)
  sigs_1_HNSC = whichSignatures(tumor.ref = sigs.input_HNSC, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_HNSC <- rbind(sigs_HNSC,sigs_1_HNSC$weights)
}

sigs.HNSC <- data.frame(sigs_HNSC)
save(sigs.HNSC, file="sigs.HNSC.RData")

# duplicating file to have a column on rowname:
sigs.HNSC.dup <- data.frame(sigs_HNSC)
sigs.HNSC.dup$Sample <- rownames(sigs.HNSC.dup)

