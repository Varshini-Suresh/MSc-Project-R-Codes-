#### DECONSTRUCT SIGS ON LUAD ####
library(TCGAbiolinks)
mutations_LUAD <- GDCquery_Maf("LUAD", pipelines = "mutect2")
snvs_LUAD <- data.frame(mutations_LUAD[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_LUAD <- snvs_LUAD[which((snvs_LUAD$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_LUAD$Reference_Allele %in%
                                  c("A","C","G","T"))),] # Filtered to 199580

library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_LUAD <- mut.to.sigs.input(mut.ref = snvs_LUAD, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_LUAD <- NULL
i <- 0
for (sample in rownames(sigs.input_LUAD)) {
  i <- i+1
  print(i)
  sigs_1_LUAD = whichSignatures(tumor.ref = sigs.input_LUAD, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_LUAD <- rbind(sigs_LUAD,sigs_1_LUAD$weights)
}

sigs.LUAD <- data.frame(sigs_LUAD)
save(sigs.LUAD, file="sigs.LUAD.RData")

# duplicating file to have a column on rowname:
sigs.LUAD.dup <- data.frame(sigs_LUAD)
sigs.LUAD.dup$Sample <- rownames(sigs.LUAD.dup)

