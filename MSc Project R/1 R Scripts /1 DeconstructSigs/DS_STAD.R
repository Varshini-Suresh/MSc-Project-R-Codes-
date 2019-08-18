#### DECONSTRUCT SIGS ON STAD ####
library(TCGAbiolinks)
mutations_STAD <- GDCquery_Maf("STAD", pipelines = "mutect2")
snvs_STAD <- data.frame(mutations_STAD[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_STAD <- snvs_STAD[which((snvs_STAD$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_STAD$Reference_Allele %in%
                                  c("A","C","G","T"))),] # filtered 213144 to 182553 obs

library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_STAD <- mut.to.sigs.input(mut.ref = snvs_STAD, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_STAD <- NULL
i <- 0
for (sample in rownames(sigs.input_STAD)) {
  i <- i+1
  print(i)
  sigs_1_STAD = whichSignatures(tumor.ref = sigs.input_STAD, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_STAD <- rbind(sigs_STAD,sigs_1_STAD$weights)
}

sigs.STAD <- data.frame(sigs_STAD)
save(sigs.STAD, file="sigs.STAD.RData")

# duplicating file to have a column on rowname:
sigs.STAD.dup <- data.frame(sigs_STAD)
sigs.STAD.dup$Sample <- rownames(sigs.STAD.dup)

