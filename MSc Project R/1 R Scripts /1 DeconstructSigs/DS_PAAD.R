#### DECONSTRUCT SIGS ON PAAD ####
library(TCGAbiolinks)
mutations_PAAD <- GDCquery_Maf("PAAD", pipelines = "mutect2")
snvs_PAAD <- data.frame(mutations_PAAD[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                      "Chromosome","Start_Position","End_Position",
                                      "Variant_Classification","Variant_Type",
                                      "Reference_Allele","Tumor_Seq_Allele1",
                                      "Tumor_Seq_Allele2")])

snvs_PAAD <- snvs_PAAD[which((snvs_PAAD$Tumor_Seq_Allele2 %in% 
                            c("A","C","G","T"))&
                           (snvs_PAAD$Reference_Allele %in%
                              c("A","C","G","T"))),] # filtered 29959 to 29267 obs

## Running DS on TCGA dataset ##
library(deconstructSigs)
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
sigs.input_PAAD <- mut.to.sigs.input(mut.ref = snvs_PAAD, 
                                   sample.id = "Tumor_Sample_Barcode", 
                                   chr = "Chromosome", 
                                   pos = "Start_Position", 
                                   ref = "Reference_Allele", 
                                   alt = "Tumor_Seq_Allele2",
                                   bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_PAAD <- NULL
i <- 0
for (sample in rownames(sigs.input_PAAD)) {
  i <- i+1
  print(i)
  sigs_1_PAAD = whichSignatures(tumor.ref = sigs.input_PAAD, 
                              signatures.ref = signatures.cosmic,
                              sample.id = sample, 
                              contexts.needed = TRUE,
                              signature.cutoff = 0,
                              tri.counts.method = 'default')
  sigs_PAAD <- rbind(sigs_PAAD,sigs_1_PAAD$weights)
}

sigs.PAAD <- data.frame(sigs_PAAD)
save(sigs.PAAD, file="sigs.PAAD.RData")

# duplicating file to have a column on rowname:
sigs.PAAD.dup <- data.frame(sigs_PAAD)
sigs.PAAD.dup$Sample <- rownames(sigs.PAAD.dup)

