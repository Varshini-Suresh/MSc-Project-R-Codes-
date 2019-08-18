#### DECONSTRUCT SIGS ON PRAD ####
library(TCGAbiolinks)
mutations_PRAD <- GDCquery_Maf("PRAD", pipelines = "mutect2")
snvs_PRAD <- data.frame(mutations_PRAD[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_PRAD <- snvs_PRAD[which((snvs_PRAD$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_PRAD$Reference_Allele %in%
                                  c("A","C","G","T"))),] # filtered 29286 to 27322 obs

## Running DS on TCGA dataset ##
library(deconstructSigs)
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

sigs.input_PRAD <- mut.to.sigs.input(mut.ref = snvs_PRAD, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_PRAD <- NULL
i <- 0
for (sample in rownames(sigs.input_PRAD)) {
  i <- i+1
  print(i)
  sigs_1_PRAD = whichSignatures(tumor.ref = sigs.input_PRAD, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_PRAD <- rbind(sigs_PRAD,sigs_1_PRAD$weights)
}

sigs.PRAD <- data.frame(sigs_PRAD)
save(sigs.PRAD, file="sigs.PRAD.RData")

# duplicating file to have a column on rowname:
sigs.PRAD.dup <- data.frame(sigs_PRAD)
sigs.PRAD.dup$Sample <- rownames(sigs.PRAD.dup)

