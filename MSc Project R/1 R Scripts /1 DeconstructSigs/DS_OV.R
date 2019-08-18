#### DECONSTRUCT SIGS ON OV ####
library(TCGAbiolinks)
mutations_OV <- GDCquery_Maf("OV", pipelines = "mutect2")
snvs_OV <- data.frame(mutations_OV[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_OV <- snvs_OV[which((snvs_OV$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_OV$Reference_Allele %in%
                                  c("A","C","G","T"))),] # filtered from 75168 to 65971 obs

library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_OV <- mut.to.sigs.input(mut.ref = snvs_OV, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_OV <- NULL
i <- 0
for (sample in rownames(sigs.input_OV)) {
  i <- i+1
  print(i)
  sigs_1_OV = whichSignatures(tumor.ref = sigs.input_OV, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_OV <- rbind(sigs_OV,sigs_1_OV$weights)
}

sigs.OV <- data.frame(sigs_OV)
save(sigs.OV, file="sigs.OV.RData")

# duplicating file to have a column on rowname:
sigs.OV.dup <- data.frame(sigs_OV)
sigs.OV.dup$Sample <- rownames(sigs.OV.dup)

