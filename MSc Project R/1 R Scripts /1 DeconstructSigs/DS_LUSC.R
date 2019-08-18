#### DECONSTRUCT SIGS ON LUSC ####
library(TCGAbiolinks)
mutations_LUSC <- GDCquery_Maf("LUSC", pipelines = "mutect2")
snvs_LUSC <- data.frame(mutations_LUSC[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_LUSC <- snvs_LUSC[which((snvs_LUSC$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_LUSC$Reference_Allele %in%
                                  c("A","C","G","T"))),] # 181116 Filtered to 173979 obs

library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_LUSC <- mut.to.sigs.input(mut.ref = snvs_LUSC, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_LUSC <- NULL
i <- 0
for (sample in rownames(sigs.input_LUSC)) {
  i <- i+1
  print(i)
  sigs_1_LUSC = whichSignatures(tumor.ref = sigs.input_LUSC, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_LUSC <- rbind(sigs_LUSC,sigs_1_LUSC$weights)
}

sigs.LUSC <- data.frame(sigs_LUSC)
save(sigs.LUSC, file="sigs.LUSC.RData")

# duplicating file to have a column on rowname:
sigs.LUSC.dup <- data.frame(sigs_LUSC)
sigs.LUSC.dup$Sample <- rownames(sigs.LUSC.dup)

