#### DECONSTRUCT SIGS ON LIHC ####
library(TCGAbiolinks)
mutations_LIHC <- GDCquery_Maf("LIHC", pipelines = "mutect2")
snvs_LIHC <- data.frame(mutations_LIHC[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                          "Chromosome","Start_Position","End_Position",
                                          "Variant_Classification","Variant_Type",
                                          "Reference_Allele","Tumor_Seq_Allele1",
                                          "Tumor_Seq_Allele2")])

snvs_LIHC <- snvs_LIHC[which((snvs_LIHC$Tumor_Seq_Allele2 %in% 
                                c("A","C","G","T"))&
                               (snvs_LIHC$Reference_Allele %in%
                                  c("A","C","G","T"))),] # 54238 filtered to 50808 obs

## Running DS on TCGA dataset ##
library(deconstructSigs)
sigs.input_LIHC <- mut.to.sigs.input(mut.ref = snvs_LIHC, 
                                     sample.id = "Tumor_Sample_Barcode", 
                                     chr = "Chromosome", 
                                     pos = "Start_Position", 
                                     ref = "Reference_Allele", 
                                     alt = "Tumor_Seq_Allele2",
                                     bsg = BSgenome.Hsapiens.UCSC.hg38) 
# Warning that some samples have fewer than 50 mutations 

sigs_LIHC <- NULL
i <- 0
for (sample in rownames(sigs.input_LIHC)) {
  i <- i+1
  print(i)
  sigs_1_LIHC = whichSignatures(tumor.ref = sigs.input_LIHC, 
                                signatures.ref = signatures.cosmic,
                                sample.id = sample, 
                                contexts.needed = TRUE,
                                signature.cutoff = 0,
                                tri.counts.method = 'default')
  sigs_LIHC <- rbind(sigs_LIHC,sigs_1_LIHC$weights)
}

sigs.LIHC <- data.frame(sigs_LIHC)
save(sigs.LIHC, file="sigs.LIHC.RData")

sigs.default_LIHC <- data.frame(sigs_LIHC)
sigs.default_LIHC$Sample <- rownames(sigs.default_LIHC)

# duplicate file for converting Sample ID into column
sigs.LIHC.dup <- sigs.LIHC
data.table::setDT(sigs.LIHC.dup, keep.rownames = TRUE) [] 
colnames(sigs.LIHC.dup)[colnames(sigs.LIHC.dup)=="rn"] <- "Sample_ID"

