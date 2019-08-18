### DECONSTRUCT SIGS FOR ESCA ###
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(TCGAbiolinks)


#### ESCA ####
mutations <- GDCquery_Maf("ESCA", pipelines = "mutect2") # Downloads data containing all columns
mutations
# Selection of columns
snvs <- data.frame(mutations[,c("Tumor_Sample_Barcode","Hugo_Symbol",
                                "Chromosome","Start_Position","End_Position",
                                "Variant_Classification","Variant_Type",
                                "Reference_Allele","Tumor_Seq_Allele1",
                                "Tumor_Seq_Allele2")])
# Selecting only Point mutations (not indels)
snvs <- snvs[which((snvs$Tumor_Seq_Allele2 %in% 
                      c("A","C","G","T"))&
                     (snvs$Reference_Allele %in%
                        c("A","C","G","T"))),]
View(snvs) 
table(snvs$Tumor_Sample_Barcode)

library(deconstructSigs)
# using mut.to.sigs.input gives mutation freq within each trinucleotide
sigs.input <- mut.to.sigs.input(mut.ref = snvs, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
dim(sigs.input) # 184 tumour samples and 96 mutations (eg 'A[A<G]C')
# Where we number the rows? and link to another file ### Need clarification###
sigs <- NULL
i <- 0
for (sample in rownames(sigs.input)) {
  i <- i+1
  print(i)
  sigs_1 = whichSignatures(tumor.ref = sigs.input, 
                           signatures.ref = signatures.cosmic,
                           sample.id = sample, 
                           contexts.needed = TRUE,
                           signature.cutoff = 0,
                           tri.counts.method = 'default')
  sigs <- rbind(sigs,sigs_1$weights)
}
sigs.ESCA <- data.frame(sigs) # converts list into a data frame 
save(sigs.ESCA, file="sigs.ESCA.RData") # saves the file 

# duplicating file to have a column on rowname:
sigs.KIRC.dup <- data.frame(sigs_KIRC)
sigs.KIRC.dup$Sample <- rownames(sigs.KIRC.dup)
