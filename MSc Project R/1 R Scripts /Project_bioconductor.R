# Installation of Bioconductor Package
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

BiocManager::available() # to find packages; shows the 999 different packages 

# To install BSgenome: 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome", version = "3.8")
browseVignettes("BSgenome") # Opens a web broswer
library(BSgenome)
head(available.genomes())

#### Reference genome #### 
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg38")
library("BSgenome.Hsapiens.UCSC.hg38")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

#### MutationalPatterns - Example data ####
# SECTION 2 # 
# To install a specific package, here, MutationalPatterns
BiocManager::install(c("MutationalPatterns"))
library(MutationalPatterns)
# Locate the VCF files
vcf_files <- list.files(system.file("extdata", package="MutationalPatterns"), 
                        pattern = ".vcf", full.names = TRUE)
# Defining names for the files: 
sample_names <- c("colon1", "colon2", "colon3",
                  "intestine1", "intestine2", "intestine3",
                  "liver1", "liver2", "liver3")
# Loading VCF files as GRangesList (stores genomic range list):
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
summary(vcfs)
# Define metdata such as tissue type on the sample:
tissue <- c(rep("colon", 3), rep("intestine", 3), rep("liver", 3))

## Start of section 3 ##
# Activity on retrieving mutation charcters (maybe this can be done later - but it works so far)
muts = mutations_from_vcf(vcfs[[1]])
head(muts, 12)

#### TCGA ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks", version = "3.8")
browseVignettes("TCGAbiolinks") # Opens a web browser showing documentation
library(TCGAbiolinks)

#### Deconstruct.Sigs ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("deconstructSigs", version = "3.8")
library(deconstructSigs)

#### COSMIC signatures ####
# loading COSMIC data and converting to similar file type as our data
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/", "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33]) # this is now a file containing 30 signatures

save(cancer_signatures, file="COSMIC.RData")

#### SURVMINER - Survival data analysis ####
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/survminer")

