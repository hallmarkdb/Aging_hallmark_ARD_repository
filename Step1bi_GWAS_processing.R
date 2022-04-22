library(tidyr)
library(stringr)
library(eply)
library(plyr)

#################################### ARD to Gene Search terms ############################################# 
# Processing the GWAS Ancestry Frame and keeping those with specific broad ancestral categories
setwd("/Code_Thesis_HCF/Aging_hallmark_ARD_repository")
GWAS_Ancestry <- read.csv("Genetic_data/Ancestry.tsv", sep="\t", header=TRUE)
GWAS_Ancestry_Europe <- GWAS_Ancestry[grepl("Euro", GWAS_Ancestry$BROAD.ANCESTRAL.CATEGORY),]
GWAS_Ancestry_Europe1 <- data.frame(GWAS_Ancestry_Europe$INITIAL.SAMPLE.DESCRIPTION,GWAS_Ancestry_Europe$REPLICATION.SAMPLE.DESCRIPTION, GWAS_Ancestry_Europe$STUDY.ACCESSION, GWAS_Ancestry_Europe$BROAD.ANCESTRAL.CATEGORY)
GWAS_Ancestry_Europe1$both <- paste(GWAS_Ancestry_Europe$INITIAL.SAMPLE.DESCRIPTION,GWAS_Ancestry_Europe$REPLICATION.SAMPLE.DESCRIPTION, sep=" ; ")
GWAS_Ancestry_Europe1$both <- gsub("[0-9]", "", GWAS_Ancestry_Europe1$both)
GWAS_Ancestry_Europe1$both <- gsub(", ", " ", GWAS_Ancestry_Europe1$both)
GWAS_Ancestry_Europe1$both <- gsub("Up", "", GWAS_Ancestry_Europe1$both)
GWAS_Ancestry_Europe1$both <- gsub("; NA", " ", GWAS_Ancestry_Europe1$both)
GWAS_Ancestry_Europe1$both <- trimws(GWAS_Ancestry_Europe1$both)
GWAS_Ancestry_Europe1 <- GWAS_Ancestry_Europe1[grepl("[A-Z]", GWAS_Ancestry_Europe1$both),]
colnames(GWAS_Ancestry_Europe1) <- c("Initial", "Replication", "Accesion", "Broad", "both")
write.csv(GWAS_Ancestry_Europe1, "Genetic_data/GWAS_Ancestry_Europe.csv")