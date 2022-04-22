library(tidyr)
library(stringr)
library(eply)
library(plyr)
setwd("/Aging_hallmark_ARD_repository-main")

#################################### ARD to Gene Search terms ############################################# 
# Import and process GWAS ancestry frame
populations <- read.csv("Genetic_data/GWASCat_Categories_update.csv")
populations <- populations[!grepl("founder", populations$Initial),]
populations <- populations[!grepl("founder", populations$Replication),]
populations <- populations[!grepl("Mylopotamos|Pomak|Sirente|Sorbian|Carlantino|Hutterite|Korcula|Split|Cilento|Talana|Vis|Sirente|Sorbian|Talana|Ogliastran|Cilento|Amish|Orcadian|Hutterite|Ashkenazi Jewish|Erasmus Rucphen", populations$Groups),]
populations$Study.Accession <- trimws(populations$Study.Accession)
colnames(populations)[5] <- "STUDY.ACCESSION"

# Import the GWAS catalog
GWASCat <- unzip("Genetic_data/GWAS_Catalog.tsv.zip", exdir = "Genetic_data")
GWASCat <- read.csv("Genetic_data/GWAS_Catalog.tsv", sep="\t")
GWASCat$STUDY.ACCESSION <- trimws(GWASCat$STUDY.ACCESSION)
GWASCat <- merge(GWASCat, populations, by = "STUDY.ACCESSION", all.x = FALSE, all.y=FALSE)

# Processing the ARDs, which are mapped to 'Mapped Traits' from GWAS
MeSH_frame <- read.csv("Genetic_data/MeSH_ARD_to_GWAS.csv")
MeSH_frame$Name_node <- trimws(MeSH_frame$Name_node)
MeSH_frame <- MeSH_frame[!grepl("Infection, Other organisms", MeSH_frame$DiseaseName_ARD_Paper),]
MeSH_frame <- MeSH_frame[!grepl("Infection, Other organs", MeSH_frame$DiseaseName_ARD_Paper),]
MeSH_frame <- MeSH_frame[!grepl("Secondary Malignancy, other", MeSH_frame$DiseaseName_ARD_Paper),]
MeSH_frame <- MeSH_frame[!grepl("Primary Malignancy, other", MeSH_frame$DiseaseName_ARD_Paper),]
paste("The number of unique ARDs are: ", length(unique(MeSH_frame$DiseaseName_ARD_Paper)))
MeSH_frame$MAPPED_TRAIT <- trimws(MeSH_frame$MAPPED_TRAIT)
MeSH_frame$MAPPED_TRAIT_URI <- trimws(MeSH_frame$MAPPED_TRAIT_URI)
Mapped_Trait_frame <- data.frame(MeSH_frame$Name_node, MeSH_frame$DiseaseName_ARD_Paper, MeSH_frame$MAPPED_TRAIT, MeSH_frame$MAPPED_TRAIT_URI)
colnames(Mapped_Trait_frame) <- c("Name_node", "Disease_ARD_Paper", "Mapped_Trait_New", "MAPPED_TRAIT_URI")
Mapped_Trait_frame[Mapped_Trait_frame==""] <- NA
Mapped_Trait_frame <- unique(drop_na(Mapped_Trait_frame))
paste("The number of ARDs included in the GWAS catalog were ", length(unique(Mapped_Trait_frame$Disease_ARD_Paper)))
paste("The number of mapped traits to these ARDs were ", length(unique(Mapped_Trait_frame$Mapped_Trait_New)))
Mapped_Trait_frame$MAPPED_TRAIT_URI <- trimws(Mapped_Trait_frame$MAPPED_TRAIT_URI)
matching_terms <- unique(data.frame(GWASCat$MAPPED_TRAIT, GWASCat$MAPPED_TRAIT_URI))
matching_terms <- unique(matching_terms[!grepl(",", matching_terms$GWASCat.MAPPED_TRAIT),])
matching_terms[matching_terms==""] <- NA
matching_terms <- unique(matching_terms[!is.na(matching_terms$GWASCat.MAPPED_TRAIT),])
write.csv(matching_terms, "matching_terms.csv")

# Processing the GWAS catalog
length(unique(GWASCat$PUBMEDID))
GWASCat <- data.frame(GWASCat$INITIAL.SAMPLE.SIZE, GWASCat$STUDY.ACCESSION, GWASCat$REPLICATION.SAMPLE.SIZE, GWASCat$MAPPED_GENE,
                      GWASCat$UPSTREAM_GENE_ID,GWASCat$DOWNSTREAM_GENE_ID, GWASCat$SNP_GENE_IDS, 
                      GWASCat$UPSTREAM_GENE_DISTANCE, GWASCat$DOWNSTREAM_GENE_DISTANCE, GWASCat$INTERGENIC, 
                      GWASCat$P.VALUE, GWASCat$OR.or.BETA, GWASCat$MAPPED_TRAIT, GWASCat$MAPPED_TRAIT_URI, GWASCat$PUBMEDID)
colnames_GWAS <- gsub(".*GWASCat.", "", colnames(GWASCat))
colnames(GWASCat) <- colnames_GWAS
GWASCat[GWASCat==""] <- NA
GWASCat <- GWASCat[!(is.na(GWASCat$MAPPED_TRAIT)),]
GWASCat <- GWASCat[!(is.na(GWASCat$MAPPED_TRAIT_URI)),]
GWASCat_copy <- GWASCat

# Keeping those significant SNPs within genes mapped to ARDs and significant SNPs in intergenic regions less than 50kbp from a gene
merge_GWAS_MeSH <- merge(GWASCat, Mapped_Trait_frame, on = "MAPPED_TRAIT_URI", all.x = FALSE, all.y=TRUE)
paste("The number of ARDs that are mappable to the GWAS catalog is ", length(unique(merge_GWAS_MeSH$Name_node)))
merge_GWAS_MeSH <- subset(merge_GWAS_MeSH, merge_GWAS_MeSH$P.VALUE<5E-8)
merge_GWAS_MeSH <- subset(merge_GWAS_MeSH, merge_GWAS_MeSH$UPSTREAM_GENE_DISTANCE<50000|merge_GWAS_MeSH$DOWNSTREAM_GENE_DISTANCE<50000|merge_GWAS_MeSH$INTERGENIC==0)
merge_GWAS_MeSH$SNP_GENE_IDS <- as.character(merge_GWAS_MeSH$SNP_GENE_IDS)
merge_GWAS_MeSH <- separate_rows(merge_GWAS_MeSH, sep=",", SNP_GENE_IDS)
merge_GWAS_MeSH$SNP_GENE_IDS <- trimws(merge_GWAS_MeSH$SNP_GENE_IDS)
paste("The number of ARDs represented in the data kept from the GWAS catalog was ", length(unique(merge_GWAS_MeSH$Disease_ARD_Paper)))
paste("The number of mapped traits to these ARDs were ", length(unique(merge_GWAS_MeSH$Mapped_Trait_New)))
paste("The number of studies represented were ", length(unique(merge_GWAS_MeSH$STUDY.ACCESSION)))
ARD_Gene_frame <- merge_GWAS_MeSH
ARD_Gene_frame <- data.frame(ARD_Gene_frame$Name_node, ARD_Gene_frame$MAPPED_TRAIT, ARD_Gene_frame$UPSTREAM_GENE_ID, ARD_Gene_frame$DOWNSTREAM_GENE_ID, 
                             ARD_Gene_frame$SNP_GENE_IDS, ARD_Gene_frame$INTERGENIC, ARD_Gene_frame$UPSTREAM_GENE_DISTANCE, ARD_Gene_frame$DOWNSTREAM_GENE_DISTANCE, ARD_Gene_frame$PUBMEDID)
ARD_Gene_frame_not_intergenic <- unique(subset(ARD_Gene_frame, ARD_Gene_frame$ARD_Gene_frame.INTERGENIC == 0))
ARD_Gene_frame_not_intergenic <- data.frame(ARD_Gene_frame_not_intergenic$ARD_Gene_frame.Name_node, ARD_Gene_frame_not_intergenic$ARD_Gene_frame.SNP_GENE_IDS, ARD_Gene_frame_not_intergenic$ARD_Gene_frame.PUBMEDID, ARD_Gene_frame_not_intergenic$ARD_Gene_frame.MAPPED_TRAIT)
colnames(ARD_Gene_frame_not_intergenic) <-  c("Name_node", "SNP_Gene_IDS", "PUBMEDID", "MAPPED_TRAIT")
ARD_Gene_frame_intergenic <- subset(ARD_Gene_frame, ARD_Gene_frame$ARD_Gene_frame.DOWNSTREAM_GENE_DISTANCE<50000)
ARD_Gene_frame_intergenic <- unique(data.frame(ARD_Gene_frame_intergenic$ARD_Gene_frame.Name_node, ARD_Gene_frame_intergenic$ARD_Gene_frame.DOWNSTREAM_GENE_ID, ARD_Gene_frame_intergenic$ARD_Gene_frame.PUBMEDID, ARD_Gene_frame_intergenic$ARD_Gene_frame.MAPPED_TRAIT))
colnames(ARD_Gene_frame_intergenic) <-  c("Name_node", "SNP_Gene_IDS", "PUBMEDID", "MAPPED_TRAIT")
ARD_Gene_frame_intergenic2 <- subset(ARD_Gene_frame, ARD_Gene_frame$ARD_Gene_frame.UPSTREAM_GENE_DISTANCE<50000)
ARD_Gene_frame_intergenic2 <- unique(data.frame(ARD_Gene_frame_intergenic2$ARD_Gene_frame.Name_node, ARD_Gene_frame_intergenic2$ARD_Gene_frame.UPSTREAM_GENE_ID, ARD_Gene_frame_intergenic2$ARD_Gene_frame.PUBMEDID, ARD_Gene_frame_intergenic2$ARD_Gene_frame.MAPPED_TRAIT))
colnames(ARD_Gene_frame_intergenic2) <-  c("Name_node", "SNP_Gene_IDS", "PUBMEDID", "MAPPED_TRAIT")
ARD_Gene_frame_all <- unique(rbind(ARD_Gene_frame_not_intergenic, ARD_Gene_frame_intergenic, ARD_Gene_frame_intergenic2))
ARD_Gene_frame_all$Name_node <- trimws(ARD_Gene_frame_all$Name_node)
ARD_Gene_frame_all$SNP_Gene_IDS <- trimws(ARD_Gene_frame_all$SNP_Gene_IDS)
ARD_Gene_frame_all <- separate_rows(ARD_Gene_frame_all, sep=",", SNP_Gene_IDS)
paste("The number of ARDs represented in the data kept from the GWAS catalog was ", length(unique(ARD_Gene_frame_all$Name_node)))
paste("The number of SNPs represented were ", length(unique(ARD_Gene_frame_all$SNP_Gene_IDS)))
paste("The number of PubMed IDs represented were", length(unique(ARD_Gene_frame_all$PUBMEDID)))
ARD_Gene_frame_all_copy <- ARD_Gene_frame_all
ARD_Gene_frame_all$PUBMEDID <- NULL
ARD_Gene_frame_all <- unique(ARD_Gene_frame_all)
colnames(ARD_Gene_frame_all) <- c("name_node", "dbXrefs", "mapped_trait")
ARD_Gene_frame_all$name_node <- trimws(ARD_Gene_frame_all$name_node)

# Assign NCBI Gene ID to genes
homo_sapiens <- read.csv("Genetic_data/Homo_sapiens_gene_info.txt", sep ="\t")
homo_sapiens[] <- lapply(homo_sapiens, gsub, pattern = "|", replacement = ";", fixed = TRUE)
homo_sapiens <- separate_rows(homo_sapiens, sep=";", dbXrefs)
homo_sapiens$dbXrefs <- gsub(".*ENSG", "ENSG", homo_sapiens$dbXrefs)
homo_sapiens <- homo_sapiens[grepl("ENSG", homo_sapiens$dbXrefs),]
ARD_Gene_frame_all <- merge(homo_sapiens, ARD_Gene_frame_all, by = "dbXrefs", all.x = FALSE, all.y=TRUE)
ARD_Gene_frame_all <- data.frame(ARD_Gene_frame_all$dbXrefs, ARD_Gene_frame_all$name_node, ARD_Gene_frame_all$GeneID, ARD_Gene_frame_all$Symbol, ARD_Gene_frame_all$mapped_trait)
colnames(ARD_Gene_frame_all) <- c("dbXrefs", "name_node", "GeneID", "Symbol", "mapped_trait")
paste("The number of Ensembl gene IDs represented were ", length(unique(ARD_Gene_frame_all$dbXrefs)))
paste("The number of NCBI codes represented were ", length(unique(ARD_Gene_frame_all$GeneID)))
paste("The number of NCBI gene symbols represented were ", length(unique(ARD_Gene_frame_all$Symbol)))
ARD_Gene_frame_all <- ARD_Gene_frame_all[!is.na(ARD_Gene_frame_all$GeneID), ]
paste("The number of ARDs with NCBI codes were ", length(unique(ARD_Gene_frame_all$name_node)))
paste("The number of NCBI codes were ", length(unique(ARD_Gene_frame_all$GeneID)))
paste("The number of mapped traits were ", length(unique(ARD_Gene_frame_all$mapped_trait)))
paste("The number of Ensembl gene IDs mapped were ", length(unique(ARD_Gene_frame_all$dbXrefs)))
write.csv(ARD_Gene_frame_all, "Genetic_data/ARD_Gene_frame_all_NCBI_50000.csv")

# For network propagation - Keeping those significant SNPs within genes mapped to ARDs and significant SNPs in intergenic regions less than 75kbp from a gene
merge_GWAS_MeSH_lower <- merge(GWASCat_copy, Mapped_Trait_frame, on = "MAPPED_TRAIT_URI", all.x = FALSE, all.y=TRUE)
merge_GWAS_MeSH_lower <- subset(merge_GWAS_MeSH_lower, merge_GWAS_MeSH_lower$P.VALUE<5E-8)
merge_GWAS_MeSH_lower <- subset(merge_GWAS_MeSH_lower, merge_GWAS_MeSH_lower$UPSTREAM_GENE_DISTANCE<75000|merge_GWAS_MeSH_lower$DOWNSTREAM_GENE_DISTANCE<75000|merge_GWAS_MeSH_lower$INTERGENIC==0)
merge_GWAS_MeSH_lower$SNP_GENE_IDS <- as.character(merge_GWAS_MeSH_lower$SNP_GENE_IDS)
merge_GWAS_MeSH_lower <- separate_rows(merge_GWAS_MeSH_lower, sep=",", SNP_GENE_IDS)
merge_GWAS_MeSH_lower$SNP_GENE_IDS <- trimws(merge_GWAS_MeSH_lower$SNP_GENE_IDS)
ARD_Gene_frame_lower <- data.frame(merge_GWAS_MeSH_lower$Name_node, merge_GWAS_MeSH_lower$MAPPED_TRAIT, merge_GWAS_MeSH_lower$UPSTREAM_GENE_ID, merge_GWAS_MeSH_lower$DOWNSTREAM_GENE_ID, 
                                   merge_GWAS_MeSH_lower$SNP_GENE_IDS, merge_GWAS_MeSH_lower$INTERGENIC, merge_GWAS_MeSH_lower$UPSTREAM_GENE_DISTANCE, merge_GWAS_MeSH_lower$DOWNSTREAM_GENE_DISTANCE,
                                   merge_GWAS_MeSH_lower$PUBMEDID)
colnames(ARD_Gene_frame_lower) <- gsub(".*merge_GWAS_MeSH_lower.", "", colnames(ARD_Gene_frame_lower))
ARD_Gene_frame_lower_not_intergenic <- unique(subset(ARD_Gene_frame_lower, ARD_Gene_frame_lower$INTERGENIC ==0))
ARD_Gene_frame_lower_not_intergenic <- data.frame(ARD_Gene_frame_lower_not_intergenic$Name_node, ARD_Gene_frame_lower_not_intergenic$SNP_GENE_IDS, ARD_Gene_frame_lower_not_intergenic$PUBMEDID)
colnames(ARD_Gene_frame_lower_not_intergenic) <-  c("Name_node", "SNP_GENE_IDS", "PUBMEDID")
ARD_Gene_frame_lower_intergenic <- subset(ARD_Gene_frame_lower, ARD_Gene_frame_lower$DOWNSTREAM_GENE_DISTANCE<75000)
ARD_Gene_frame_lower_intergenic <- unique(data.frame(ARD_Gene_frame_lower_intergenic$Name_node, ARD_Gene_frame_lower_intergenic$DOWNSTREAM_GENE_ID, ARD_Gene_frame_lower_intergenic$PUBMEDID))
colnames(ARD_Gene_frame_lower_intergenic) <-  c("Name_node", "SNP_GENE_IDS", "PUBMEDID")
ARD_Gene_frame_lower_intergenic2 <- subset(ARD_Gene_frame_lower, ARD_Gene_frame_lower$UPSTREAM_GENE_DISTANCE<75000)
ARD_Gene_frame_lower_intergenic2 <- unique(data.frame(ARD_Gene_frame_lower_intergenic2$Name_node, ARD_Gene_frame_lower_intergenic2$UPSTREAM_GENE_ID, ARD_Gene_frame_lower_intergenic2$PUBMEDID))
colnames(ARD_Gene_frame_lower_intergenic2) <-  c("Name_node", "SNP_GENE_IDS", "PUBMEDID")
ARD_Gene_frame_lower_all <- unique(rbind(ARD_Gene_frame_lower_not_intergenic, ARD_Gene_frame_lower_intergenic, ARD_Gene_frame_lower_intergenic2))
ARD_Gene_frame_lower_all$Name_node <- trimws(ARD_Gene_frame_lower_all$Name_node)
ARD_Gene_frame_lower_all$SNP_GENE_IDS <- trimws(ARD_Gene_frame_lower_all$SNP_GENE_IDS)
write.csv(ARD_Gene_frame_lower_all$PUBMEDID, "data/human_aging_corpus/PubMedID_list.csv")
setdiff(ARD_Gene_frame_lower_all$PUBMEDID, ARD_Gene_frame_all_copy$PUBMEDID)
colnames(ARD_Gene_frame_lower_all) <- c("name_node", "dbXrefs", "PUBMEDID")
ARD_Gene_frame_lower_all <- merge(homo_sapiens,ARD_Gene_frame_lower_all, by = "dbXrefs", all.x = FALSE, all.y=TRUE)
ARD_Gene_frame_lower_all <- unique(data.frame(ARD_Gene_frame_lower_all$dbXrefs, ARD_Gene_frame_lower_all$name_node, ARD_Gene_frame_lower_all$GeneID, ARD_Gene_frame_lower_all$Symbol))
colnames(ARD_Gene_frame_lower_all) <- c("dbXrefs", "name_node", "GeneID", "Symbol")
write.csv(ARD_Gene_frame_lower_all, "Genetic_data/ARD_Gene_frame_all_lower_NCBI_50000.csv")

