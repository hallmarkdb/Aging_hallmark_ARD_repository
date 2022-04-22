library(tidyr)
library(stringr)
library(eply)
library(plyr)

################################### Development of ARD dictionaries ###################################
# Importing data
setwd("/Code_Thesis_HCF/Aging_hallmark_ARD_repository")
MeSH_Terms <- unique(read.csv("data/ARD_dictionary/MeSH_frames_original.csv", sep =","))
MeSH_Terms <- unique(data.frame(MeSH_Terms$DiseaseName_ARD_Paper, MeSH_Terms$MESH_Term, MeSH_Terms$MESH_ID))
colnames(MeSH_Terms) <- c("DiseaseName_ARD_Paper", "MESH_Term", "DiseaseID")
MeSH_Terms$DiseaseID <- trimws(MeSH_Terms$DiseaseID)
MeSH_Terms <- MeSH_Terms[grepl("MESH:D", MeSH_Terms$DiseaseID), ]
paste("Of 207 ARDs, the number of ARDs mapped to MESH:D terms was", length(unique(MeSH_Terms$DiseaseName_ARD_Paper)))
MeSH_Terms <- MeSH_Terms[!grepl("Infection, Other organisms", MeSH_Terms$DiseaseName_ARD_Paper),]
MeSH_Terms <- MeSH_Terms[!grepl("Infection, Other organs", MeSH_Terms$DiseaseName_ARD_Paper),]
MeSH_Terms <- MeSH_Terms[!grepl("Secondary Malignancy, other", MeSH_Terms$DiseaseName_ARD_Paper),]
MeSH_Terms <- MeSH_Terms[!grepl("Primary Malignancy, other", MeSH_Terms$DiseaseName_ARD_Paper),]
paste("Of the 203 ARDs included, the number of ARDs mapped to MESH:D terms was", length(unique(MeSH_Terms$DiseaseName_ARD_Paper)))

# Creating a hierarchical tree of MeSH terms
CTD_diseases <- unique(read.csv("data/ARD_dictionary/CTD_diseases_New.csv", sep =",", header=TRUE))
CTD_diseases <- CTD_diseases[grepl("MESH:D", CTD_diseases$DiseaseID), ]
CTD_diseases$SlimMappings <- tolower(CTD_diseases$SlimMappings)
CTD_diseases <- CTD_diseases[!grepl("animal disease", CTD_diseases$SlimMappings), ]
CTD_diseases$SlimMappings <- NULL
CTD_diseases <- data.frame(CTD_diseases$ParentIDs, CTD_diseases$DiseaseID, CTD_diseases$DiseaseName)
colnames(CTD_diseases) <- c("ParentIDs", "DiseaseID", "DiseaseName")
CTD_diseases[] <- lapply(CTD_diseases, gsub, pattern = "|", replacement = ";", fixed = TRUE)
CTD_diseases <- separate_rows(CTD_diseases, sep=";", ParentIDs)
CTD_diseases[CTD_diseases==" "] <- NA
CTD_diseases$DiseaseName <- trimws(CTD_diseases$DiseaseName)
CTD_diseases$DiseaseID <- trimws(CTD_diseases$DiseaseID)
CTD_diseases$ParentIDs <- trimws(CTD_diseases$ParentIDs)
CTD_diseases[CTD_diseases==""] <- NA
CTD_diseases <- unique(CTD_diseases)
paste("The file contains", dim(CTD_diseases)[1], "rows")

CTDdisease_1 <- data.frame(CTD_diseases$DiseaseID, CTD_diseases$DiseaseName)
colnames(CTDdisease_1) <- c("DiseaseID", "DiseaseName_CTD")
CTDdisease_1 <- unique(CTDdisease_1)
paste("The number of diseases mapped to MESH:D terms in the CTD database is", length(unique(CTDdisease_1$DiseaseName_CTD)))
New_ARD_Tables <- merge(MeSH_Terms, CTDdisease_1, by = "DiseaseID", all.x =TRUE, all.y=FALSE)
New_ARD_Tables[New_ARD_Tables==""] <- NA

IDsinIDs <- CTD_diseases[CTD_diseases$ParentIDs %in% New_ARD_Tables$DiseaseID,]
IDsinIDs2 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs$DiseaseID, ]
IDsinIDs3 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs2$DiseaseID, ]
IDsinIDs4 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs3$DiseaseID, ]
IDsinIDs5 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs4$DiseaseID, ]
IDsinIDs6 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs5$DiseaseID, ]
IDsinIDs7 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs6$DiseaseID, ]
IDsinIDs8 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs7$DiseaseID, ]
IDsinIDs9 <- CTD_diseases[CTD_diseases$ParentIDs %in% IDsinIDs8$DiseaseID, ]

colnames(IDsinIDs) <- c("DiseaseID", "DiseaseID2", "DiseaseName2")
colnames(IDsinIDs2) <- c("DiseaseID2", "DiseaseID3", "DiseaseName3")
colnames(IDsinIDs3) <- c("DiseaseID3", "DiseaseID4", "DiseaseName4")
colnames(IDsinIDs4) <- c("DiseaseID4", "DiseaseID5", "DiseaseName5")
colnames(IDsinIDs5) <- c("DiseaseID5", "DiseaseID6", "DiseaseName6")
colnames(IDsinIDs6) <- c("DiseaseID6", "DiseaseID7", "DiseaseName7")
colnames(IDsinIDs7) <- c("DiseaseID7", "DiseaseID8", "DiseaseName8")
colnames(IDsinIDs8) <- c("DiseaseID8", "DiseaseID9", "DiseaseName9")

New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs, by="DiseaseID", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs2, by="DiseaseID2", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs3, by="DiseaseID3", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs4, by="DiseaseID4", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs5, by="DiseaseID5", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs6, by="DiseaseID6", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs7, by="DiseaseID7", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- merge(New_ARD_Tables, IDsinIDs8, by="DiseaseID8", all.x=TRUE, all.y=TRUE)
New_ARD_Tables <- unique(New_ARD_Tables)
New_ARD_Tables <- New_ARD_Tables[, c(19, 1:18, 20)]
paste("The exploded MeSH terms has", dim(New_ARD_Tables)[1], "rows")
New_ARD_Tables <- New_ARD_Tables[order(New_ARD_Tables$DiseaseName_ARD_Paper),]
write.csv(New_ARD_Tables, "data/ARD_dictionary/MeSH_tree_table.csv")

# Using the hierarchical tree of MeSH terms, we form a subclasses of the original MeSH terms
MeSH_Dataframe_final <- unique(read.csv("data/ARD_dictionary/MeSH_frames.csv", sep =","))
MeSH_Dataframe_final <- MeSH_Dataframe_final[!grepl("Infection, Other organisms", MeSH_Dataframe_final$DiseaseName_ARD_Paper),]
MeSH_Dataframe_final <- MeSH_Dataframe_final[!grepl("Infection, Other organs", MeSH_Dataframe_final$DiseaseName_ARD_Paper),]
MeSH_Dataframe_final <- MeSH_Dataframe_final[!grepl("Secondary Malignancy, other", MeSH_Dataframe_final$DiseaseName_ARD_Paper),]
MeSH_Dataframe_final <- MeSH_Dataframe_final[!grepl("Primary Malignancy, other", MeSH_Dataframe_final$DiseaseName_ARD_Paper),]
MeSH_Dataframe_final <- MeSH_Dataframe_final[grepl("MESH:D", MeSH_Dataframe_final$All_MESH),]
paste("We mapped ", length(unique(MeSH_Dataframe_final$All_MESH_Terms)), "subclasses to the ARDs")
paste("The subclasses were mapped to ", length(unique(MeSH_Dataframe_final$DiseaseName_ARD_Paper)), "ARDs")

# We then derive synonyms to ARDs for the ARD dictionary
ctd_database <- read.csv("data/ARD_dictionary/CTD_diseases_New.csv", sep =",", header=TRUE)
ctd_database <- data.frame(ctd_database$DiseaseName, ctd_database$DiseaseID, ctd_database$AltDiseaseIDs, ctd_database$Synonyms, ctd_database$SlimMappings)
colnames(ctd_database) <- c("DiseaseName", "DiseaseID", "AltDiseaseIDs", "Synonyms", "SlimMappings")
ctd_database$DiseaseName <- trimws(ctd_database$DiseaseName)
ctd_database$DiseaseID <- trimws(ctd_database$DiseaseID)
ctd_database$AltDiseaseIDs <- trimws(ctd_database$AltDiseaseIDs)
ctd_database$Synonyms <- trimws(ctd_database$Synonyms)
ctd_database$SlimMappings <- trimws(ctd_database$SlimMappings)
ctd_database[ctd_database==""] <- NA
ctd_database <- unique(ctd_database)
paste("The number of unique diseases is", length(unique(ctd_database$DiseaseName))) 

# Identify rows where the “DiseaseID” contains the “MESH:D” term, remove terms assigned to animal diseases
ctd_database <- ctd_database[grepl("MESH:D", ctd_database$DiseaseID), ]
paste("The number of diseases with a MeSH:D term is", length(unique(ctd_database$DiseaseName)))
ctd_Alt_MESH <- ctd_database[grep("MESH", ctd_database$AltDiseaseID), ]
paste("The number of disease IDs where the Alt Disease ID contains MESH is", length(unique(ctd_Alt_MESH$DiseaseName)))
ctd_database$AltDiseaseIDs <- NULL
rm(ctd_Alt_MESH)
ctd_database$SlimMappings <- tolower(ctd_database$SlimMappings)
ctd_database <- ctd_database[!grepl("animal disease", ctd_database$SlimMappings), ]
ctd_database$SlimMappings <- NULL
paste("The number of unique diseases with a MESH:D term that are not animal diseases is", length(unique(ctd_database$DiseaseName)))

# Processing synonyms: Separate the synonyms into different rows, use "DiseaseName" column to create synonyms, and clean synonyms
ctd_database[] <- lapply(ctd_database, gsub, pattern = "|", replacement = ";", fixed = TRUE)
ctd_database <- separate_rows(ctd_database, sep=";", Synonyms)
ctd_database$Synonyms <- trimws(ctd_database$Synonyms)
colnames(ctd_database) <- c("DiseaseName", "DiseaseID", "Synonyms") 
paste("After separating synonyms there were", length(unique(ctd_database$Synonyms)), "synonyms")

synonyms_extra <- ctd_database[,c("DiseaseName", "DiseaseID", "DiseaseName")]
colnames(synonyms_extra) <- c("DiseaseName", "DiseaseID", "Synonyms")
synonyms_extra <- unique(synonyms_extra)
ctd_database <- unique(rbind(ctd_database, synonyms_extra))
ctd_database[ctd_database==""] <- NA
ctd_database <- ctd_database[!is.na(ctd_database$Synonyms),]
colnames(ctd_database)[3] <- c("AllSynonyms")
paste("There are", length(unique(ctd_database$AllSynonyms)), "unique synonyms after adding disease names termed 'All Synonyms'")

ctd_database[] <- lapply(ctd_database, str_squish)
ctd_database$AllSynonyms <- trimws(ctd_database$AllSynonyms)
ctd_database$AllSynonyms <- as.character(ctd_database$AllSynonyms)
ctd_database$AllSynonyms <- unquote(ctd_database$AllSynonyms)
ctd_database[] <- lapply(ctd_database, gsub, pattern = "-", replacement = " ", fixed = TRUE)
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, "\\(.*\\)", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, "'", "")
ctd_database <- unique(ctd_database)
paste("There are", length(unique(ctd_database$AllSynonyms)), "unique synonyms after cleaning")

# Keep words without full stop, but keep all decimals in numbers
ctd_remainder <- ctd_database[!grepl("[a-zA-Z]\\.", ctd_database$AllSynonyms),]
ctd_subset_dot <- ctd_database[grepl("[a-zA-Z]\\.", ctd_database$AllSynonyms),]
ctd_subset_dot$AllSynonyms <- str_replace(ctd_subset_dot$AllSynonyms, "\\.", "")
ctd_subset_dot$AllSynonyms <- str_replace(ctd_subset_dot$AllSynonyms, "\\.", "")
ctd_database <- unique(rbind(ctd_subset_dot, ctd_remainder))
rm(ctd_subset_dot, ctd_remainder)
paste("There are", length(unique(ctd_database$AllSynonyms)), "unique synonyms after examining letters followed by decimal/ full stop etc.")

# Remove the terms "included", "formerly, "susceptibility to" under certain conditions
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", INCLUDED", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", Included", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", included", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", FORMERLY", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", Formerly", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", formerly", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", SUSCEPTIBILITY TO, [0-9]", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", Susceptibility To, [0-9]", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", susceptibility to, [0-9]", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", SUSCEPTIBILITY TO$", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", Susceptibility To$", "")
ctd_database$AllSynonyms <- str_replace(ctd_database$AllSynonyms, ", susceptibility to$", "")
ctd_database <- unique(ctd_database)
paste("There are", length(unique(ctd_database$AllSynonyms)), "unique synonyms after removing terms with 'included, formerly, susceptibility to' under certain conditions")

# Examine & process synonyms with commas
ctd_subset_nocommas <- ctd_database[!grepl(",", ctd_database$AllSynonyms),]
ctd_subset_commas <- ctd_database[grepl(",", ctd_database$AllSynonyms),]
ctd_subset_commas$count <- str_count(ctd_subset_commas$AllSynonyms, ",")
ctd_subset_commas <- subset(ctd_subset_commas, ctd_subset_commas$count == 1)

ctd_subset_commas_09_first <- ctd_subset_commas[grepl("[0-9],", ctd_subset_commas$AllSynonyms),]
ctd_subset_commas_09_second <- ctd_subset_commas[grepl("[a-zA-Z], [0-9]", ctd_subset_commas$AllSynonyms),]
ctd_subset_commas_09 <- rbind(ctd_subset_commas_09_first, ctd_subset_commas_09_second)

ctd_subset_commas_09_nosep1 <-  ctd_subset_commas_09[grepl("[0-9],[0-9]", ctd_subset_commas_09$AllSynonyms),]
ctd_subset_commas_09X <-  ctd_subset_commas_09[grepl("[0-9],[XX]", ctd_subset_commas_09$AllSynonyms),] 
ctd_subset_commas_09X2 <-  ctd_subset_commas_09[grepl("[0-9], [XX]", ctd_subset_commas_09$AllSynonyms),]
ctd_subset_commas_09XY <-  ctd_subset_commas_09[grepl("[0-9],[XY]", ctd_subset_commas_09$AllSynonyms),] 
ctd_subset_commas_09XY2 <-  ctd_subset_commas_09[grepl("[0-9], [XY]", ctd_subset_commas_09$AllSynonyms),]
numerics_commas <- rbind(ctd_subset_commas_09_nosep1, ctd_subset_commas_09X, ctd_subset_commas_09X2, ctd_subset_commas_09XY, ctd_subset_commas_09XY2)
paste("There are", length(unique(numerics_commas$AllSynonyms)), "terms with numerics and commas kept")

ctd_subset_reverse <- ctd_subset_commas[!(ctd_subset_commas$AllSynonyms %in% ctd_subset_commas_09$AllSynonyms),]
rm(ctd_subset_commas_09, ctd_subset_commas_09_first, ctd_subset_commas_09_second, ctd_subset_commas_09_nosep1, ctd_subset_commas_09X, ctd_subset_commas_09X2, ctd_subset_commas_09XY, ctd_subset_commas_09XY2)
ctd_subset_reverse$AllSynonyms <- lapply(strsplit(ctd_subset_reverse$AllSynonyms, split = ","), function(x) {paste(rev(trimws(x)), collapse = " ")})
paste("There are", length(unique(ctd_subset_reverse$AllSynonyms)), "terms with letters and commas reversed")

ctd_subset_commas <- rbind(numerics_commas, ctd_subset_reverse)
ctd_subset_commas$count <- NULL
ctd_database <- rbind(ctd_subset_commas, ctd_subset_nocommas)
ctd_database[] <- lapply(ctd_database, gsub, pattern = ",", replacement = " ", fixed = TRUE)
paste("There are", length(unique(ctd_database$AllSynonyms)), "synonyms")

# Remove any synonyms that are uppercase and less than 10 letters and convert the remaining uppercase synonyms to lower case
ctd_database$length <- nchar(ctd_database$AllSynonyms)
ctd_database$Upper <- toupper(ctd_database$AllSynonyms)
ctd_database <- ctd_database[!duplicated(ctd_database$AllSynonyms),]
ctd_Upper <- subset(ctd_database, ctd_database$AllSynonyms==ctd_database$Upper)
ctd_Upper <- subset(ctd_Upper, ctd_Upper$length<10)
ctd_database <- ctd_database[!(ctd_database$AllSynonyms %in% ctd_Upper$AllSynonyms),]
ctd_database$AllSynonyms <- tolower(ctd_database$AllSynonyms)
ctd_database$length <- NULL
ctd_database$Upper <- NULL
ctd_database <- unique(ctd_database)
ctd_database[] <- lapply(ctd_database, str_squish)
ctd_database <- ctd_database[!duplicated(ctd_database$AllSynonyms),]
paste("There are", length(unique(ctd_database$AllSynonyms)), "synonyms once those of less than length 10 were removed and duplicated, dashes, double spaces and synonyms converted to lowercase")

# Creation of the search terms for Entrez esearch and exact string match in Python
ctd_database$DiseaseName <- trimws(ctd_database$DiseaseName)
ctd_database$DiseaseID <- trimws(ctd_database$DiseaseID)
ctd_names <- ctd_database[,c(1,2)]
ctd_database$AllSynonyms <- trimws(ctd_database$AllSynonyms)
ctd_database <- ctd_database[,c(2,3)]
ctd_database <- as.data.frame(ctd_database)
ctd_database$DiseaseID <- as.character(ctd_database$DiseaseID)
ctd_database$AllSynonyms <- as.character(ctd_database$AllSynonyms)
length(unique(ctd_database$DiseaseID))
ctd_database <- unique(ctd_database)
ctd_database <- ctd_database[!duplicated(ctd_database$AllSynonyms), ]
ctd_database$New <- paste("r'\\b", ctd_database$AllSynonyms, sep="","\\b'")
ctd_database <- unique(ctd_database)
ctd_database <- ddply(ctd_database, .(DiseaseID), summarize, Xc=paste(unique(New), collapse =", "))
ctd_database$length <- nchar(ctd_database$Xc)
ctd_database <- unique(ctd_database)
colnames(ctd_database) <- c("All_MESH", "filter_frames", "length")
merge_frames <- merge(ctd_database, MeSH_Dataframe_final, by = "All_MESH", all.x=FALSE, all.y=TRUE)
merge_frames$length <- NULL
merge_frames <- merge_frames[order(merge_frames$DiseaseName_ARD_Paper),]
write.csv(merge_frames, "data/ARD_dictionary/merge_frames.csv")

# Manually edit the search terms from merge_frames.csv giving merge_frames_ARD.csv and then append them using code below
search_terms_final <- read.csv("data/ARD_dictionary/merge_frames_ARD.csv")
length(unique(search_terms_final$DiseaseName_ARD_Paper))
length(unique(search_terms_final$MESH_ID))
search_terms_final <- search_terms_final[!grepl("Infection, Other organisms", search_terms_final$DiseaseName_ARD_Paper),]
search_terms_final <- search_terms_final[!grepl("Infection, Other organs", search_terms_final$DiseaseName_ARD_Paper),]
search_terms_final <- search_terms_final[!grepl("Secondary Malignancy, other", search_terms_final$DiseaseName_ARD_Paper),]
search_terms_final <- search_terms_final[!grepl("Primary Malignancy, other", search_terms_final$DiseaseName_ARD_Paper),]
search_terms_final <- unique(search_terms_final)
search_terms_final[search_terms_final==""] <- NA
search_terms_final$Likelihood.disease.is.ageing.related <- trimws(search_terms_final$Likelihood.disease.is.ageing.related)
search_terms_final1 <- subset(search_terms_final, search_terms_final$Likelihood.disease.is.ageing.related=="High")
length(unique(search_terms_final1$DiseaseName_ARD_Paper))
search_terms_final2 <- data.frame(search_terms_final1$MESH_ID, search_terms_final1$DiseaseName_ARD_Paper, search_terms_final1$Search_term)
search_terms_final2$search_terms_final1.Search_term <- trimws(search_terms_final2$search_terms_final1.Search_term)
search_terms_final2$search_terms_final1.MESH_ID <- trimws(search_terms_final2$search_terms_final1.MESH_ID)
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- trimws(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper)
search_terms_final2 <- drop_na(search_terms_final2)
length(unique(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper))
length(unique(search_terms_final2$search_terms_final1.MESH_ID))
search_terms_final2$search_terms_final1.MESH_ID <- NULL
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- trimws(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper)
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_") 
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, " ", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, ",", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, "-", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, "__", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, "/", "_")
search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper <- str_replace(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, "'", "")
search_terms_final2 <- ddply(search_terms_final2, .(search_terms_final1.DiseaseName_ARD_Paper), summarize, Xc=paste(unique(search_terms_final1.Search_term), collapse =", "))
search_terms_final2$New <- paste(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper,"2", ' = ', "search_synonyms([", search_terms_final2$Xc, '], "', search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, '"', sep="",")")
search_terms_final2$per_disease <- paste(search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, ' = ', "per_disease2([", search_terms_final2$Xc, '], "', search_terms_final2$search_terms_final1.DiseaseName_ARD_Paper, '"', sep="",")")
write.csv(search_terms_final2, "data/ARD_dictionary/search_terms_final_new.csv")
new_file <- read.csv("data/ARD_dictionary/merge_frames_ARD.csv")
length(unique(new_file$DiseaseName_ARD_Paper))
new_file <- new_file[!grepl("Infection, Other organisms", new_file$DiseaseName_ARD_Paper),]
new_file <- new_file[!grepl("Infection, Other organs", new_file$DiseaseName_ARD_Paper),]
new_file <- new_file[!grepl("Secondary Malignancy, other", new_file$DiseaseName_ARD_Paper),]
new_file <- new_file[!grepl("Primary Malignancy, other", new_file$DiseaseName_ARD_Paper),]
new_file <- new_file[grepl("MESH:D", new_file$MESH_ID),]
new_file <- new_file[grepl("MESH:D", new_file$All_MESH),]
paste("The number of ARDs mapped to MeSH terms were", length(unique(new_file$DiseaseName_ARD_Paper)))
paste("The number of mapped MeSH terms was", length(unique(new_file$All_MESH_Terms)))
