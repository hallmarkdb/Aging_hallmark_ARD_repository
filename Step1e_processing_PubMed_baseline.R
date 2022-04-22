################################### Converting the PubMed Baseline XML to CSV files ###################################
# We created a folder named "baseline" for download of the 2019 PubMed baseline in XML format from "https://www.nlm.nih.gov/databases/download/pubmed_medline.html"
# We created a folder called "baseline_csv" file to store the csv file output
library(XML)
setwd("/Aging_hallmark_ARD_repository-main/")

extract_xml <- function(theFile){
  # This converts the PubMed baseline XML files into .csv files
  # This is a modified version of the extract_xml function in the pubmedXML package by Christopher Belter (i.e., https://github.com/christopherBelter/pubmedXML.R)
  newData <- xmlParse(theFile)
  records <- getNodeSet(newData, "//PubmedArticle")
  pmid <- xpathSApply(newData,"//MedlineCitation/PMID", xmlValue)
  year <- lapply(records, xpathSApply, ".//PubDate/Year", xmlValue) 
  year[sapply(year, is.list)] <- NA
  year[which(sapply(year, is.na) == TRUE)] <- lapply(records[which(sapply(year, is.na) == TRUE)], xpathSApply, ".//PubDate/MedlineDate", xmlValue)
  year <- gsub(" .+", "", year)
  year <- gsub("-.+", "", year)
  articletitle <- lapply(records, xpathSApply, ".//ArticleTitle", xmlValue) 
  articletitle[sapply(articletitle, is.list)] <- NA
  articletitle <- unlist(articletitle)
  abstract <- lapply(records, xpathSApply, ".//Abstract/AbstractText", xmlValue)
  abstract[sapply(abstract, is.list)] <- NA
  abstract <- sapply(abstract, paste, collapse = "|")
  meshHeadings <- lapply(records, xpathSApply, ".//DescriptorName", xmlValue)
  meshHeadings[sapply(meshHeadings, is.list)] <- NA
  meshHeadings <- sapply(meshHeadings, paste, collapse = "|")
  ptype <- lapply(records, xpathSApply, ".//PublicationType", xmlValue)
  ptype[sapply(ptype, is.list)] <- NA
  ptype <- sapply(ptype, paste, collapse = "|")
  theDF <- data.frame(pmid, year, articletitle, abstract, meshHeadings, ptype, stringsAsFactors = FALSE)
  return(theDF)
}

pubmed_csv <- function(number, route, location){
  # Saves csv files derived from the PubMed XML file
  file_string <- paste(route, number, ".xml", sep="")
  PubmedDoc <- extract_xml(file_string)
  file_string2 <- paste(location, number, ".csv", sep="")
  write.csv(PubmedDoc, file_string2)
  print(number)
  rm(PubmedDoc)
}

for (number in c(1:9)){
  print(pubmed_csv(number, route = "pubmed_data/baseline/pubmed19n000", location= "pubmed_data/baseline_csv/pubmed19n000"))
}

for (number in c(10:99)){
  print(pubmed_csv(number, route = "pubmed_data/baseline/pubmed19n00", location= "pubmed_data/baseline_csv/pubmed19n00"))
}

for (number in c(100:972)){
  print(pubmed_csv(number, route = "pubmed_data/baseline/pubmed19n0", location= "pubmed_data/baseline_csv/pubmed19n0"))
}
