### Code_Thesis_HCF/Aging_Hallmark_ARD_repository

Analysis in Chapters 3 and 4
- The custom code used for the analysis described in Chapters 3 and 4 of the thesis are in this folder: “Aging_hallmark_ARD_repository”. 
- The input datasets are available from the relevant subfolders. 
- The labelled sentences from manual curation are in the subfolder, “benchmark_sentences”. 
- The multimorbidity networks for the "Aging_hallmark_ARD_repository/multimorbidity networks" folder are available on request from Dr. Val Kuan (v.kuan@ucl.ac.uk). 
- Publicly available datasets used in the analysis are available for download from the respective website. 
- The following should be downloaded and added to the corresponding folder:
(i) 2019 PubMed baseline (https://www.nlm.nih.gov/databases/download/pubmed_medline.html) --> "Aging_hallmark_ARD_repository/pubmed_data/baseline"
(ii) CTD database's 'merged disease vocabulary' (http://ctdbase.org) --> "Aging_hallmark_ARD_repository/data/ARD_dictionary" and saved as "CTD_diseases_New.csv"
(iii) NCBI Gene Database of Homo sapiens (https://ftp.ncbi.nih.gov/gene/DATA/) --> "Aging_hallmark_ARD_repository/Genetic_data" and saved as "Homo_sapiens_gene_info.txt"
(iv) The STRING database (https://string-db.org/mapping_files/geneontology/, the human.GO_2_string.2018.tsv.gz file) --> "data/gene_ontology" and saved as "human.GO_2_string.2018.tsv"
(v) NHGRI-EBI GWAS catalog (https://www.ebi.ac.uk/gwas/) --> "Aging_hallmark_ARD_repository/Genetic_data" and saved as "GWAS_Catalog.tsv.zip"
(vi) NHGRI-EBI GWAS Ancestry data (https://www.ebi.ac.uk/gwas/) --> "Aging_hallmark_ARD_repository/Genetic_data" and saved as "Ancestry.tsv". 
(vii) The BioNetSmooth package was kindly provided by Prof. Andreas Beyer and is available from http://cellnet-sb.cecad.uni-koeln.de/resources/network-propagation/ and compatible with R version 3.3.0.


### Code_Thesis_HCF/Aging_Hallmark_ARD_repository
The code should be run in the following order: 

(i) Step1a_development_of_ARD_dictionary.R, R version 3.6.3
- Development of the ARD dictionary. 

(ii) Step1bi_GWAS_processing.R, R version 3.6.3
- Initial GWAS Ancestry data processing. 

(iii) Step1bii_converting_GWAS_ancestry.py, Python 3.7.0
- Continued GWAS Ancestry data processing. 

(iv) Step1c_ARD_GWAS_mappings.R, R version 3.6.3
- Mapping ARDs to genes using the GWAS catalog. 

(v) Step_1d_the_human_aging_corpus_PMID.py, Python 3.7.0
- Identifying the PubMed IDs representing the human aging corpus.

(vi) Step1e_processing_PubMed_baseline.R, R version 3.6.3
- Converting the PubMed baseline from XML format to CSV format. 

(vii) Step_2a_the_human_aging_corpus_data_and_metadata.py, Python 3.7.0
- Extracting titles, abstracts, and metadata for the human aging corpus. 

(viii) Step_2b_EDA_pubmed_baseline.py, Python 3.7.0
- Exploratory data analysis of the PubMed baseline. 

(ix) Step3_ochiai_and_manual_curation.py, Python 3.7.0
- Identifying sentences mentioning ARDs and aging hallmarks and co-mentioning ARDs and aging hallmarks. 
- Identifying the co-occurrence scores between ARDs and aging hallmarks. 
- Identifying the co-mentioning sentences and verifying aging hallmark-ARD associations using manual curation.

(x) Step4_heatmaps_of_Ochiai.R, R version 3.6.3
- Calculating the updated Ochiai coefficient.
- Identifying the top 30 ARDs associated with each aging hallmark.

(xi) Step5_create_heatmaps_genetic_data_updates.R, R version 3.6.3
- Identifying whether there was significant enrichment of gene ontology terms related to the same aging hallmark as was associated with the ARDs via literature mining. 
- Identifying whether any signalling pathways were significantly enriched across all aging hallmarks. Plotting circos plots to associate underlying genes to ARDs.

(xii) Step6_Multimorbidity_and_Network_Prop.R, R version 3.3.0
- Identifying whether aging hallmarks predict ARD multimorbidities. 
- Using network propagation to identify ARDs newly ranked in the top 30 ARDs associated with an aging hallmark. 

