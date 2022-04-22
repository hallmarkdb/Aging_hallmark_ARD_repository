library(ComplexHeatmap)
library(GOSim)
library(tidyr)
library(topGO)
library("hgu95av2.db")
library(dplyr)
library(plyr)
library(stringr)
library("ggpubr")
library(RColorBrewer)
library(GO.db)
library(circlize)
setwd("/Aging_hallmark_ARD_repository-main/")
# With thanks to Jan Grossbach from the Beyer lab for his advice on these analyses

#################################### Heatmaps of Genetic data and Circos Plots ############################################# 
top_30 <- function(a_frame, rep_number) {
  # This returns the top 30 ARDs for a given aging hallmark 
  # a_frame_col = "netmap_0.5_adjacency.disease_names", or "name_node"
  a <- a_frame[order(-a_frame[, rep_number]),]
  ochiai_frame_ARDs_AH <- data.frame(a[1:30,c("name_node")])
  colnames(ochiai_frame_ARDs_AH) <-c("name_node")
  return(ochiai_frame_ARDs_AH)
}

merge_frames <- function(top_frame) {
  # merge dataframes by "name_node" and assign genes to the top 30 ranked ARDs
  colnames(top_frame) <- c("name_node")
  edge_list_merge <- merge(ARD_Gene_frame_all, top_frame, by = "name_node", all.x =FALSE, all.y =FALSE)
  return(edge_list_merge)
}

plot_heatmap_pathway <- function(name, cluster_cols_true, number, hallmark, col_list){
  # Plot a heatmap for each aging hallmark
  name[is.na(name)] <- NA
  name <- as.data.frame(name)
  rownames(name) <- name$Term
  name$Term <- NULL
  AH <- name$AH
  name$AH <- NULL
  ha <- rowAnnotation("Pathway GO Term" = AH, col = col_list)
  name <- as.matrix((name))
  shape = dim(name)[1]/2 - number
  Heatmap(name, na_col = "white", cluster_rows=FALSE,cluster_columns =cluster_cols_true,
          col=rev(col2), left_annotation = ha,
          row_names_side = c("left"), row_names_gp = gpar(fontsize = 9),
          heatmap_width = unit(15, "cm"),rect_gp = gpar(col = "black", lwd = 0.7),
          width = NULL,
          heatmap_height = unit(shape, "cm"),
          height = NULL, heatmap_legend_param = list(
            title = "p-value",
            legend_height = unit(3, "cm"),
            title_position = "lefttop-rot"
          ))
}

plot_heatmap <- function(name, cluster_cols_true, number, hallmark, col_list){
  # Plot a heatmap for "pathway" and "cascade"
  name[is.na(name)] <- NA
  name <- as.data.frame(name)
  rownames(name) <- name$Term
  name$Term <- NULL
  AH <- name$AH
  name$AH <- NULL
  ha <- rowAnnotation("AH GO Term" = AH, col = col_list)
  name <- as.matrix((name))
  shape = dim(name)[1]/2 - number
  Heatmap(name, na_col = "white", cluster_rows=FALSE,cluster_columns =cluster_cols_true,
          col=rev(col2), left_annotation = ha,
          row_names_side = c("left"), row_names_gp = gpar(fontsize = 9),
          heatmap_width = unit(15, "cm"),rect_gp = gpar(col = "black", lwd = 0.7),
          width = NULL,
          heatmap_height = unit(shape, "cm"),
          height = NULL, heatmap_legend_param = list(
            title = "p-value",
            legend_height = unit(3, "cm"),
            title_position = "lefttop-rot"
          ))
}

draw_circos <- function(hallmark_df){
  # drawing circos plots
  hallmark_df <- data.frame(hallmark_df$Name, hallmark_df$Symbol)
  colnames(hallmark_df) <- c("to", "from")
  circos.clear()
  circos.par(gap.degree=1)
  chordDiagram(hallmark_df, grid.col= grid.col, annotationTrack = "grid", preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + 1.5, sector.name, facing = "clockwise", niceFacing = FALSE, adj = c(0, 0.5), cex = 0.9)
    circos.axis(h = "top", labels.cex = 0.1, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
  }, track.height = 0.15, bg.border = NA)
  return(hallmark_df)
}

star_rating <- function(data2) {
  # star ratings
  if (data2$weightFisher < 0.0001) {
    a <- "****"
  } else if (data2$weightFisher < 0.001) {
    a <- "***"
  } else if (data2$weightFisher < 0.01) {
    a <- "**"
  } else if (data2$weightFisher < 0.05){
    a <- "*"
  } else {
    a <- "null"
  }
  return(a)
}

list_to_matrix = function(lt, universal_set = NULL) {
  # list to matrix function developed previously and derived from the ComplexHeatmap package
  if(!is.null(universal_set)) {
    lt = lapply(lt, function(x) intersect(x, universal_set))
  } else {
    universal_set = unique(unlist(lt))
  }
  mat = matrix(0, nrow = length(universal_set), ncol = length(lt))
  rownames(mat) = sort(universal_set)
  colnames(mat) = names(lt)
  for(i in seq_along(lt)) {
    mat[as.character(unique(lt[[i]])), i] = 1
  }
  return(mat)
}

AH_search <- function(name) {
  new_frame <- enrichment_update[grepl(name, enrichment_update$Term, ignore.case = TRUE),]
  return (new_frame)
}

GI_list <- list("genome instability", "genomes instability", "genome instabilities", "genomic instability", 
                "genetic instability", "genetic instabilities", "damage dna", "dna damage", "damages dna", 
                "dna damages", "damage to dna", "dna injuries", "dna injury", "injury to dna", "injuries to dna", 
                "injures dna", "damage deoxyribonucleic acid", "damage to deoxyribonucleic acid", 
                "deoxyribonucleic acid damage", "damages deoxyribonucleic acid", "deoxyribonucleic acid damages", 
                "deoxyribonucleic acid injuries", "deoxyribonucleic acid injury", "injury to deoxyribonucleic acid", 
                "injuries to deoxyribonucleic acid", "injures deoxyribonucleic acid", "somatic mutation", "somatic mutations", 
                "genetic alteration", "genetic alterations", "dna break", "deoxyribonucleic acid break", "dna breaks", 
                "deoxyribonucleic acid breaks", "strand breaks", "strand break", "chromosomal breakage", "chromosome breakage", 
                "chromosomal breakages", "chromosome breakages", "chromosomal breaks", "chromosome breaks", "chromosomal break", 
                "chromosome break", "dna double strand break", "deoxyribonucleic acid double strand break", "double strand break", 
                "double stranded break", "double strand breaks", "double stranded breaks", "dna nick", "dna nicks", 
                "dna single strand break", "deoxyribonucleic acid nick", "deoxyribonucleic acid nicks", 
                "deoxyribonucleic acid single strand break", "single strand break", "single stranded break", 
                "single strand breaks", "single stranded breaks", "deficient dna repair", "deficient dna repairs", 
                "dna repair deficiency", "dna repair deficiencies", "dna repair disorder", "deficiency of dna repair", 
                "deficient deoxyribonucleic acid repair", "deficient deoxyribonucleic acid repairs", 
                "deoxyribonucleic acid repair deficiency", "deoxyribonucleic acid repair deficiencies", 
                "deoxyribonucleic acid repair disorder", "deficiency of deoxyribonucleic acid repair", "chromosome m damage", 
                "mitochondrial dna damage", "m chromosome damage", "mitochondrial chromosome damage", "dna mitochondrial damage", 
                "mdna damage", "mtdna damage", "damages chromosome m", "damages mitochondrial dna", "damages m chromosome", 
                "damages mitochondrial chromosome", "damages mitochondria dna", "damages dna mitochondrial", "damages mdna", 
                "damages mtdna", "damage chromosome m", "damage mitochondrial dna", "damage m chromosome", "damage mitochondrial chromosome", 
                "damage mitochondria dna", "damage dna mitochondrial", "damage mdna", "damage mtdna", "chromosome m injuries", 
                "m chromosome injuries", "mitochondrial chromosome injuries", "mdna injuries", "mtdna injuries", "chromosome m injury", 
                "m chromosome injury", "mitochondrial chromosome injury", "mitochondria dna injury", "mitochondrial dna injury", 
                "dna mitochondrial injury", "mdna injury", "mtdna injury", "injury to chromosome m", "injury to mitochondrial dna", 
                "injury to m chromosome", "injury to mitochondrial chromosome", "injury to mitochondria dna", "injury to dna mitochondrial", 
                "injury to mdna", "injury to mtdna", "injuries to chromosome m", "injuries to mitochondrial dna", "injuries to m chromosome", 
                "injuries to mitochondrial chromosome", "injuries to mitochondria dna", "injuries to dna mitochondrial", 
                "injuries to mdna", "injuries to mtdna", "injures chromosome m", "injures mitochondrial dna", "injures m chromosome", 
                "injures mitochondrial chromosome", "injures mitochondria dna", "injures dna mitochondrial", "injures mdna", 
                "injures mtdna", "transposable element", "transposable elements")

TA_list <- list("telomere attrition", "telomere erosion", "telomere exhaustion", "telomere loss", "telomere dysfunction", 
                "exhaustion of telomeres", "attrition of telomeres", "erosion of telomeres", "loss of telomeres", 
                "dysfunctional telomeres", "dysfunction of telomeres", "telomere shortening", "telomere shortenings", 
                "shortened telomeres", "short telomeres", "eroded telomeres", "shorter telomeres", "critically shortened telomeres", 
                "reduced ltl", "decreased ltl", "shortened ltl", "shorter ltl", "reduced mean ltl", "decreased mean ltl", "shortened mean ltl", 
                "shorter mean ltl", "telomere length", "length of telomeres", "reduced tl", "decreased tl", "shortened tl", "shorter tl", 
                "reduced mean tl", "decreased mean tl", "shortened mean tl", "shorter mean tl", "telomere", "telomeric")

EA_list <- list("dna methylation", "dna methylations", "methylation of dna", "deoxyribonucleic acid methylation", 
                "deoxyribonucleic acid methylations", "methylation of deoxyribonucleic acid", "methylated regions of dna", 
                "methylated region of deoxyribonucleic acid", "histone modification", "histone modifications", 
                "modification of histones", "modifications of histones", "histone acetylation", "acetylation of histones", 
                "histone acetylations", "acetylations of histones", "histone methylation", "methylation of histones", 
                "histone methylations", "methylations of histones", "epimutations", "epimutation", "epigenetics", "epigenetic", 
                "epigenetic process", "epigenetic processes", "epigenetic mechanism", "epigenetic alteration", "epigenetic alterations", 
                "epigenetic defect", "epigenetic defects", "epigenetic perturbation", "epigenetic perturbations", "epigenetic changes", 
                "epigenetic change", "changes in epigenetics", "alterations in epigenetics", "defects in epigenetics", "gene transcription", 
                "genetic transcription", "transcriptional alterations", "rna transcription", "transcriptional changes", 
                "transcriptional modifications", "ribonucleic acid transcription", "transcription of rna", "transcription of genes", 
                "transcription of ribonucleic acids", "transcription of ribonucleinicum", "non coding rna", "non peptide coding rna", 
                "noncoding rna", "nontranslated rna", "npcrna", "non protein coding rna", "untranslated rna", "functional rna", "ncrna", 
                "non coding rnas", "non peptide coding rnas", "noncoding rnas", "nontranslated rnas", "non protein coding rnas", 
                "untranslated rnas", "functional rnas", "non coding ribonucleic acid", "non peptide coding ribonucleic acid", 
                "noncoding ribonucleic acid", "nontranslated ribonucleic acid", "non protein coding ribonucleic acid", 
                "untranslated ribonucleic acid", "functional ribonucleic acid", "non coding ribonucleic acids", 
                "non peptide coding ribonucleic acids", "noncoding ribonucleic acids", "nontranslated ribonucleic acids", 
                "non protein coding ribonucleic acids", "untranslated ribonucleic acids", "functional ribonucleic acids", 
                "non coding ribonucleinicum acidum", "non peptide coding ribonucleinicum acidum", "noncoding ribonucleinicum acidum", 
                "nontranslated ribonucleinicum acidum", "non protein coding ribonucleinicum acidum", "untranslated ribonucleinicum acidum", 
                "functional ribonucleinicum acidum", "microrna", "mirna", "micro rna", "micrornas", "mirnas", "micro rnas", 
                "micro ribonucleic acid or micro ribonucleic acids", "micro ribonucleinicum acidum", "coding rna", "translated rna", 
                "coding rnas", "translated rnas", "coding ribonucleic acid", "translated ribonucleic acid", "coding ribonucleic acids", 
                "translated ribonucleic acids", "coding ribonucleinicum acidum", "translated ribonucleinicum acidum")

LOP_list <- list("protein homeostasis", "proteostasis", "proteostases", "endoplasmic reticulum stress", 
                 "endoplasmic reticulum stresses", "ergastoplasm stress", "ergastoplasm stresses", "er stress", 
                 "er stresses", "response to unfolded protein", "unfolded protein response", "unfolded protein responses", 
                 "proteolysis", "protein degradation", "protein degradations", "protein cleavage", "proteolytic", 
                 "proteasome activity", "autophagy", "macro autophagy", "macroautophagy", "autophagocytosis", 
                 "autophagies", "macroautophagies", "protein aggregates", "protein aggregation", "misfolding of proteins", 
                 "aggregation of proteins", "aggregates of proteins", "aggregated proteins", "misfolded proteins", "chaperone", 
                 "chaperones", "proteasome", "proteasomes")

DNS_list <- list("insulin resistance", "resistance to insulin", "intolerance to glucose", "glucose intolerance", 
                 "dyslipidemia", "dyslipidemias", "dyslipidaemia", "dyslipidaemias", "5 adenosine monophosphate activated kinase", 
                 "5 adenosine monophosphate activated kinases", "5 adenosine monophosphate activated protese", 
                 "5 adenosin kinaine monophosphate activated protein kinases", "5adenosine monophosphate activated kinase", 
                 "5adenosine monophosphate activated kinases", "5adenosine monophosphate activated protese", 
                 "5adenosin kinaine monophosphate activated protein kinases", "amp activated kinase", "amp activated kinases", 
                 "amp activated protein kinases", "amp activated protein kinase", "adenosine monophosphate activated kinase", 
                 "adenosine monophosphate activated kinases", "adenosine monophosphate activated protein kinase", 
                 "adenosine monophosphate activated protein kinases", "ampk", "prkaa", "sirt1", "sirtuin 1", 
                 "silent mating type information regulation 2 homolog 1", "mtorc1", "mtor complex 1", "torc1", 
                 "tor complex 1", "tor complex location 1", "torc 1", "rapamycin and nutrient sensitive tor complex", 
                 "target of rapamycin complex 1", "target of rapamycin complex", "insulin igf 1 signaling", 
                 "insulin igf 1 signalling", "insulin insulin like growth factor signalling", "insulin insulin like growth factor signaling", 
                 "insulin igf signaling", "insulin igf signalling", "insulin insulin like growth factor 1 signalling", 
                 "insulin insulin like growth factor 1 signaling", "insulin insulin like growth factor1 signalling", 
                 "insulin insulin like growth factor1 signaling", "insulin igf1 signaling", "insulin igf1 signalling", 
                 "iis signaling", "iis signalling", "nutrient detection", "nutrient sensing", "nutrient perception", 
                 "perception of nutrients", "detection of nutrients", "sensing of nutrients", "nutrient detection pathways", 
                 "nutrient sensing pathways", "nutrient perception pathways", "nutrient detection pathway", "nutrient sensing pathway", 
                 "nutrient perception pathway", " tor", " mtor", "insulin", "glucose", " lipid", "5 adenosine monophosphate activated kinase",
                 "5 adenosine monophosphate activated kinases", "5 adenosine monophosphate activated protese",
                 "5 adenosin kinaine monophosphate activated protein kinases", "5adenosine monophosphate activated kinase",
                 "5adenosine monophosphate activated kinases", "5adenosine monophosphate activated protese", 
                 "5adenosin kinaine monophosphate activated protein kinases", "amp activated kinase", "amp activated kinases", 
                 "amp activated protein kinases", "amp activated protein kinase", "adenosine monophosphate activated kinase", 
                 "adenosine monophosphate activated kinases", "adenosine monophosphate activated protein kinase", 
                 "adenosine monophosphate activated protein kinases", "ampk", "prkaa", "sirt1", "sirtuin 1", 
                 "silent mating type information regulation 2 homolog 1", "mtorc1", "mtor complex 1", "torc1",
                 "tor complex 1", "tor complex location 1", "torc 1", "rapamycin and nutrient sensitive tor complex",
                 "target of rapamycin complex 1", "target of rapamycin complex", "insulin igf 1 signaling", 
                 "insulin igf 1 signalling", "insulin insulin like growth factor signalling", "insulin insulin like growth factor signaling", 
                 "insulin igf signaling", "insulin igf signalling", "insulin insulin like growth factor 1 signalling", 
                 "insulin insulin like growth factor 1 signaling", "insulin insulin like growth factor1 signalling", 
                 "insulin insulin like growth factor1 signaling", "insulin igf1 signaling", "insulin igf1 signalling", 
                 "iis signaling", "iis signalling", "nutrient")

AIC_list <- list("inflammaging", "inflammageing", "inflamm aging", "inflamm ageing", "inflammation", "inflammations", 
                 "inflammatory", "intercellular communication", "cell cell signaling", "cell cell signalling", 
                 "cell cell communication", "cell to cell signaling", "cell to cell signalling", "cell to cell communication", 
                 "short range multicellular communication", "signalling between cells", "signalling between cells", 
                 "communication between cells", "hormone signal", "endocrine signal", "hormonal signal", "hormone signals", 
                 "endocrine signals", "hormonal signals", "hormone signaling", "endocrine signaling", "hormonal signaling", 
                 "hormone signalling", "endocrine signalling", "hormonal signalling", "inflammatory signal", "inflammatory conditions", 
                 "inflammatory signals", "inflammatory signalling", "inflammatory signaling", "neural signal", "neural signals", 
                 "neural signalling", "neural signaling", "neuronal signal", "neuronal signals", "neuronal signalling", 
                 "neuronal signaling", "nerve signal", "nerve signals", "nerve signalling", "nerve signaling", "neurotransmission", 
                 "synaptic transmission", "nerve transmission", "nerve impulse transmission", "neural transmission", "neurotransmissions", 
                 "conduction of nerve impulse", "signal transmission across a synapse", "signal transmission along a neuron", 
                 "transmission of nerve impulse", "neuronal transmission", "neural signal transduction", "nerve signal transduction", 
                 "nerve transmitter substances", "nerve transmitter substance", "neurotransmitte", "neurotransmitters", "neurosteroids", 
                 "neurosteroid", "neuroregulators", "neuroregulato", "neurohumors", "neurohumo", "neuromodulaters", "neuromodulato", 
                 "neurohormones", "neurohormone", "endocrine gland secretion", "endocrine gland secretions", "hormone", "hormones", "hormonal")

MD_list <- list("mitochondri", "impaired mitochondrial function", "mitochondrial dysfunction", "dysfunction of mitochondria", 
                "dysfunctional mitochondria", "mitochondrial deterioration", "deterioration of the mitochondria", 
                "mitochondrial degeneration", "degeneration of mitcohondria", "mitochondrial deficiencies", "mitochondrial decline", 
                "decline in mitochondrial function", "mitochondrial toxicity", "mitochondrial damage", "toxic to mitochondria", 
                "damages mitochondria", "damage to mitochondria", "oxygen radicals", "oxygen radical", "oxygen reactive species", 
                "pro oxidant", "reactive oxygen species", "pro oxidants", "active oxygen", "mitochondrial bioenergetics", 
                "mitochondrial bioenergetic", "mitochondria bioenergetics", "mitochondria bioenergetic", "mitochondrial turnover", 
                "turnover of mitochondria", "mitochondria turnover", "mitochondrial biogenesis", "biogenesis of mitochondria", 
                "mitochondria biogenesis", "mitochondrial degradation", "degradation of mitochondria", "mitochondria degradation", 
                "mitophagy", "autophagy of mitochondrion", "macromitophagy", "mitochondrion autophagy", "mitochondria autophagy", 
                "mitochondrial dynamics", "mitochondrial dynamic", "mitochondrial fission", "mitochondrial division", 
                "mitochondrion in division", "mitochondrial fissions", "mitochondrial fusion", "mitochondrial fusions", 
                "mitochondria dynamics", "mitochondria fussion", "mitochondria division", "mitochondria in division", 
                "mitochondria fission", "mitochondria fusion", "mitochondria fusions", "electron transport chain", 
                "respiratory chain", "citric acid cycle", "citric acid cycles", "krebs cycle", "tricarboxylic acid cycle", 
                "tricarboxylic acid cycles", "tca cycle", "acid citric cycle", "cycle kreb", "cycle krebs", "cycles krebs", "kreb cycle")

CS_list <- list("cell ageing", "cellular ageing", "cell aging", "cellular aging", "cell senescence", 
                "aging cells", "cellular senescence", "replicative senescence", "age cells", "ageing cell", 
                "aging cell", "cell age", "senescent cell", "senescent cells", "immunosenescence", "immunosenescent", 
                "senescence associated secretory phenotype", "sasp", "senescence marker", "senescence markers", 
                "markers of cellular senescence", "marker of cellular senescence", "markers of senescence", 
                "marker of senescence", "markers of cell senescence", "marker of cell senescence")

grid.col = c("AF" = "gold", "Breast Cancer" = "green", "Bladder Cancer"="yellow",
             "Cirrhosis" = "cyan", "Dementia" = "red", "Hypertension" = "purple", 
             "Giant Cell Arteritis" = "blue", "CHD" = "powderblue", "COPD" = "skyblue", 
             "Glomerulonephritis" = "darkseagreen", "LRTI" = "black", "Low HDL" = "plum", 
             "Leukaemia" = "blue", "Myelodysplasia" = "indianred1", "CHF" = "hotpink",
             "NHL" = "navy", "Plasma Cell Cancer" = "turquoise1", "Lung Cancer" = "brown", 
             "Melanoma" = "darkorange", "Prostate Cancer" = "yellow", "Skin Cancer" = "blue", 
             "Ovarian Cancer" = "ivory2", "Pancreatic Cancer" = "deepskyblue",
             "Uterine Cancer" = "plum1", "Glaucoma" = "darkgreen", "AMD" = "thistle", 
             "Myasthenia Gravis" = "maroon1", "Osteoporosis" = "chocolate1", "PAD" = "cyan", 
             "Osteoarthritis" = "navyblue", "Parkinson's" = "brown", "Stroke" = "blue4", 
             "Raised LDL" = "orange", "High Cholesterol" = "black", "RA" = "purple", 
             "Scleroderma" = "black", "Sjogren's" = "cyan", "T2DM" = "deeppink", 
             "Thyroid Disease" = "darkgreen", "T2DM" = "hotpink", "Ischaemic Stroke"="darkblue")

# Import the relevant data frames
ochiai_frame_ARDs <- read.csv('data/aging_hallmark_subnetworks/ochiai_frame_0321.csv', header=TRUE)
ochiai_frame_ARDs$X.1 <- NULL
mappings_ARD_abbreviation <- read.csv("data/aging_hallmark_subnetworks/mappings_ARD_abbreviation.csv")
colnames(mappings_ARD_abbreviation) <- c("X", "name_node", "Name")
ochiai_frame_ARDs <- merge(ochiai_frame_ARDs, mappings_ARD_abbreviation, by = "X")
ochiai_frame_ARDs$name_node <- trimws(ochiai_frame_ARDs$name_node)
ochiai_frame_ARDs <- ochiai_frame_ARDs[,c(1, 48, 49, 2:47)]
ochiai_frame_ARDs_copy <- ochiai_frame_ARDs

ARD_Gene_frame_all <- read.csv("Genetic_data/ARD_Gene_frame_all_NCBI_50000.csv")
ARD_Gene_frame_all$X <- NULL
ARD_Gene_frame_all$name_node <- trimws(ARD_Gene_frame_all$name_node)

mappings_file <- read.csv("Genetic_data/human.entrez_2_string.2018.tsv", sep="\t")
mappings_file <- subset(mappings_file, mappings_file$NCBI.taxid=="9606")
mappings_file$NCBI.taxid <- NULL
mappings_file[] <- lapply(mappings_file, gsub, pattern = "|", replacement = ";", fixed = TRUE)
mappings_file <- separate_rows(mappings_file, sep=";", entrez)
colnames(mappings_file) <- c("GeneID", "stringId")
mappings_file$GeneID <- as.integer(mappings_file$GeneID)

# Creating the gene ontology (GO) background set
gene_ontology_string <- read.csv("data/gene_ontology/human.GO_2_string.2018.tsv", sep="\t", header=FALSE)
colnames(gene_ontology_string) <- c("NCBI_taxid", "category", "geneontology", "STRING")
go <- subset(gene_ontology_string, gene_ontology_string$NCBI_taxid =="9606")
go$NCBI_taxid <- NULL
paste("There were ", length(unique(go$STRING)), "stringIds in the GO download")
go <- subset(go, go$category=="Process")
go$category <- NULL
go <- unique(go)
go <- unique(na.omit(go))
go <- ddply(go, .(STRING), summarize, Xc=paste(unique(geneontology), collapse=", "))
write.table(go, "data/gene_ontology/go_ontology.map", sep="\t", row.names = F, col.names=F) 
#Remove quote marks from the "go_ontology.map" file in a text editor
geneID2GO <- readMappings("data/gene_ontology/go_ontology.map",sep="\t", IDsep=",")
geneUniverse <- names(geneID2GO)
paste("There are", length(unique(geneUniverse)), "genes in the gene universe mapped to biological processes")
geneID2GO_2 <- list_to_matrix(geneID2GO)
geneUniverse2 <- data.frame(geneUniverse)
colnames(geneUniverse2) <- "stringId"

# Protein list for the top 30 ARDs for genomic instability
top_GI <- top_30(ochiai_frame_ARDs, "GI_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_GI$name_node)))
GI <- unique(merge_frames(top_GI))
GI_copy <- GI
GI <- data.frame(unique(GI_copy$GeneID))
GI_copy$GeneID <- as.integer(GI_copy$GeneID)
GI_string <- merge(mappings_file, GI_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the genomic instability gene set is", length(unique(GI_string$name_node)))
gi_count <- merge(GI_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(gi_count$stringId))
length(unique(gi_count$name_node))
length(unique(gi_count$Symbol))

# Protein list for the top 30 ARDs for telomere attrition 
top_TA <- top_30(ochiai_frame_ARDs, "TA_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_TA$name_node)))
TA <- unique(merge_frames(top_TA))
TA_copy <- TA
TA <- data.frame(unique(TA$GeneID))
TA_copy$GeneID <- as.integer(TA_copy$GeneID)
TA_string <- merge(mappings_file, TA_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the telomere attrition gene set is", length(unique(TA_string$name_node)))
ta_count <- merge(TA_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(ta_count$stringId))
length(unique(ta_count$name_node))
setdiff(unique(ta_count$name_node), unique(gi_count$name_node))

# Protein list for the top 30 ARDs for epigenetic alterations 
top_EA <- top_30(ochiai_frame_ARDs, "EA_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_EA$name_node)))
EA <- unique(merge_frames(top_EA))
EA_copy <- EA
EA <- data.frame(unique(EA$GeneID))
EA_copy$GeneID <- as.integer(EA_copy$GeneID)
EA_string <- merge(mappings_file, EA_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the epigenetic alterations gene set is", length(unique(EA_string$name_node)))
ea_count <- merge(EA_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(ea_count$stringId))
length(unique(ea_count$name_node))

# Protein list for the top 30 ARDs for loss of proteostasis 
top_LOP <- top_30(ochiai_frame_ARDs, "LOP_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_LOP$name_node)))
LOP <- unique(merge_frames(top_LOP))
LOP_copy <- LOP
LOP <- data.frame(unique(LOP$GeneID))
LOP_copy$GeneID <- as.integer(LOP_copy$GeneID)
LOP_string <- merge(mappings_file, LOP_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the loss of proteostasis gene set is", length(unique(LOP_string$name_node)))
lop_count <- merge(LOP_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(lop_count$stringId))
length(unique(lop_count$name_node))

# Protein list for the top 30 ARDs for deregulated nutrient sensing 
top_DNS <- top_30(ochiai_frame_ARDs, "DNS_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_DNS$name_node)))
DNS <- unique(merge_frames(top_DNS))
DNS_copy <- DNS
DNS <- data.frame(unique(DNS$GeneID))
DNS_copy$GeneID <- as.integer(DNS_copy$GeneID)
DNS_string <- merge(mappings_file, DNS_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the deregulated nutrient sensing gene set is", length(unique(DNS_string$name_node)))
dns_count <- merge(DNS_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(dns_count$stringId))
length(unique(dns_count$name_node))

# Protein list for the top 30 ARDs for mitochondrial dysfunction 
top_MD <- top_30(ochiai_frame_ARDs, "MD_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_MD$name_node)))
MD <- unique(merge_frames(top_MD))
MD_copy <- MD
MD <- data.frame(unique(MD$GeneID))
MD_copy$GeneID <- as.integer(MD_copy$GeneID)
MD_string <- merge(mappings_file, MD_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the mitochondrial dysfunction gene set is", length(unique(MD_string$name_node)))
md_count <- merge(MD_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(md_count$stringId))
length(unique(md_count$name_node))

# Protein list for the top 30 ARDs for stem cell exhaustion
top_SCE <- top_30(ochiai_frame_ARDs, "SCE_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_SCE$name_node)))
SCE <- unique(merge_frames(top_SCE))
SCE_copy <- SCE
SCE <- data.frame(unique(SCE$GeneID))
SCE_copy$GeneID <- as.integer(SCE_copy$GeneID)
SCE_string <- merge(mappings_file, SCE_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the stem cell exhaustion gene set is", length(unique(SCE_string$name_node)))
sce_count <- merge(SCE_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(sce_count$stringId))
length(unique(sce_count$name_node))

# Protein list for the top 30 ARDs for cellular senescence 
top_CS <- top_30(ochiai_frame_ARDs, "CS_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_CS$name_node)))
CS <- unique(merge_frames(top_CS))
CS_copy <- CS
CS <- data.frame(unique(CS$GeneID))
CS_copy$GeneID <- as.integer(CS_copy$GeneID)
CS_string <- drop_na(merge(mappings_file, CS_copy, by = "GeneID", all.x = FALSE, all.y = FALSE))
paste("The number of ARDs represented in the cellular senescence gene set is", length(unique(CS_string$name_node)))
cs_count <- merge(CS_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(cs_count$stringId))
length(unique(cs_count$name_node))
length(unique(cs_count$Symbol))
cs_count2 <- unique(data.frame(cs_count$stringId, cs_count$Symbol))
cs_count2[duplicated(cs_count2$cs_count.stringId),]

# Protein list for the top 30 ARDs for altered intercellular communication
top_AIC <- top_30(ochiai_frame_ARDs, "AIC_multiple")
paste("The number of ARDs included in the initial analysis was ", length(unique(top_AIC$name_node)))
AIC <- unique(merge_frames(top_AIC))
AIC_copy <- AIC
AIC <- data.frame(unique(AIC$GeneID))
AIC_copy$GeneID <- as.integer(AIC_copy$GeneID)
AIC_string <- merge(mappings_file, AIC_copy, by = "GeneID", all.x = FALSE, all.y = FALSE)
paste("The number of ARDs represented in the altered intercellular communication gene set is", length(unique(AIC_string$name_node)))
aic_count <- merge(AIC_string, geneUniverse2, by = "stringId", all.x=FALSE, all.y=FALSE)
length(unique(aic_count$stringId))
length(unique(aic_count$name_node))
length(unique(aic_count$Symbol))
aic_count2 <- unique(data.frame(aic_count$stringId, aic_count$Symbol))
aic_count2[duplicated(aic_count2$aic_count.stringId),]

full_set <- union(union(union(union(union(union(union(union(top_GI, top_TA), top_EA), top_LOP), top_DNS), top_CS), top_SCE), top_AIC), top_MD)
full_set2 <- union(union(union(union(union(union(union(union(GI_string, TA_string), EA_string), LOP_string), DNS_string), CS_string), SCE_string), AIC_string), MD_string)
full_set3 <- union(union(union(union(union(union(union(union(gi_count, ta_count), ea_count), lop_count), dns_count), cs_count), sce_count), aic_count), md_count)

genes_all <- rbind(GI_string, TA_string, EA_string, LOP_string, DNS_string, CS_string, SCE_string, MD_string, AIC_string)
paste("The number of top 30 ARDs represented were", length(unique(full_set2$name_node)))
paste("The number of Gene IDs represented were", length(unique(full_set2$GeneID)))
paste("The number of stringIds represented were", length(unique(full_set2$stringId)))
paste("The total number of top 30 ARDs were ", length(unique(full_set$name_node)))
paste("The number of top 30 ARDs represented were (only those also in background set)", length(unique(full_set3$name_node)))
paste("The number of Gene IDs represented were (only those also in background set)", length(unique(full_set3$GeneID)))
paste("The number of stringIds represented were (only those also in background set)", length(unique(full_set3$stringId)))

gi_genesOfInterest <- unique(as.character(GI_string$stringId))
ta_genesOfInterest <- unique(as.character(TA_string$stringId))
ea_genesOfInterest <- unique(as.character(EA_string$stringId))
lop_genesOfInterest <- unique(as.character(LOP_string$stringId))
cs_genesOfInterest <- unique(as.character(CS_string$stringId))
aic_genesOfInterest <- unique(as.character(AIC_string$stringId))
md_genesOfInterest <- unique(as.character(MD_string$stringId))
sce_genesOfInterest <- unique(as.character(SCE_string$stringId))
dns_genesOfInterest <- unique(as.character(DNS_string$stringId))

# Genomic instability based on text mining
gi_geneList <- as.factor(as.integer(geneUniverse %in% gi_genesOfInterest))
sum(as.integer(geneUniverse %in% gi_genesOfInterest))
names(gi_geneList) <- geneUniverse
myGOdata_gi <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = gi_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_GI.weight <- runTest(myGOdata_gi, algorithm="weight01", statistic="fisher")
allRes_GI <- GenTable(myGOdata_gi, weightFisher = resultFisher_GI.weight, orderBy = "weightFisher", topNodes = 5000)
allRes_GI <- subset(allRes_GI, allRes_GI$weightFisher < 0.05)
allRes_GI$weightFisher <- as.numeric(allRes_GI$weightFisher)
allRes_IFN_gamma_GI <- allRes_GI[grep("GO:0060333", allRes_GI$GO.ID),]
allRes_T_Cell_GI <- allRes_GI[grep("GO:0050852", allRes_GI$GO.ID),] 
allRes_Pos_T_Cell_GI  <- allRes_GI[grep("GO:0050862", allRes_GI$GO.ID),] 
allRes_GI_ERK <- allRes_GI[grep("GO:0070374", allRes_GI$GO.ID),]
allRes_Intrinsic_Apoptotic_GI <- allRes_GI[grep("GO:0042771", allRes_GI$GO.ID),] 
GI_table <- rbind(allRes_IFN_gamma_GI, allRes_T_Cell_GI, allRes_Pos_T_Cell_GI, allRes_GI_ERK, allRes_Intrinsic_Apoptotic_GI)
star_rating(allRes_T_Cell_GI)
star_rating(allRes_IFN_gamma_GI)
star_rating(allRes_Pos_T_Cell_GI)
star_rating(allRes_GI_ERK)
star_rating(allRes_Intrinsic_Apoptotic_GI)

# Telomere attrition based on text mining
ta_geneList <- as.factor(as.integer(geneUniverse %in% ta_genesOfInterest))
sum(as.integer(geneUniverse %in% ta_genesOfInterest))
names(ta_geneList) <- geneUniverse
myGOdata_ta <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = ta_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_TA.weight <- runTest(myGOdata_ta, algorithm="weight01", statistic="fisher")
allRes_TA <- GenTable(myGOdata_ta, weightFisher = resultFisher_TA.weight, orderBy="weightFisher", topNodes = 5000)
allRes_TA$weightFisher <- as.numeric(allRes_TA$weightFisher)
allRes_TA <- subset(allRes_TA, allRes_TA$weightFisher < 0.05)
allRes_IFN_gamma_TA <- allRes_TA[grep("GO:0060333", allRes_TA$GO.ID),]
allRes_T_Cell_TA <- allRes_TA[grep("GO:0050852", allRes_TA$GO.ID),] 
allRes_Pos_T_Cell_TA  <- allRes_TA[grep("GO:0050862", allRes_TA$GO.ID),] 
allRes_TA_ERK <- allRes_TA[grep("GO:0070374", allRes_TA$GO.ID),]
allRes_Intrinsic_Apoptotic_TA <- allRes_TA[grep("GO:0042771", allRes_TA$GO.ID),] 
TA_table <- rbind(allRes_IFN_gamma_TA, allRes_T_Cell_TA, allRes_Pos_T_Cell_TA, allRes_TA_ERK, allRes_Intrinsic_Apoptotic_TA)
star_rating(allRes_IFN_gamma_TA)
star_rating(allRes_T_Cell_TA)
star_rating(allRes_Pos_T_Cell_TA)
star_rating(allRes_TA_ERK)
star_rating(allRes_Intrinsic_Apoptotic_TA)

# Epigenetic alterations based on text mining
ea_geneList <- as.factor(as.integer(geneUniverse %in% ea_genesOfInterest))
sum(as.integer(geneUniverse %in% ea_genesOfInterest))
names(ea_geneList) <- geneUniverse
myGOdata_ea <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = ea_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_EA.weight <- runTest(myGOdata_ea, algorithm="weight01", statistic="fisher")
allRes_EA <- GenTable(myGOdata_ea, weightFisher = resultFisher_EA.weight, orderBy="weightFisher", topNodes = 5000)
allRes_EA$weightFisher <- as.numeric(allRes_EA$weightFisher)
allRes_EA <- subset(allRes_EA, allRes_EA$weightFisher < 0.05)
allRes_IFN_gamma_EA <- allRes_EA[grep("GO:0060333", allRes_EA$GO.ID),]
allRes_T_Cell_EA <- allRes_EA[grep("GO:0050852", allRes_EA$GO.ID),] 
allRes_Pos_T_Cell_EA  <- allRes_EA[grep("GO:0050862", allRes_EA$GO.ID),] 
allRes_EA_ERK <- allRes_EA[grep("GO:0070374", allRes_EA$GO.ID),]
allRes_Intrinsic_Apoptotic_EA <- allRes_EA[grep("GO:0042771", allRes_EA$GO.ID),] 
EA_table <- rbind(allRes_IFN_gamma_EA, allRes_T_Cell_EA, allRes_Pos_T_Cell_EA, allRes_EA_ERK, allRes_Intrinsic_Apoptotic_EA)
star_rating(allRes_IFN_gamma_EA)
star_rating(allRes_T_Cell_EA)
star_rating(allRes_Pos_T_Cell_EA)
star_rating(allRes_EA_ERK)
star_rating(allRes_Intrinsic_Apoptotic_EA)

# Loss of proteostasis based on text mining
lop_geneList <- as.factor(as.integer(geneUniverse %in% lop_genesOfInterest))
sum(as.integer(geneUniverse %in% lop_genesOfInterest))
names(lop_geneList) <- geneUniverse
myGOdata_lop <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = lop_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_LOP.weight <- runTest(myGOdata_lop, algorithm="weight01", statistic="fisher")
allRes_LOP <- GenTable(myGOdata_lop, weightFisher = resultFisher_LOP.weight, orderBy="weightFisher", topNodes = 5000)
allRes_LOP$weightFisher <- as.numeric(allRes_LOP$weightFisher)
allRes_LOP <- subset(allRes_LOP, allRes_LOP$weightFisher < 0.05)
allRes_IFN_gamma_LOP <- allRes_LOP[grep("GO:0060333", allRes_LOP$GO.ID),]
allRes_T_Cell_LOP <- allRes_LOP[grep("GO:0050852", allRes_LOP$GO.ID),] 
allRes_Pos_T_Cell_LOP  <- allRes_LOP[grep("GO:0050862", allRes_LOP$GO.ID),] 
allRes_LOP_ERK <- allRes_LOP[grep("GO:0070374", allRes_LOP$GO.ID),]
allRes_Intrinsic_Apoptotic_LOP <- allRes_LOP[grep("GO:0042771", allRes_LOP$GO.ID),] 
LOP_table <- rbind(allRes_IFN_gamma_LOP, allRes_T_Cell_LOP, allRes_Pos_T_Cell_LOP, allRes_LOP_ERK, allRes_Intrinsic_Apoptotic_LOP)
star_rating(allRes_IFN_gamma_LOP)
star_rating(allRes_T_Cell_LOP)
star_rating(allRes_Pos_T_Cell_LOP)
star_rating(allRes_LOP_ERK)
star_rating(allRes_Intrinsic_Apoptotic_LOP)

# Deregulated nutrient sensing based on text mining
dns_geneList <- as.factor(as.integer(geneUniverse %in% dns_genesOfInterest))
sum(as.integer(geneUniverse %in% dns_genesOfInterest))
names(dns_geneList) <- geneUniverse
myGOdata_dns <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = dns_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_DNS.weight <- runTest(myGOdata_dns, algorithm="weight01", statistic="fisher")
allRes_DNS <- GenTable(myGOdata_dns, weightFisher = resultFisher_DNS.weight, orderBy="weightFisher", topNodes=5000)
allRes_DNS$weightFisher <- as.numeric(allRes_DNS$weightFisher)
allRes_DNS <- subset(allRes_DNS, allRes_DNS$weightFisher < 0.05)
allRes_IFN_gamma_DNS <- allRes_DNS[grep("GO:0060333", allRes_DNS$GO.ID),]
allRes_T_Cell_DNS <- allRes_DNS[grep("GO:0050852", allRes_DNS$GO.ID),] 
allRes_Pos_T_Cell_DNS  <- allRes_DNS[grep("GO:0050862", allRes_DNS$GO.ID),] 
allRes_DNS_ERK <- allRes_DNS[grep("GO:0070374", allRes_DNS$GO.ID),]
allRes_Intrinsic_Apoptotic_DNS <- allRes_DNS[grep("GO:0042771", allRes_DNS$GO.ID),]
DNS_table <- rbind(allRes_IFN_gamma_DNS, allRes_T_Cell_DNS, allRes_Pos_T_Cell_DNS, allRes_DNS_ERK, allRes_Intrinsic_Apoptotic_DNS)
star_rating(allRes_IFN_gamma_DNS)
star_rating(allRes_T_Cell_DNS)
star_rating(allRes_Pos_T_Cell_DNS)
star_rating(allRes_DNS_ERK)
star_rating(allRes_Intrinsic_Apoptotic_DNS)

# Mitochondrial dysfunction based on text mining
md_geneList <- as.factor(as.integer(geneUniverse %in% md_genesOfInterest))
sum(as.integer(geneUniverse %in% md_genesOfInterest))
names(md_geneList) <- geneUniverse
myGOdata_md <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = md_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_MD.weight <- runTest(myGOdata_md, algorithm="weight01", statistic="fisher")
allRes_MD <- GenTable(myGOdata_md, weightFisher = resultFisher_MD.weight, orderBy="weightFisher", topNodes=5000)
allRes_MD$weightFisher <- as.numeric(allRes_MD$weightFisher)
allRes_MD <- subset(allRes_MD, allRes_MD$weightFisher < 0.05)
allRes_IFN_gamma_MD <- allRes_MD[grep("GO:0060333", allRes_MD$GO.ID),]
allRes_T_Cell_MD <- allRes_MD[grep("GO:0050852", allRes_MD$GO.ID),] 
allRes_Pos_T_Cell_MD  <- allRes_MD[grep("GO:0050862", allRes_MD$GO.ID),] 
allRes_MD_ERK <- allRes_MD[grep("GO:0070374", allRes_MD$GO.ID),]
allRes_Intrinsic_Apoptotic_MD <- allRes_MD[grep("GO:0042771", allRes_MD$GO.ID),]
MD_table <- rbind(allRes_IFN_gamma_MD, allRes_T_Cell_MD, allRes_Pos_T_Cell_MD, allRes_MD_ERK, allRes_Intrinsic_Apoptotic_MD)
star_rating(allRes_IFN_gamma_MD)
star_rating(allRes_T_Cell_MD)
star_rating(allRes_Pos_T_Cell_MD)
star_rating(allRes_MD_ERK)
star_rating(allRes_Intrinsic_Apoptotic_MD)

# Cellular senescence based on text mining
cs_geneList <- as.factor(as.integer(geneUniverse %in% cs_genesOfInterest))
sum(as.integer(geneUniverse %in% cs_genesOfInterest))
names(cs_geneList) <- geneUniverse
myGOdata_cs <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = cs_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_CS.weight <- runTest(myGOdata_cs, algorithm="weight01", statistic="fisher")
allRes_CS <- GenTable(myGOdata_cs, weightFisher = resultFisher_CS.weight, orderBy="weightFisher", topNodes=5000)
allRes_CS$weightFisher <- as.numeric(allRes_CS$weightFisher)
allRes_CS <- subset(allRes_CS, allRes_CS$weightFisher < 0.05)
allRes_IFN_gamma_CS <- allRes_CS[grep("GO:0060333", allRes_CS$GO.ID),]
allRes_T_Cell_CS <- allRes_CS[grep("GO:0050852", allRes_CS$GO.ID),] 
allRes_Pos_T_Cell_CS  <- allRes_CS[grep("GO:0050862", allRes_CS$GO.ID),] 
allRes_CS_ERK <- allRes_CS[grep("GO:0070374", allRes_CS$GO.ID),]
allRes_Intrinsic_Apoptotic_CS <- allRes_CS[grep("GO:0042771", allRes_CS$GO.ID),]
CS_table <- rbind(allRes_IFN_gamma_CS, allRes_T_Cell_CS, allRes_Pos_T_Cell_CS, allRes_CS_ERK, allRes_Intrinsic_Apoptotic_CS)
star_rating(allRes_IFN_gamma_CS)
star_rating(allRes_T_Cell_CS)
star_rating(allRes_Pos_T_Cell_CS)
star_rating(allRes_CS_ERK)
star_rating(allRes_Intrinsic_Apoptotic_CS)

# Stem cell exhaustion based on text mining
sce_geneList <- as.factor(as.integer(geneUniverse %in% sce_genesOfInterest))
sum(as.integer(geneUniverse %in% sce_genesOfInterest))
names(sce_geneList) <- geneUniverse
myGOdata_sce <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = sce_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_SCE.weight <- runTest(myGOdata_sce, algorithm="weight01", statistic="fisher")
allRes_SCE <- GenTable(myGOdata_sce, weightFisher = resultFisher_SCE.weight, orderBy="weightFisher", topNodes=5000)
allRes_SCE$weightFisher <- as.numeric(allRes_SCE$weightFisher)
allRes_SCE <- subset(allRes_SCE, allRes_SCE$weightFisher < 0.05)
allRes_IFN_gamma_SCE <- allRes_SCE[grep("GO:0060333", allRes_SCE$GO.ID),]
allRes_T_Cell_SCE <- allRes_SCE[grep("GO:0050852", allRes_SCE$GO.ID),] 
allRes_Pos_T_Cell_SCE  <- allRes_SCE[grep("GO:0050862", allRes_SCE$GO.ID),] 
allRes_SCE_ERK <- allRes_SCE[grep("GO:0070374", allRes_SCE$GO.ID),]
allRes_Intrinsic_Apoptotic_SCE <- allRes_SCE[grep("GO:0042771", allRes_SCE$GO.ID),]
SCE_table <- rbind(allRes_IFN_gamma_SCE, allRes_T_Cell_SCE, allRes_Pos_T_Cell_SCE, allRes_SCE_ERK, allRes_Intrinsic_Apoptotic_SCE)
star_rating(allRes_IFN_gamma_SCE)
star_rating(allRes_T_Cell_SCE)
star_rating(allRes_Pos_T_Cell_SCE)
star_rating(allRes_SCE_ERK)
star_rating(allRes_Intrinsic_Apoptotic_SCE)

# Altered intercellular communication based on text mining
aic_geneList <- as.factor(as.integer(geneUniverse %in% aic_genesOfInterest))
sum(as.integer(geneUniverse %in% aic_genesOfInterest))
names(aic_geneList) <- geneUniverse
myGOdata_aic <- new("topGOdata", description="My project", ontology=c("BP"), allGenes = aic_geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO, nodeSize=5)
resultFisher_AIC.weight <- runTest(myGOdata_aic, algorithm="weight01", statistic="fisher")
allRes_AIC <- GenTable(myGOdata_aic, weightFisher = resultFisher_AIC.weight, orderBy="weightFisher", topNodes=5000)
allRes_AIC$weightFisher <- as.numeric(allRes_AIC$weightFisher)
allRes_AIC <- subset(allRes_AIC, allRes_AIC$weightFisher < 0.05)
allRes_IFN_gamma_AIC <- allRes_AIC[grep("GO:0060333", allRes_AIC$GO.ID),]
allRes_T_Cell_AIC <- allRes_AIC[grep("GO:0050852", allRes_AIC$GO.ID),] 
allRes_Pos_T_Cell_AIC  <- allRes_AIC[grep("GO:0050862", allRes_AIC$GO.ID),] 
allRes_AIC_ERK <- allRes_AIC[grep("GO:0070374", allRes_AIC$GO.ID),]
allRes_Intrinsic_Apoptotic_AIC <- allRes_AIC[grep("GO:0042771", allRes_AIC$GO.ID),]
AIC_table <- rbind(allRes_IFN_gamma_AIC, allRes_T_Cell_AIC, allRes_Pos_T_Cell_AIC, allRes_AIC_ERK, allRes_Intrinsic_Apoptotic_AIC)
star_rating(allRes_IFN_gamma_AIC)
star_rating(allRes_T_Cell_AIC)
star_rating(allRes_Pos_T_Cell_AIC)
star_rating(allRes_AIC_ERK)
star_rating(allRes_Intrinsic_Apoptotic_AIC)

allRes_GI <- allRes_GI[,c(1,2,6)]
colnames(allRes_GI) <- c("GO_ID", "Term", "weightFisher")
allRes_GI$Hallmark <- "GI"
allRes_GI2 <- allRes_GI[,c(1,3,4)]
allRes_TA <- allRes_TA[,c(1,2,6)]
colnames(allRes_TA) <- c("GO_ID", "Term", "weightFisher")
allRes_TA$Hallmark <- "TA"
allRes_TA2 <- allRes_TA[,c(1,3,4)]
allRes_EA <- allRes_EA[,c(1,2,6)]
colnames(allRes_EA) <- c("GO_ID", "Term", "weightFisher")
allRes_EA$Hallmark <- "EA"
allRes_EA2 <- allRes_EA[,c(1,3,4)]
allRes_LOP <- allRes_LOP[,c(1,2,6)]
colnames(allRes_LOP) <- c("GO_ID", "Term", "weightFisher")
allRes_LOP$Hallmark <- "LOP"
allRes_LOP2 <- allRes_LOP[,c(1,3,4)]
allRes_DNS <- allRes_DNS[,c(1,2,6)]
colnames(allRes_DNS) <- c("GO_ID", "Term", "weightFisher")
allRes_DNS$Hallmark <- "DNS"
allRes_DNS2 <- allRes_DNS[,c(1,3,4)]
allRes_MD <- allRes_MD[,c(1,2,6)]
colnames(allRes_MD) <- c("GO_ID", "Term", "weightFisher")
allRes_MD$Hallmark <- "MD"
allRes_MD2 <- allRes_MD[,c(1,3,4)]
allRes_CS <- allRes_CS[,c(1,2,6)]
colnames(allRes_CS) <- c("GO_ID", "Term", "weightFisher")
allRes_CS$Hallmark <- "CS"
allRes_CS2 <- allRes_CS[,c(1,3,4)]
allRes_SCE <- allRes_SCE[,c(1,2,6)]
colnames(allRes_SCE) <- c("GO_ID", "Term", "weightFisher")
allRes_SCE$Hallmark <- "SCE"
allRes_SCE2 <- allRes_SCE[,c(1,3,4)]
allRes_AIC <- allRes_AIC[,c(1,2,6)]
colnames(allRes_AIC) <- c("GO_ID", "Term", "weightFisher")
allRes_AIC$Hallmark <- "AIC"
allRes_AIC2 <- allRes_AIC[,c(1,3,4)]

enrichment <- unique(as.data.frame(rbind(allRes_GI2, allRes_TA2,allRes_EA2,allRes_LOP2, allRes_DNS2, allRes_MD2, allRes_CS2, allRes_SCE2, allRes_AIC2)))
length(unique(enrichment$GO_ID))
goterms = unlist(Term(GOTERM))
goterms = data.frame(goterms)
goterms$GO <- rownames(goterms)
rownames(goterms) <- NULL
colnames(goterms) <- c("goterms", "GO_ID")
enrichment_update <- unique(merge(enrichment, goterms, by="GO_ID", all.x=TRUE, all.y =FALSE))
colnames(enrichment_update) <- c("GO_ID", "weightFisher", "Hallmark", "Term")
col2=colorRampPalette((brewer.pal(8, "YlGnBu")))(325)
#enrichment_update <- enrichment_update[grepl("GI|TA|EA|LOP|CS|MD|AIC|DNS|SCE", enrichment_update$Hallmark), ]

# Search for significant enrichment of terms related to genomic instability
new_AH_frame <- data.frame()
new_GI_frame <- AH_search(GI_list[1])
for (i in 2:length(GI_list)) {
  new_frame <- AH_search(GI_list[i])
  new_GI_frame <- rbind(new_GI_frame, new_frame)
}
enrichment_GI <- unique(new_GI_frame)
enrichment_GI$GO_ID <- NULL
enrichment_GI <- data.frame(pivot_wider(enrichment_GI, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_GI <- enrichment_GI[with(enrichment_GI, order(GI)), ]
levels(enrichment_GI$Term) <- c(levels(enrichment_GI$Term), "intrinsic apoptotic signaling pathway response to DNA damage (via p53)")
levels(enrichment_GI$Term) <- c(levels(enrichment_GI$Term), "DNA damage response (p53 class mediator, leading to p21 transcription)")
levels(enrichment_GI$Term) <- c(levels(enrichment_GI$Term), "DNA damage response (by p53 mediator leading to cell cycle arrest)")
levels(enrichment_GI$Term) <- c(levels(enrichment_GI$Term), "positive regulation of DNA damage response (via p53 class mediator)")
levels(enrichment_GI$Term) <- c(levels(enrichment_GI$Term), "negative regulation of DNA damage response (via p53 class mediator)")
enrichment_GI$Term[enrichment_GI$Term=="intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator"] <- "intrinsic apoptotic signaling pathway response to DNA damage (via p53)"
enrichment_GI$Term[enrichment_GI$Term=="DNA damage response, signal transduction by p53 class mediator resulting in transcription of p21 class mediator"] <- "DNA damage response (p53 class mediator, leading to p21 transcription)"
enrichment_GI$Term[enrichment_GI$Term=="DNA damage response, signal transduction by p53 class mediator resulting in cell cycle arrest"] <- "DNA damage response (by p53 mediator leading to cell cycle arrest)"
enrichment_GI$Term[enrichment_GI$Term=="positive regulation of DNA damage response, signal transduction by p53 class mediator"] <- "positive regulation of DNA damage response (via p53 class mediator)"
enrichment_GI$Term[enrichment_GI$Term=="negative regulation of DNA damage response, signal transduction by p53 class mediator"] <- "negative regulation of DNA damage response (via p53 class mediator)"

# Search for significant enrichment of terms related to telomere attrition
new_TA_frame <- AH_search(TA_list[1])
for (i in 2:length(TA_list)) {
  new_frame <- AH_search(TA_list[i])
  new_TA_frame <- rbind(new_TA_frame, new_frame)
}
enrichment_TA <- unique(new_TA_frame)
enrichment_TA$GO_ID <- NULL
enrichment_TA <- data.frame(pivot_wider(enrichment_TA, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_TA <- enrichment_TA[with(enrichment_TA, order(TA)), ]

# Search for significant enrichment of terms related to epigenetic alterations
new_EA_frame <- AH_search(EA_list[1])
for (i in 2:length(EA_list)) {
  new_frame <- AH_search(EA_list[i])
  new_EA_frame <- rbind(new_EA_frame, new_frame)
}
enrichment_EA <- unique(new_EA_frame)
enrichment_EA$GO_ID <- NULL
enrichment_EA <- data.frame(pivot_wider(enrichment_EA, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_EA <- enrichment_EA[with(enrichment_EA, order(EA)), ]
levels(enrichment_EA$Term) <- c(levels(enrichment_EA$Term), "positive regulation of pri-miRNA transcription")
levels(enrichment_EA$Term) <- c(levels(enrichment_EA$Term), "negative regulation of miRNA production (gene silencing)")
enrichment_EA$Term[enrichment_EA$Term=="positive regulation of pri-miRNA transcription from RNA polymerase II promoter"] <- "positive regulation of pri-miRNA transcription"
enrichment_EA$Term[enrichment_EA$Term=="negative regulation of production of miRNAs involved in gene silencing by miRNA"] <- "negative regulation of miRNA production (gene silencing)"

# Search for significant enrichment of terms related to loss of proteostasis
new_LOP_frame <- AH_search(LOP_list[1])
for (i in 2:length(LOP_list)) {
  new_frame <- AH_search(LOP_list[i])
  new_LOP_frame <- rbind(new_LOP_frame, new_frame)
}
enrichment_LOP <- unique(new_LOP_frame)
enrichment_LOP$GO_ID <- NULL
enrichment_LOP <- data.frame(pivot_wider(enrichment_LOP, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_LOP <- enrichment_LOP[with(enrichment_LOP, order(LOP)), ]
levels(enrichment_LOP$Term) <- c(levels(enrichment_LOP$Term), "intrinsic apoptotic signaling pathway in response to ERS")
levels(enrichment_LOP$Term) <- c(levels(enrichment_LOP$Term), "negative regulation of ERS-induced intrinsic apoptotic signaling pathway")
enrichment_LOP$Term[enrichment_LOP$Term=="intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress"] <- "intrinsic apoptotic signaling pathway in response to ERS"
enrichment_LOP$Term[enrichment_LOP$Term=="negative regulation of endoplasmic reticulum stress-induced intrinsic apoptotic signaling pathway"] <- "negative regulation of ERS-induced intrinsic apoptotic signaling pathway"

# Search for significant enrichment of terms related to dergeulated nutrient sensing
new_DNS_frame <- AH_search(DNS_list[1])
for (i in 2:length(DNS_list)) {
  new_frame <- AH_search(DNS_list[i])
  new_DNS_frame <- rbind(new_DNS_frame, new_frame)
}
enrichment_DNS <- unique(new_DNS_frame)
enrichment_DNS$GO_ID <- NULL
enrichment_DNS <- data.frame(pivot_wider(enrichment_DNS, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_DNS <- enrichment_DNS[with(enrichment_DNS, order(DNS)), ]
levels(enrichment_DNS$Term) <- c(levels(enrichment_DNS$Term), "regulation of insulin secretion (glucose stimulus)")
levels(enrichment_DNS$Term) <- c(levels(enrichment_DNS$Term), "negative regulation of insulin secretion (glucose stimulus)")
enrichment_DNS$Term[enrichment_DNS$Term=="regulation of insulin secretion involved in cellular response to glucose stimulus"] <- "regulation of insulin secretion (glucose stimulus)"
enrichment_DNS$Term[enrichment_DNS$Term=="negative regulation of insulin secretion involved in cellular response to glucose stimulus"] <- "negative regulation of insulin secretion (glucose stimulus)"

# Search for significant enrichment of terms related to mitochondrial dysfunction
new_MD_frame <- AH_search(MD_list[1])
for (i in 2:length(MD_list)) {
  new_frame <- AH_search(MD_list[i])
  new_MD_frame <- rbind(new_MD_frame, new_frame)
}
enrichment_MD <- unique(new_MD_frame)
enrichment_MD$GO_ID <- NULL
enrichment_MD <- data.frame(pivot_wider(enrichment_MD, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_MD <- enrichment_MD[with(enrichment_MD, order(MD)), ]
levels(enrichment_MD$Term) <- c(levels(enrichment_MD$Term), "protein insertion into mitochondrial membrane (apoptosis)")
levels(enrichment_MD$Term) <- c(levels(enrichment_MD$Term), "positive regulation of protein insertion into mitochondrial membrane (apoptosis)")
enrichment_MD$Term[enrichment_MD$Term=="protein insertion into mitochondrial membrane involved in apoptotic signaling pathway"] <- "protein insertion into mitochondrial membrane (apoptosis)"
enrichment_MD$Term[enrichment_MD$Term=="positive regulation of protein insertion into mitochondrial membrane involved in apoptotic signaling pathway"] <- "positive regulation of protein insertion into mitochondrial membrane (apoptosis)"

# Search for significant enrichment of terms related to cellular senescence
new_CS_frame <- AH_search(CS_list[1])
for (i in 2:length(CS_list)) {
  new_frame <- AH_search(CS_list[i])
  new_CS_frame <- rbind(new_CS_frame, new_frame)
}
enrichment_CS <- unique(new_CS_frame)
enrichment_CS$GO_ID <- NULL
enrichment_CS <- data.frame(pivot_wider(enrichment_CS, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_CS <- enrichment_CS[with(enrichment_CS, order(CS)), ]

# Search for significant enrichment of terms related to stem cell exhaustion
SCE_list <- list("stem cell", "stem cells", "cells stems", "progenitor cell", "progenitor cells", "cell progenitor", "cell progenitors")
new_SCE_frame <- AH_search(SCE_list[1])
for (i in 2:length(SCE_list)) {
  new_frame <- AH_search(SCE_list[i])
  new_SCE_frame <- rbind(new_SCE_frame, new_frame)
}
enrichment_SCE <- unique(new_SCE_frame)
enrichment_SCE$GO_ID <- NULL
enrichment_SCE <- data.frame(pivot_wider(enrichment_SCE, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_SCE <- enrichment_SCE[with(enrichment_SCE, order(SCE)), ]

# Search for significant enrichment of terms related to altered intercellular communication
new_AIC_frame <- AH_search(AIC_list[1])
for (i in 2:length(AIC_list)) {
  new_frame <- AH_search(AIC_list[i])
  new_AIC_frame <- rbind(new_AIC_frame, new_frame)
}
enrichment_AIC <- unique(new_AIC_frame)
enrichment_AIC$GO_ID <- NULL
enrichment_AIC <- data.frame(pivot_wider(enrichment_AIC, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_AIC <- enrichment_AIC[with(enrichment_AIC, order(AIC)), ]
levels(enrichment_AIC$Term) <- c(levels(enrichment_AIC$Term), "inflammatory response wound healing")
enrichment_AIC$Term[enrichment_AIC$Term=="connective tissue replacement involved in inflammatory response wound healing"] <- "inflammatory response wound healing"

# Search for significant enrichment of terms related to pathway/cascade
enrichment_Pathway <- enrichment_update[grepl("pathway|cascade", enrichment_update$Term, ignore.case = TRUE),]
enrichment_Pathway$GO_ID <- NULL
enrichment_Pathway <- data.frame(pivot_wider(enrichment_Pathway, names_from = "Hallmark", values_from ="weightFisher"))
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(GI)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(TA)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(EA)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(LOP)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(DNS)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(SCE)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(CS)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(AIC)), ]
enrichment_Pathway <- enrichment_Pathway[with(enrichment_Pathway, order(MD)), ]
enrichment_Pathway <- enrichment_Pathway[!is.na(enrichment_Pathway$Term),]
enrichment_Pathway <- drop_na(enrichment_Pathway)
levels(enrichment_Pathway$Term) <- c(levels(enrichment_Pathway$Term), "intrinsic apoptotic signaling pathway (DNA damage response)")
enrichment_Pathway$Term[enrichment_Pathway$Term=="intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator"] <- "intrinsic apoptotic signaling pathway (DNA damage response)"

# Heatmaps of primary aging hallmarks
enrichment_GI$AH <- rep("GI", 12)
enrichment_TA$AH <- rep("TA", 4)
enrichment_EA$AH <- rep("EA", 9)
enrichment_LOP$AH <- rep("LOP", 10)
primary_AH <- rbind(enrichment_GI, enrichment_TA, enrichment_EA, enrichment_LOP)
plot_heatmap(primary_AH, FALSE, 1, col_list = list("AH GO Term" = c("GI" = "black", "TA" = "#D3D3D3", "EA" = "#FEE23e", "LOP" = "#D21F3C")))

# Heatmaps of antagonistic aging hallmarks
enrichment_DNS$AH <- rep("DNS", 19)
enrichment_CS$AH <- rep("CS", 5)
enrichment_MD$AH <- rep("MD", 13) 
antagonistic_AH <- rbind(enrichment_DNS, enrichment_MD, enrichment_CS)
plot_heatmap(antagonistic_AH, FALSE, 1, col_list = list("AH GO Term" = c("DNS" = "#6A1B9A", "MD" = "#FF9800", "CS" = "#FF4081")))

# Heatmaps of integrative aging hallmarks
enrichment_AIC$AH <- rep("AIC", 24)
enrichment_SCE$AH <- rep("SCE", 9)
integrative_AH <- rbind(enrichment_SCE, enrichment_AIC)
plot_heatmap(integrative_AH, FALSE, 1, col_list = list("AH GO Term" = c("SCE" = "green", "AIC" = "blue")))

# Heatmaps of pathways/ cascades
enrichment_Pathway$AH <- rep("Pathway GO Term", 5)
plot_heatmap_pathway(enrichment_Pathway, TRUE, -2, col_list = list("Pathway GO Term" = c("Pathway GO Term" = "limegreen")))

# Creating the circos plots
# For the Ras-ERK circos plot: GO ontology term - GO: 0070374
ochiai_frame_ARDs <- ochiai_frame_ARDs_copy
ochiai_frame_ARDs <- ochiai_frame_ARDs[,c(2,3)]
my_ERK <- genesInTerm(myGOdata_cs, "GO:0070374")
my_ERK <- data.frame(my_ERK)
colnames(my_ERK) <- c("stringId")
gi_ERK <- merge(GI_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
ta_ERK <- merge(TA_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
ea_ERK <- merge(EA_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
lop_ERK <- merge(LOP_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
dns_ERK <- merge(DNS_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
md_ERK <- merge(MD_string, my_ERK, by = "stringId", all.x =FALSE, all.y = FALSE)
cs_ERK <- merge(CS_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
aic_ERK <- merge(AIC_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
sce_ERK <- merge(SCE_string, my_ERK, by="stringId", all.x=FALSE, all.y=FALSE)
ERK_all <- unique(rbind(gi_ERK, ta_ERK, ea_ERK, lop_ERK, dns_ERK, md_ERK, cs_ERK, sce_ERK, aic_ERK))
ERK_all <- unique(data.frame(ERK_all$name_node, ERK_all$Symbol, ERK_all$stringId))
colnames(ERK_all) <- c("name_node", "Symbol", "stringId")
ERK_all <- merge(ERK_all, ochiai_frame_ARDs, by = "name_node", all.x = TRUE, all.y = FALSE)
paste("The number of symbols linked to the ERK pathway was", length(unique(ERK_all$Symbol)))
paste("The number of ARDs linked to the ERK pathway was", length(unique(ERK_all$name_node)))
paste("The number of ARDs linked to the ERK pathway was", length(unique(ERK_all$stringId)))
ERK_all <- ERK_all[,c(2,4)]
colnames(ERK_all) <- c("Symbol", "Name")
draw_circos(ERK_all)
length(unique(ERK_all$Symbol))
length(unique(ERK_all$Name))

# For the IFN-gamma circos plot: GO ontology term - GO:0060333
my_IFN_genes <- genesInTerm(myGOdata_cs, "GO:0060333")
my_IFN_genes <- data.frame(my_IFN_genes)
colnames(my_IFN_genes) <- c("stringId")
gi_IFN <- merge(GI_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
ta_IFN <- merge(TA_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
ea_IFN <- merge(EA_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
lop_IFN <- merge(LOP_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
dns_IFN <- merge(DNS_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
md_IFN <- merge(MD_string, my_IFN_genes, by = "stringId", all.x =FALSE, all.y = FALSE)
cs_IFN <- merge(CS_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
aic_IFN <- merge(AIC_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
sce_IFN <- merge(SCE_string, my_IFN_genes, by="stringId", all.x=FALSE, all.y=FALSE)
IFN_all <- unique(rbind(gi_IFN, ta_IFN, ea_IFN, lop_IFN, dns_IFN, md_IFN, cs_IFN, sce_IFN, aic_IFN))
IFN_all <- unique(data.frame(IFN_all$name_node, IFN_all$Symbol, IFN_all$stringId))
colnames(IFN_all) <- c("name_node", "Symbol", "stringId")
IFN_all <- merge(IFN_all, ochiai_frame_ARDs, by = "name_node", all.x = TRUE, all.y = FALSE)
paste("The number of symbols linked to the IFN pathway was", length(unique(IFN_all$Symbol)))
paste("The number of ARDs linked to the IFN pathway was", length(unique(IFN_all$Name)))
paste("The number of ARDs linked to the IFN pathway was", length(unique(IFN_all$stringId)))
IFN_all <- IFN_all[,c(2,4)]
draw_circos(IFN_all)
length(unique(IFN_all$Symbol))

# For the T cell circos plot: GO ontology term - GO:0050852
my_T_cell_genes <- genesInTerm(myGOdata_cs, "GO:0050852")
my_T_cell_genes <- data.frame(my_T_cell_genes)
colnames(my_T_cell_genes) <- c("stringId")
gi_T_cell <- merge(GI_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
ta_T_cell <- merge(TA_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
ea_T_cell <- merge(EA_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
lop_T_cell <- merge(LOP_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
dns_T_cell <- merge(DNS_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
md_T_cell <- merge(MD_string, my_T_cell_genes, by = "stringId", all.x =FALSE, all.y = FALSE)
cs_T_cell <- merge(CS_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
aic_T_cell <- merge(AIC_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
sce_T_cell <- merge(SCE_string, my_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
T_cell_all <- unique(rbind(gi_T_cell, ta_T_cell, ea_T_cell, lop_T_cell, dns_T_cell, md_T_cell, cs_T_cell, sce_T_cell, aic_T_cell))
T_cell_all <- unique(data.frame(T_cell_all$name_node, T_cell_all$Symbol, T_cell_all$stringId))
colnames(T_cell_all) <- c("name_node", "Symbol", "stringId")
T_cell_all <- merge(T_cell_all, ochiai_frame_ARDs, by = "name_node", all.x = TRUE, all.y = FALSE)
paste("The number of symbols linked to the T-cell pathway was", length(unique(T_cell_all$Symbol)))
paste("The number of ARDs linked to the T-cell pathway was", length(unique(T_cell_all$name_node)))
T_cell_all <- T_cell_all[, c(2,4)]
draw_circos(T_cell_all)

# For the Positive T cell circos plot: GO ontology term - GO:0050862
my_pos_T_cell_genes <- genesInTerm(myGOdata_cs, "GO:0050862")
my_pos_T_cell_genes <- data.frame(my_pos_T_cell_genes)
colnames(my_pos_T_cell_genes) <- c("stringId")
gi_pos_T_cell <- merge(GI_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
ta_pos_T_cell <- merge(TA_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
ea_pos_T_cell <- merge(EA_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
lop_pos_T_cell <- merge(LOP_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
dns_pos_T_cell <- merge(DNS_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
md_pos_T_cell <- merge(MD_string, my_pos_T_cell_genes, by = "stringId", all.x =FALSE, all.y = FALSE)
cs_pos_T_cell <- merge(CS_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
aic_pos_T_cell <- merge(AIC_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
sce_pos_T_cell <- merge(SCE_string, my_pos_T_cell_genes, by="stringId", all.x=FALSE, all.y=FALSE)
Pos_T_cell_all <- unique(rbind(gi_pos_T_cell, ta_pos_T_cell, ea_pos_T_cell, lop_pos_T_cell, dns_pos_T_cell, 
                               md_pos_T_cell, cs_pos_T_cell, sce_pos_T_cell, aic_pos_T_cell))
Pos_T_cell_all <- unique(data.frame(Pos_T_cell_all$name_node, Pos_T_cell_all$Symbol, Pos_T_cell_all$stringId))
colnames(Pos_T_cell_all) <- c("name_node", "Symbol", "stringId")
Pos_T_cell_all <- merge(Pos_T_cell_all, ochiai_frame_ARDs, by = "name_node", all.x = TRUE, all.y = FALSE)
paste("The number of symbols linked to the positive T-cell pathway was", length(unique(Pos_T_cell_all$Symbol)))
paste("The number of ARDs linked to the positive T-cell pathway was", length(unique(Pos_T_cell_all$name_node)))
paste("The number of ARDs linked to the positive T-cell pathway was", length(unique(Pos_T_cell_all$stringId)))
draw_circos(Pos_T_cell_all)

# For the intrinsic apoptotic pathway circos plot: GO ontology term - GO:0042771
my_apoptotic <- genesInTerm(myGOdata_cs, "GO:0042771")
my_apoptotic <- data.frame(my_apoptotic)
colnames(my_apoptotic) <- c("stringId")
gi_apoptotic <- merge(GI_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
ta_apoptotic <- merge(TA_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
ea_apoptotic <- merge(EA_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
lop_apoptotic <- merge(LOP_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
dns_apoptotic <- merge(DNS_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
md_apoptotic <- merge(MD_string, my_apoptotic, by = "stringId", all.x =FALSE, all.y = FALSE)
cs_apoptotic <- merge(CS_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
aic_apoptotic <- merge(AIC_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
sce_apoptotic <- merge(SCE_string, my_apoptotic, by="stringId", all.x=FALSE, all.y=FALSE)
apoptotic_all <- unique(rbind(gi_apoptotic, ta_apoptotic, ea_apoptotic, lop_apoptotic, dns_apoptotic, 
                              md_apoptotic, cs_apoptotic, sce_apoptotic, aic_apoptotic))
apoptotic_all <- unique(data.frame(apoptotic_all$name_node, apoptotic_all$Symbol, apoptotic_all$stringId))
colnames(apoptotic_all) <- c("name_node", "Symbol", "stringId")
apoptotic_all <- merge(apoptotic_all, ochiai_frame_ARDs, by = "name_node", all.x = TRUE, all.y = FALSE)
paste("The number of symbols linked to the apoptotic pathway was", length(unique(apoptotic_all$Symbol)))
paste("The number of ARDs linked to the apoptotic pathway was", length(unique(apoptotic_all$name_node)))
paste("The number of ARDs linked to the apoptotic pathway was", length(unique(apoptotic_all$stringId)))
apoptotic_all <- apoptotic_all[,c(2,4)]
colnames(apoptotic_all) <- c("Symbol", "Name")
draw_circos(apoptotic_all)
