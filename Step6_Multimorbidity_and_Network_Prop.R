library("devtools")
library(igraph)
library(RColorBrewer)
library(BioNetSmooth)  #This is available from the Beyer Lab website and requires R 3.3.0 or less.
library("tidyverse")
library(dplyr)
library(plyr)
library(permute)
setwd("/Aging_hallmark_ARD_repository-main/")
# With thanks to the Prof. Andreas Beyer and Ronja Johnen for providing several of the functions for network analyses, which were adapted for this project, and the BioNetSmooth package 
# The shuffles were run in R 3.3.0 - set.seed()

#################################### Subnetworks from multimorbidity networks and Network Propagation ############################################# 
retrieve_node_ARDs <- function(frame) {
  # Import the multimorbidity networks
  # Retrieve the 184 of 207 ARD nodes in from the definitions paper with significant and partial correlation >= 0 
  # As the graphs are undirected: the Bonferroni correction is based on half the number of nodes
  v1 <- frame$DzA %in% ochiai_frame_ARDs$name_node
  frame <- frame[v1,]
  v2 <- frame$DzB %in% ochiai_frame_ARDs$name_node
  frame <- frame[v2,]
  frame$partialcor[frame$partialcor < 0] <- NA
  frame$partialcor.pval[frame$partialcor.pval > 0.05/(184*183/2)] <- NA
  frame <- frame[!is.na(frame$partialcor),]
  frame <- frame[!is.na(frame$partialcor.pval),]
  return(frame)
}

create_shuffled_frame <- function(frame_column, number_of_shuffles) {
  # Returns a shuffled data frame of ARDs mapped to different Ochiai coefficients
  shuffled_list <- data.frame()
  set.seed(12345)
  shuffled_list[1:184,1] <- ochiai_frame_ARDs$name_node
  colnames(shuffled_list) <- c('name_node')
  for (i in seq(1, number_of_shuffles)) {
    shuffled_list[1:184,i+1] <- frame_column[shuffle(frame_column)]
  }
  return(shuffled_list)
}

normalized_adj_matrix <- function(frame, ochiai_frame){ 
  # Produce a normalized adjacency matrix
  netmap_act_adjacency1 <- network_mapping_weighted(network=frame[,1:2], weights=frame[,4], expr_mat=ochiai_frame, type="adjacency", merge.by="name_node", global=T)
  names(netmap_act_adjacency1) <- c("disease_network", "aging_hallmark_OC", "disease_names")
  Adjacency_Matrix_weighted <- netmap_act_adjacency1$disease_network
  degrees <- colSums(as.matrix(Adjacency_Matrix_weighted))
  Adjacency_Matrix_weighted_normalized <- lapply(1:length(degrees), function(x) Adjacency_Matrix_weighted[x,]/degrees[x])
  Adjacency_Matrix_weighted_normalized <- do.call(rbind, Adjacency_Matrix_weighted_normalized)
  AdjMatCALIBERDisease <- data.matrix(Adjacency_Matrix_weighted_normalized)
  colnames(AdjMatCALIBERDisease) <- netmap_act_adjacency1$disease_names
  rownames(AdjMatCALIBERDisease) <- netmap_act_adjacency1$disease_names
  networkCALIBERDisease <- graph_from_adjacency_matrix(AdjMatCALIBERDisease,weighted=TRUE, mode="undirected")
  return(networkCALIBERDisease)
}

ochiai_frame_labels <- function(b_frame, hallmark_name){
  # Identify the top 30 ranked ARDs associated with each aging hallmark
  hallmark_name = hallmark_name + 1
  b_frame <- b_frame[order(-b_frame[,hallmark_name]),]
  b_frame <- b_frame[1:30,1:2]
  b_frame <- b_frame[order(b_frame$name_node),]
  b_frame$Name <- str_replace(b_frame$Name, " ", "\n")
  b_frame$Name <- str_replace(b_frame$Name, " ", "\n")
  return(b_frame$Name) 
}

order_frame <- function(a_frame, a_graph, AH_number, AH_title, Age_title, no_vertices, layout2) {
  # Plot a subnetwork of the top 30 nodes which is colour coded with red (1st to 10th), orange (11th to 20th), yellow (21st to 30th)
  a_graph_frame <- data.frame(vertex_attr(a_graph))
  a_frame1 <- a_frame[a_frame[,1] %in% a_graph_frame$name,]
  a <- a_frame1[order(-a_frame1[,AH_number]),]
  a1 <- as.vector(a[1:10,1])
  a2 <- as.vector(a[11:20,1])
  a3 <- as.vector(a[21:30,1])
  V(a_graph)$color <- "grey"
  V(a_graph)[a1]$color <- "red"
  V(a_graph)[a2]$color <- "orange"
  V(a_graph)[a3]$color <- "yellow"
  z1 <- edge_density(a_graph)
  a4 <- as.vector(a[1:no_vertices,1])
  top30subgraph <- induced_subgraph(a_graph, a4, impl = "create_from_scratch")
  g <- plot.igraph(top30subgraph, layout=layout2, vertex.size= 18, vertex.label.font=1.7, vertex.label.cex=0.8, 
                   vertex.label.color="black", edge.width=E(top30subgraph)$weight*25, edge.color="black")
  z2 <- edge_density(top30subgraph)
  return(top30subgraph)
}

order_frame2 <- function(a_frame, a_graph, b_frame, AH_number, AH_title, Age_title, no_vertices, layout2) {
  # Plot a subnetwork of the top 30 nodes which is colour coded with red (1st to 10th), orange (11th to 20th), yellow (21st to 30th)
  # Additional b_frame argument to change node names
  a_graph_frame <- data.frame(vertex_attr(a_graph))
  a_frame1 <- a_frame[a_frame[,1] %in% a_graph_frame$name,]
  a <- a_frame1[order(-a_frame1[,AH_number]),]
  a1 <- as.vector(a[1:10,1])
  a2 <- as.vector(a[11:20,1])
  a3 <- as.vector(a[21:30,1])
  V(a_graph)$color <- "grey"
  V(a_graph)[a1]$color <- "red"
  V(a_graph)[a2]$color <- "orange"
  V(a_graph)[a3]$color <- "yellow"
  z1 <- edge_density(a_graph)
  a4 <- as.vector(a[1:no_vertices,1])
  top30subgraph <- induced_subgraph(a_graph, a4, impl = "create_from_scratch")
  V(top30subgraph)$name <- ochiai_frame_labels(b_frame, AH_number)
  g <- plot.igraph(top30subgraph, layout=layout2, vertex.size= 18, vertex.label.font=1.7, vertex.label.cex=0.8, 
                   vertex.label.color="black", edge.width=E(top30subgraph)$weight*25, edge.color="black")
  z2 <- edge_density(top30subgraph)
  return(top30subgraph)
}

retrieve_edge_density <- function(scores_frame50, scores_frame60, scores_frame70, scores_frame80, n){
  # Returns the network density for every top 30 ARD subnetwork across all age categories
  data_m = data.frame()
  for (i in c(2,3,4,5,6,7,8,9,10)){
    z1 = data.frame(edge_density(order_frame(scores_frame50, g_50, i, "all", "age 50 to 59", n, layout2 = layout_with_fr)))
    data_m[1,1] = "age_50_to_59"
    data_m[i,1] = z1}
  for (i in c(2,3,4,5,6,7,8,9,10)){
    z2 = data.frame(edge_density(order_frame(scores_frame60, g_60, i, "all", "age 60 to 69", n, layout2 = layout_with_fr)))
    data_m[1,2] = "age_60_to_69"
    data_m[i,2] = z2}
  for (i in c(2,3,4,5,6,7,8,9,10)){
    z3 = data.frame(edge_density(order_frame(scores_frame70, g_70, i, "all", "age 70 to 79", n, layout2 = layout_with_fr)))
    data_m[1,3] = "age_70_to_79"
    data_m[i,3] = z3}
  for (i in c(2,3,4,5,6,7,8,9,10)){
    z4 = data.frame(edge_density(order_frame(scores_frame80, g_80, i, "all", "age 80+", n, layout2 = layout_with_fr)))
    data_m[1,4] = "age_80_to_89"
    data_m[i,4] = z4}
  colnames(data_m) = c(data_m[1,])
  data_m = data_m[2:10,]
  for (i in seq(1,4,1)){
    data_m[,i] = as.numeric(data_m[,i])}
  data_m = t(data_m)
  data_m = data.frame(data_m)
  return(data_m)
}

average_edge_density <- function(a_frame, a_graph, rep_number, no_vertices) {
  # Returns the network density for subnetworks of the top 30 nodes
  a_graph_frame <- data.frame(vertex_attr(a_graph))
  a_frame1 <- a_frame[a_frame[,1] %in% a_graph_frame$name,]
  a <- a_frame1[order(-a_frame1[,rep_number]),]
  a4 <- as.vector(a[1:no_vertices,1])
  top30subgraph <- induced_subgraph(a_graph, a4, impl = "create_from_scratch")
  z2 <- edge_density(top30subgraph)
  return(z2)
}

shuffling_frame <- function(shuffled_data50, shuffled_data60, shuffled_data70, shuffled_data80, n, a){
  # Retrieves the edge density for every set of top 30 ARDs nodes in the shuffled data and age category
  shuffled_frame = data.frame()
  for (i in c(seq(2,a))){
    z1 = average_edge_density(shuffled_data50, g_50, i, n)
    shuffled_frame[i-1,1] = z1}
  for (i in c(seq(2,a))){
    z2 = average_edge_density(shuffled_data60, g_60, i, n)
    shuffled_frame[i-1,2] = z2}
  for (i in c(seq(2,a))){
    z3 = average_edge_density(shuffled_data70, g_70, i, n)
    shuffled_frame[i-1,3] = z3}
  for (i in c(seq(2,a))){
    z4 = average_edge_density(shuffled_data80, g_80, i, n)
    shuffled_frame[i-1,4] = z4}
  colnames(shuffled_frame) <- c("50_to_59", "60_to_69", "70_to_79", "80_to_89")
  mean_n_shuffled <- t(data.frame(lapply(shuffled_frame, mean, 2)))
  return(shuffled_frame)
}

network_smoothed_values <- function(frame, ochiai_frame){
  # Scores are smoothed out over the network and the posterior scores of the nodes are returned in a data frame
  netmap_act_adjacency1 <- network_mapping_weighted(network=frame[,1:2], weights=frame[,4], expr_mat=ochiai_frame, type="adjacency", merge.by="name_node", global=T)
  names(netmap_act_adjacency1) <- c("disease_network", "aging_hallmark_OC", "disease_names")
  Adjacency_Matrix_weighted <- netmap_act_adjacency1$disease_network
  degrees <- colSums(as.matrix(Adjacency_Matrix_weighted))
  Adjacency_Matrix_weighted_normalized <- lapply(1:length(degrees), function(x) Adjacency_Matrix_weighted[x,]/degrees[x])
  Adjacency_Matrix_weighted_normalized <- do.call(rbind, Adjacency_Matrix_weighted_normalized)
  rowSums(Adjacency_Matrix_weighted_normalized)
  netmap_act_adjacency1$disease_network <- Adjacency_Matrix_weighted_normalized
  netmap_0.5_adjacency <- netmap_act_adjacency1
  netmap_0.5_adjacency$SmoothMats <- lapply(1:30, function(x) network_smoothing(net= netmap_0.5_adjacency$disease_network, mat_intensities=netmap_0.5_adjacency$aging_hallmark_OC,
                                                                                conditions = colnames(netmap_0.5_adjacency$aging_hallmark_OC), iter=x, alpha=0.5,network_type = "adjacency"))
  norm(netmap_0.5_adjacency$SmoothMats[[24]]- netmap_0.5_adjacency$SmoothMats[[25]])
  netmap_0.5_adjacency$ConvergedSmoothMat <- netmap_0.5_adjacency$SmoothMats[[25]]
  a <- data.frame(netmap_0.5_adjacency$disease_names, netmap_0.5_adjacency$ConvergedSmoothMat)
  return(a)
}

# Import the data
ochiai_frame_ARDs <- read.csv('data/aging_hallmark_subnetworks/ochiai_frame_0321.csv', header=TRUE)
ochiai_frame_ARDs$X.1 <- NULL
multiple_abbreviations <- read.csv("data/aging_hallmark_subnetworks/mappings_ARD_abbreviation.csv", header=TRUE)
ochiai_frame_ARDs <- merge(ochiai_frame_ARDs, multiple_abbreviations, by = "X")
ochiai_frame_ARDs <- ochiai_frame_ARDs[,c(1, 48, 49, 2:47)]
ochiai_frame_ARDs3 <- ochiai_frame_ARDs[, c("name_node", "Name", "GI_multiple", "TA_multiple", "EA_multiple", "LOP_multiple", "DNS_multiple", 
                                            "MD_multiple", "SCE_multiple", "CS_multiple", "AIC_multiple")]
ochiai_frame_ARDs <- ochiai_frame_ARDs[, c("name_node", "GI_multiple", "TA_multiple", "EA_multiple", "LOP_multiple", "DNS_multiple", 
                                           "MD_multiple", "SCE_multiple", "CS_multiple", "AIC_multiple")]
colnames(ochiai_frame_ARDs) <- c("name_node", "Disease_Ochiai_GI", "Disease_Ochiai_TA", "Disease_Ochiai_EA", 
                                 "Disease_Ochiai_LOP", "Disease_Ochiai_DNS", "Disease_Ochiai_MD", 
                                 "Disease_Ochiai_SCE", "Disease_Ochiai_CS", "Disease_Ochiai_AIC")

# Import the multimorbidity networks
age_50 <- retrieve_node_ARDs(read.table(file="multimorbidity_networks/link_age_50.dms"))
age_60 <- retrieve_node_ARDs(read.table(file="multimorbidity_networks/link_age_60.dms"))
age_70 <- retrieve_node_ARDs(read.table(file="multimorbidity_networks/link_age_70.dms"))
age_80 <- retrieve_node_ARDs(read.table(file="multimorbidity_networks/link_age_80.dms"))

# Normalized Adjaceny Matrix
g_50 <- normalized_adj_matrix(age_50, ochiai_frame_ARDs)
g_60 <- normalized_adj_matrix(age_60, ochiai_frame_ARDs)
g_70 <- normalized_adj_matrix(age_70, ochiai_frame_ARDs)
g_80 <- normalized_adj_matrix(age_80, ochiai_frame_ARDs)
paste("The edge density of the network for 50-59 years: ", edge_density(g_50))
paste("There are", gorder(g_50), "nodes in the network for 50-59 years")
paste("The edge density of the network for 60-69 years: ", edge_density(g_60))
paste("There are", gorder(g_60), "nodes in the network for 60-69 years")
paste("The edge density of the network for 70-79 years: ", edge_density(g_70))
paste("There are", gorder(g_70), "nodes in the network for 70-79 years")
paste("The edge density of the network for 80+ years: ", edge_density(g_80))
paste("There are", gorder(g_80), "nodes in the network for 80+ years")

# Network density based on the top 30 ARDs derived using the Ochiai coefficient - the true network density for GI, TA, EA, LOP, DNS, MD, CS, SCE, and AIC hallmarks
m_originals <- retrieve_edge_density(ochiai_frame_ARDs, ochiai_frame_ARDs, ochiai_frame_ARDs, ochiai_frame_ARDs, 30)
colnames(m_originals) <- colnames(ochiai_frame_ARDs)[2:10]

# Permutation test
# GI shuffles
shuffled_GI_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_GI, 20001)
write.csv(shuffled_GI_20000, "data/aging_hallmark_subnetworks/shuffled_GI_20000_ochiai_new1.csv")
shuffled_GI_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_GI_20000_ochiai_new1.csv")
shuffled_GI_20000$X <- NULL

# TA shuffles
shuffled_TA_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_TA, 20001)
write.csv(shuffled_TA_20000, "data/aging_hallmark_subnetworks/shuffled_TA_20000_ochiai_new1.csv")
shuffled_TA_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_TA_20000_ochiai_new1.csv")
shuffled_TA_20000$X <- NULL

# EA shuffles
shuffled_EA_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_EA, 20001)
write.csv(shuffled_EA_20000, "data/aging_hallmark_subnetworks/shuffled_EA_20000_ochiai_new1.csv")
shuffled_EA_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_EA_20000_ochiai_new1.csv")
shuffled_EA_20000$X <- NULL

# LOP shuffles
shuffled_LOP_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_LOP, 20001)
write.csv(shuffled_LOP_20000, "data/aging_hallmark_subnetworks/shuffled_LOP_20000_ochiai_new1.csv")
shuffled_LOP_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_LOP_20000_ochiai_new1.csv")
shuffled_LOP_20000$X <- NULL

# DNS shuffles
shuffled_DNS_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_DNS, 20001)
write.csv(shuffled_DNS_20000, "data/aging_hallmark_subnetworks/shuffled_DNS_20000_ochiai_new1.csv")
shuffled_DNS_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_DNS_20000_ochiai_new1.csv")
shuffled_DNS_20000$X <- NULL

# MD shuffles
shuffled_MD_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_MD, 20001)
write.csv(shuffled_MD_20000, "data/aging_hallmark_subnetworks/shuffled_MD_20000_ochiai_new1.csv")
shuffled_MD_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_MD_20000_ochiai_new1.csv")
shuffled_MD_20000$X <- NULL

# CS shuffles
shuffled_CS_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_CS, 20001)
write.csv(shuffled_CS_20000, "data/aging_hallmark_subnetworks/shuffled_CS_20000_ochiai_new1.csv")
shuffled_CS_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_CS_20000_ochiai_new1.csv") 
shuffled_CS_20000$X <- NULL

# SCE shuffles
shuffled_SCE_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_SCE, 20001)
write.csv(shuffled_SCE_20000, "data/aging_hallmark_subnetworks/shuffled_SCE_20000_ochiai_new1.csv")
shuffled_SCE_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_SCE_20000_ochiai_new1.csv")
shuffled_SCE_20000$X <- NULL

# AIC shuffles
shuffled_AIC_20000 <- create_shuffled_frame(ochiai_frame_ARDs$Disease_Ochiai_AIC, 20001)
write.csv(shuffled_AIC_20000, "data/aging_hallmark_subnetworks/shuffled_AIC_20000_ochiai_new1.csv")
shuffled_AIC_20000 <- read.csv("data/aging_hallmark_subnetworks/shuffled_AIC_20000_ochiai_new1.csv")
shuffled_AIC_20000$X <- NULL

# Permutation tests- This calculates the network density for each of the 20,000 shuffles
# GI
mean_n_GI_shuffled <- shuffling_frame(shuffled_GI_20000, shuffled_GI_20000, shuffled_GI_20000, shuffled_GI_20000, 30, a = 20001)
saveRDS(mean_n_GI_shuffled, "shuffled_ochiai/mean_n_GI_shuffled_20000_prior_new1.rds")
mean_n_GI_shuffled <- readRDS("shuffled_ochiai/mean_n_GI_shuffled_20000_prior_new1.rds")

# TA
mean_n_TA_shuffled <- shuffling_frame(shuffled_TA_20000, shuffled_TA_20000, shuffled_TA_20000, shuffled_TA_20000, 30, a = 20001)
saveRDS(mean_n_TA_shuffled, "shuffled_ochiai/mean_n_TA_shuffled_20000_prior_new1.rds")
mean_n_TA_shuffled <- readRDS("shuffled_ochiai/mean_n_TA_shuffled_20000_prior_new1.rds")

# EA
mean_n_EA_shuffled <- shuffling_frame(shuffled_EA_20000, shuffled_EA_20000, shuffled_EA_20000, shuffled_EA_20000, 30, a = 20001)
saveRDS(mean_n_EA_shuffled, "shuffled_ochiai/mean_n_EA_shuffled_20000_prior_new1.rds")
mean_n_EA_shuffled <- readRDS("shuffled_ochiai/mean_n_EA_shuffled_20000_prior_new1.rds")

# LOP
mean_n_LOP_shuffled <- shuffling_frame(shuffled_LOP_20000, shuffled_LOP_20000, shuffled_LOP_20000, shuffled_LOP_20000, 30, a = 20001)
saveRDS(mean_n_LOP_shuffled, "shuffled_ochiai/mean_n_LOP_shuffled_20000_prior_new1.rds")
mean_n_LOP_shuffled <- readRDS("shuffled_ochiai/mean_n_LOP_shuffled_20000_prior_new1.rds")

# DNS
mean_n_DNS_shuffled <- shuffling_frame(shuffled_DNS_20000, shuffled_DNS_20000, shuffled_DNS_20000, shuffled_DNS_20000, 30, a = 20001)
saveRDS(mean_n_DNS_shuffled, "shuffled_ochiai/mean_n_DNS_shuffled_20000_prior_new1.rds")
mean_n_DNS_shuffled <- readRDS("shuffled_ochiai/mean_n_DNS_shuffled_20000_prior_new1.rds")

# MD
mean_n_MD_shuffled <- shuffling_frame(shuffled_MD_20000, shuffled_MD_20000, shuffled_MD_20000, shuffled_MD_20000, 30, a = 20001)
saveRDS(mean_n_MD_shuffled, "shuffled_ochiai/mean_n_MD_shuffled_20000_prior_new1.rds")
mean_n_MD_shuffled <- readRDS("shuffled_ochiai/mean_n_MD_shuffled_20000_prior_new1.rds")

# CS
mean_n_CS_shuffled <- shuffling_frame(shuffled_CS_20000, shuffled_CS_20000, shuffled_CS_20000, shuffled_CS_20000, 30, a = 20001)
saveRDS(mean_n_CS_shuffled, "shuffled_ochiai/mean_n_CS_shuffled_20000_prior_new1.rds")
mean_n_CS_shuffled <- readRDS("shuffled_ochiai/mean_n_CS_shuffled_20000_prior_new1.rds")

# SCE
mean_n_SCE_shuffled <- shuffling_frame(shuffled_SCE_20000, shuffled_SCE_20000, shuffled_SCE_20000, shuffled_SCE_20000, 30, a = 20001)
saveRDS(mean_n_SCE_shuffled, "shuffled_ochiai/mean_n_SCE_shuffled_20000_prior_new1.rds")
mean_n_SCE_shuffled <- readRDS("shuffled_ochiai/mean_n_SCE_shuffled_20000_prior_new1.rds")

# AIC
mean_n_AIC_shuffled <- shuffling_frame(shuffled_AIC_20000, shuffled_AIC_20000, shuffled_AIC_20000, shuffled_AIC_20000, 30, a = 20001)
saveRDS(mean_n_AIC_shuffled, "shuffled_ochiai/mean_n_AIC_shuffled_20000_prior_new1.rds")
mean_n_AIC_shuffled <- readRDS("shuffled_ochiai/mean_n_AIC_shuffled_20000_prior_new1.rds")

# Deriving the p-values
GI_50 <- sum(mean_n_GI_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_GI[1])/20000
GI_60 <- sum(mean_n_GI_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_GI[2])/20000
GI_70 <- sum(mean_n_GI_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_GI[3])/20000
GI_80 <- sum(mean_n_GI_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_GI[4])/20000
p_GI <- c(GI_50, GI_60, GI_70, GI_80)
p_GI <- p.adjust(p_GI, method = "fdr", n = length(p_GI))
paste("The p values for the GI 50-59 year subnetworks are: ", p_GI[1])
paste("The p values for the GI 60-69 year subnetworks are: ", p_GI[2])
paste("The p values for the GI 70-79 year subnetworks are: ", p_GI[3])
paste("The p values for the GI 80-89 year subnetworks are: ", p_GI[4])

TA_50 <- sum(mean_n_TA_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_TA[1])/20000
TA_60 <- sum(mean_n_TA_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_TA[2])/20000
TA_70 <- sum(mean_n_TA_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_TA[3])/20000
TA_80 <- sum(mean_n_TA_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_TA[4])/20000
p_TA <- c(TA_50, TA_60, TA_70, TA_80)
p_TA <- p.adjust(p_TA, method = "fdr", n = length(p_TA))
paste("The p values for the TA 50-59 year subnetworks are: ", p_TA[1])
paste("The p values for the TA 60-69 year subnetworks are: ", p_TA[2])
paste("The p values for the TA 70-79 year subnetworks are: ", p_TA[3])
paste("The p values for the TA 80-89 year subnetworks are: ", p_TA[4])

EA_50 <- sum(mean_n_EA_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_EA[1])/20000
EA_60 <- sum(mean_n_EA_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_EA[2])/20000
EA_70 <- sum(mean_n_EA_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_EA[3])/20000
EA_80 <- sum(mean_n_EA_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_EA[4])/20000
p_EA <- c(EA_50, EA_60, EA_70, EA_80)
p_EA <- p.adjust(p_EA, method = "fdr", n = length(p_EA))
paste("The p values for the EA 50-59 year subnetworks are: ", p_EA[1])
paste("The p values for the EA 60-69 year subnetworks are: ", p_EA[2])
paste("The p values for the EA 70-79 year subnetworks are: ", p_EA[3])
paste("The p values for the EA 80-89 year subnetworks are: ", p_EA[4])

LOP_50 <- sum(mean_n_LOP_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_LOP[1])/20000
LOP_60 <- sum(mean_n_LOP_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_LOP[2])/20000
LOP_70 <- sum(mean_n_LOP_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_LOP[3])/20000
LOP_80 <- sum(mean_n_LOP_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_LOP[4])/20000
p_LOP <- c(LOP_50, LOP_60, LOP_70, LOP_80)
p_LOP <- p.adjust(p_LOP, method = "fdr", n = length(p_LOP))
paste("The p values for the LOP 50-59 year subnetworks are: ", p_LOP[1])
paste("The p values for the LOP 60-69 year subnetworks are: ", p_LOP[2])
paste("The p values for the LOP 70-79 year subnetworks are: ", p_LOP[3])
paste("The p values for the LOP 80-89 year subnetworks are: ", p_LOP[4])

DNS_50 <- sum(mean_n_DNS_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_DNS[1])/20000
DNS_60 <- sum(mean_n_DNS_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_DNS[2])/20000
DNS_70 <- sum(mean_n_DNS_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_DNS[3])/20000
DNS_80 <- sum(mean_n_DNS_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_DNS[4])/20000
p_DNS <- c(DNS_50,DNS_60,DNS_70,DNS_80)
p_DNS <- p.adjust(p_DNS, method = "fdr", n = length(p_DNS))
paste("The p values for the DNS 50-59 year subnetworks are: ", p_DNS[1])
paste("The p values for the DNS 60-69 year subnetworks are: ", p_DNS[2])
paste("The p values for the DNS 70-79 year subnetworks are: ", p_DNS[3])
paste("The p values for the DNS 80-89 year subnetworks are: ", p_DNS[4])

MD_50 <- sum(mean_n_MD_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_MD[1])/20000
MD_60 <- sum(mean_n_MD_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_MD[2])/20000
MD_70 <- sum(mean_n_MD_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_MD[3])/20000
MD_80 <- sum(mean_n_MD_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_MD[4])/20000
p_MD <- c(MD_50, MD_60, MD_70, MD_80)
p_MD <- p.adjust(p_MD, method = "fdr", n = length(p_MD))
paste("The p values for the MD 50-59 year subnetworks are: ", p_MD[1])
paste("The p values for the MD 60-69 year subnetworks are: ", p_MD[2])
paste("The p values for the MD 70-79 year subnetworks are: ", p_MD[3])
paste("The p values for the MD 80-89 year subnetworks are: ", p_MD[4])

CS_50 <- sum(mean_n_CS_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_CS[1])/20000
CS_60 <- sum(mean_n_CS_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_CS[2])/20000
CS_70 <- sum(mean_n_CS_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_CS[3])/20000
CS_80 <- sum(mean_n_CS_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_CS[4])/20000
p_CS <- c(CS_50,CS_60,CS_70, CS_80)
p_CS <- p.adjust(p_CS, method = "fdr", n = length(p_CS))
paste("The p values for the CS 50-59 year subnetworks are: ", p_CS[1])
paste("The p values for the CS 60-69 year subnetworks are: ", p_CS[2])
paste("The p values for the CS 70-79 year subnetworks are: ", p_CS[3])
paste("The p values for the CS 80-89 year subnetworks are: ", p_CS[4])

SCE_50 <- sum(mean_n_SCE_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_SCE[1])/20000
SCE_60 <- sum(mean_n_SCE_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_SCE[2])/20000
SCE_70 <- sum(mean_n_SCE_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_SCE[3])/20000
SCE_80 <- sum(mean_n_SCE_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_SCE[4])/20000
p_SCE <- c(SCE_50, SCE_60, SCE_70, SCE_80)
p_SCE <- p.adjust(p_SCE, method = "fdr", n = length(p_SCE))
paste("The p values for the SCE 50-59 year subnetworks are: ", p_SCE[1])
paste("The p values for the SCE 60-69 year subnetworks are: ", p_SCE[2])
paste("The p values for the SCE 70-79 year subnetworks are: ", p_SCE[3])
paste("The p values for the SCE 80-89 year subnetworks are: ", p_SCE[4])

AIC_50 <- sum(mean_n_AIC_shuffled$`50_to_60` >= m_originals$Disease_Ochiai_AIC[1])/20000
AIC_60 <- sum(mean_n_AIC_shuffled$`60_to_70` >= m_originals$Disease_Ochiai_AIC[2])/20000
AIC_70 <- sum(mean_n_AIC_shuffled$`70_to_80` >= m_originals$Disease_Ochiai_AIC[3])/20000
AIC_80 <- sum(mean_n_AIC_shuffled$`80_to_90` >= m_originals$Disease_Ochiai_AIC[4])/20000
p_AIC <- c(AIC_50,AIC_60,AIC_70,AIC_80)
p_AIC <- p.adjust(p_AIC, method = "fdr", n = length(p_AIC))
paste("The p values for the AIC 50-59 year subnetworks are: ", p_AIC[1])
paste("The p values for the AIC 60-69 year subnetworks are: ", p_AIC[2])
paste("The p values for the AIC 70-79 year subnetworks are: ", p_AIC[3])
paste("The p values for the AIC 80-89 year subnetworks are: ", p_AIC[4])

# Summarizing the p-values
m_originals <- round(m_originals, digits = 4)
m_originals1 <- t(m_originals)
m_originals_p_vals <- data.frame(p_GI, p_TA, p_EA, p_LOP, p_DNS, p_MD, p_CS, p_SCE, p_AIC)
m_originals_p_vals <- round(m_originals_p_vals, digits = 4)
m_originals_p_vals1 <- t(m_originals_p_vals)

# Plot the subgraphs for the top 30 ARDs in the 50 - 59 year age category
plot.new()
old.par <- par(mar = c(0, 0, 0, 0))
par(mfrow = c(1,1))
# GI
edge_GI_50 <- order_frame2(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 2, "GI", "Age 50 to 59 years:", 30, layout2=layout_with_fr)
# TA
edge_TA_50 <- order_frame2(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 3, "TA", "Age 50 to 59 years:", 30, layout2=layout_with_fr)
# EA
edge_EA_50 <- order_frame2(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 4, "EA", "Age 50 to 59 years:", 30, layout2=layout_with_fr)
# LOP
edge_LOP_50 <- order_frame2(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 5, "LOP", "Age 50 to 59 years:", 30, layout2=layout_with_fr)
# DNS
edge_DNS_50 <- order_frame2(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 6, "DNS", "Age 50 to 59 years:", 30, layout2=layout_with_graphopt)
# MD
edge_MD_50 <- order_frame(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 7, "MD", "Age 50 to 60 years:", 30, layout2=layout_with_graphopt)
# CS
edge_CS_50 <- order_frame(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 9, "CS", "Age 60 to 70 years:", 30, layout2=layout_with_fr)
# SCE
edge_SCE_50 <- order_frame(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 8, "SCE", "Age 50 to 60 years:", 30, layout2=layout_with_fr)
# AIC
edge_AIC_50 <- order_frame(ochiai_frame_ARDs, g_50, ochiai_frame_ARDs3, 10, "AIC", "Age 50 to 60 years:", 30, layout2=layout_with_fr)

# Network propagation - Posterior scores 
a_50 <- network_smoothed_values(age_50, ochiai_frame_ARDs)
a_60 <- network_smoothed_values(age_60, ochiai_frame_ARDs)
a_70 <- network_smoothed_values(age_70, ochiai_frame_ARDs)
a_80 <- network_smoothed_values(age_80, ochiai_frame_ARDs)

# Network density for the true subnetworks based on top 30 ARDs after network propagation
m = data.frame()
m <- retrieve_edge_density(a_50, a_60, a_70, a_80, 30)
colnames(m) <- colnames(ochiai_frame_ARDs)[2:10]

# Network density for shuffled data based on posterior scores
# GI
gi_a_50 <- network_smoothed_values(age_50, shuffled_GI_20000)
gi_a_60 <- network_smoothed_values(age_60, shuffled_GI_20000)
gi_a_80 <- network_smoothed_values(age_80, shuffled_GI_20000)
gi_a_70 <- network_smoothed_values(age_70, shuffled_GI_20000)
saveRDS(gi_a_50, "shuffled_ochiai/gi_a_50_ARD_20000_ochiai.rds")
saveRDS(gi_a_60, "shuffled_ochiai/gi_a_60_ARD_20000_ochiai.rds")
saveRDS(gi_a_70, "shuffled_ochiai/gi_a_70_ARD_20000_ochiai.rds")
saveRDS(gi_a_80, "shuffled_ochiai/gi_a_80_ARD_20000_ochiai.rds")
gi_a_50 <- readRDS("shuffled_ochiai/gi_a_50_ARD_20000_ochiai.rds")
gi_a_60 <- readRDS("shuffled_ochiai/gi_a_60_ARD_20000_ochiai.rds")
gi_a_70 <- readRDS("shuffled_ochiai/gi_a_70_ARD_20000_ochiai.rds")
gi_a_80 <- readRDS("shuffled_ochiai/gi_a_80_ARD_20000_ochiai.rds")
mean_n_GI <- shuffling_frame(gi_a_50, gi_a_60, gi_a_70, gi_a_80, 30, a = 20001)
saveRDS(mean_n_GI, "shuffled_ochiai/mean_n_GI_shuffled_posterior.rds")
mean_n_GI <- readRDS("shuffled_ochiai/mean_n_GI_shuffled_posterior.rds")

# TA
ta_a_50 <- network_smoothed_values(age_50, shuffled_TA_20000)
ta_a_60 <- network_smoothed_values(age_60, shuffled_TA_20000)
ta_a_70 <- network_smoothed_values(age_70, shuffled_TA_20000)
ta_a_80 <- network_smoothed_values(age_80, shuffled_TA_20000)
saveRDS(ta_a_50, "shuffled_ochiai/ta_a_50_ARD_20000_ochiai.rds")
saveRDS(ta_a_60, "shuffled_ochiai/ta_a_60_ARD_20000_ochiai.rds")
saveRDS(ta_a_70, "shuffled_ochiai/ta_a_70_ARD_20000_ochiai.rds")
saveRDS(ta_a_80, "shuffled_ochiai/ta_a_80_ARD_20000_ochiai.rds")
ta_a_50 <- readRDS("shuffled_ochiai/ta_a_50_ARD_20000_ochiai.rds")
ta_a_60 <- readRDS("shuffled_ochiai/ta_a_60_ARD_20000_ochiai.rds")
ta_a_70 <- readRDS("shuffled_ochiai/ta_a_70_ARD_20000_ochiai.rds")
ta_a_80 <- readRDS("shuffled_ochiai/ta_a_80_ARD_20000_ochiai.rds")
mean_n_TA <- shuffling_frame(ta_a_50, ta_a_60, ta_a_70, ta_a_80, 30, a = 20001)
saveRDS(mean_n_TA, "shuffled_ochiai/mean_n_TA_shuffled_posterior.rds")
mean_n_TA <- readRDS("shuffled_ochiai/mean_n_TA_shuffled_posterior.rds")

# EA
ea_a_50 <- network_smoothed_values(age_50, shuffled_EA_20000)
ea_a_60 <- network_smoothed_values(age_60, shuffled_EA_20000)
ea_a_70 <- network_smoothed_values(age_70, shuffled_EA_20000)
ea_a_80 <- network_smoothed_values(age_80, shuffled_EA_20000)
saveRDS(ea_a_50, "shuffled_ochiai/ea_a_50_ARD_20000_ochiai.rds")
saveRDS(ea_a_60, "shuffled_ochiai/ea_a_60_ARD_20000_ochiai.rds")
saveRDS(ea_a_70, "shuffled_ochiai/ea_a_70_ARD_20000_ochiai.rds")
saveRDS(ea_a_80, "shuffled_ochiai/ea_a_80_ARD_20000_ochiai.rds")
ea_a_50 <- readRDS("shuffled_ochiai/ea_a_50_ARD_20000_ochiai.rds")
ea_a_60 <- readRDS("shuffled_ochiai/ea_a_60_ARD_20000_ochiai.rds")
ea_a_70 <- readRDS("shuffled_ochiai/ea_a_70_ARD_20000_ochiai.rds")
ea_a_80 <- readRDS("shuffled_ochiai/ea_a_80_ARD_20000_ochiai.rds")
mean_n_EA <- shuffling_frame(ea_a_50, ea_a_60, ea_a_70, ea_a_80, 30, a = 20001)
saveRDS(mean_n_EA, "shuffled_ochiai/mean_n_EA_shuffled_posterior.rds")
mean_n_EA <- readRDS("shuffled_ochiai/mean_n_EA_shuffled_posterior.rds")

# LOP 
lop_a_50 <- network_smoothed_values(age_50, shuffled_LOP_20000)
lop_a_60 <- network_smoothed_values(age_60, shuffled_LOP_20000)
lop_a_70 <- network_smoothed_values(age_70, shuffled_LOP_20000)
lop_a_80 <- network_smoothed_values(age_80, shuffled_LOP_20000)
saveRDS(lop_a_50, "shuffled_ochiai/lop_a_50_ARD_20000_ochiai.rds")
saveRDS(lop_a_60, "shuffled_ochiai/lop_a_60_ARD_20000_ochiai.rds")
saveRDS(lop_a_70, "shuffled_ochiai/lop_a_70_ARD_20000_ochiai.rds")
saveRDS(lop_a_80, "shuffled_ochiai/lop_a_80_ARD_20000_ochiai.rds")
lop_a_50 <- readRDS("shuffled_ochiai/lop_a_50_ARD_20000_ochiai.rds")
lop_a_60 <- readRDS("shuffled_ochiai/lop_a_60_ARD_20000_ochiai.rds")
lop_a_70 <- readRDS("shuffled_ochiai/lop_a_70_ARD_20000_ochiai.rds")
lop_a_80 <- readRDS("shuffled_ochiai/lop_a_80_ARD_20000_ochiai.rds")
mean_n_LOP <- shuffling_frame(lop_a_50, lop_a_60, lop_a_70, lop_a_80, 30, a = 20001)
saveRDS(mean_n_LOP, "shuffled_ochiai/mean_n_LOP_shuffled_posterior.rds")
mean_n_LOP <- readRDS("shuffled_ochiai/mean_n_LOP_shuffled_posterior.rds")

# DNS
dns_a_50 <- network_smoothed_values(age_50, shuffled_DNS_20000)
dns_a_60 <- network_smoothed_values(age_60, shuffled_DNS_20000)
dns_a_70 <- network_smoothed_values(age_70, shuffled_DNS_20000)
dns_a_80 <- network_smoothed_values(age_80, shuffled_DNS_20000)
saveRDS(dns_a_50, "shuffled_ochiai/dns_a_50_ARD_20000_ochiai.rds")
saveRDS(dns_a_60, "shuffled_ochiai/dns_a_60_ARD_20000_ochiai.rds")
saveRDS(dns_a_70, "shuffled_ochiai/dns_a_70_ARD_20000_ochiai.rds")
saveRDS(dns_a_80, "shuffled_ochiai/dns_a_80_ARD_20000_ochiai.rds")
dns_a_50 <- readRDS("shuffled_ochiai/dns_a_50_ARD_20000_ochiai.rds")
dns_a_60 <- readRDS("shuffled_ochiai/dns_a_60_ARD_20000_ochiai.rds")
dns_a_70 <- readRDS("shuffled_ochiai/dns_a_70_ARD_20000_ochiai.rds")
dns_a_80 <- readRDS("shuffled_ochiai/dns_a_80_ARD_20000_ochiai.rds")
mean_n_DNS <- shuffling_frame(dns_a_50, dns_a_60, dns_a_70, dns_a_80, 30, a = 20001)
saveRDS(mean_n_DNS, "shuffled_ochiai/mean_n_DNS_shuffled_posterior.rds")
mean_n_DNS <- readRDS("shuffled_ochiai/mean_n_DNS_shuffled_posterior.rds")

# MD
md_a_50 <- network_smoothed_values(age_50, shuffled_MD_5000)
md_a_60 <- network_smoothed_values(age_60, shuffled_MD_5000)
md_a_70 <- network_smoothed_values(age_70, shuffled_MD_5000)
md_a_80 <- network_smoothed_values(age_80, shuffled_MD_5000)
saveRDS(md_a_50, "shuffled_ochiai/md_a_50_ARD_20000_ochiai.rds")
saveRDS(md_a_60, "shuffled_ochiai/md_a_60_ARD_20000_ochiai.rds")
saveRDS(md_a_70, "shuffled_ochiai/md_a_70_ARD_20000_ochiai.rds")
saveRDS(md_a_80, "shuffled_ochiai/md_a_80_ARD_20000_ochiai.rds")
md_a_50 <- readRDS("shuffled_ochiai/md_a_50_ARD_20000_ochiai.rds")
md_a_60 <- readRDS("shuffled_ochiai/md_a_60_ARD_20000_ochiai.rds")
md_a_70 <- readRDS("shuffled_ochiai/md_a_70_ARD_20000_ochiai.rds")
md_a_80 <- readRDS("shuffled_ochiai/md_a_80_ARD_20000_ochiai.rds")
mean_n_MD <- shuffling_frame(md_a_50, md_a_60, md_a_70, md_a_80, 30, a = 1001)
saveRDS(mean_n_MD, "shuffled_ochiai/mean_n_MD_shuffled_posterior.rds")
mean_n_MD <- readRDS("shuffled_ochiai/mean_n_MD_shuffled_posterior.rds")

# CS
cs_a_50 <- network_smoothed_values(age_50, shuffled_CS_20000)
cs_a_60 <- network_smoothed_values(age_60, shuffled_CS_20000)
cs_a_70 <- network_smoothed_values(age_70, shuffled_CS_20000)
cs_a_80 <- network_smoothed_values(age_80, shuffled_CS_20000)
saveRDS(cs_a_50, "shuffled_ochiai/cs_a_50_ARD_20000_ochiai.rds")
saveRDS(cs_a_60, "shuffled_ochiai/cs_a_60_ARD_20000_ochiai.rds")
saveRDS(cs_a_70, "shuffled_ochiai/cs_a_70_ARD_20000_ochiai.rds")
saveRDS(cs_a_80, "shuffled_ochiai/cs_a_80_ARD_20000_ochiai.rds")
cs_a_50 <- readRDS("shuffled_ochiai/cs_a_50_ARD_20000_ochiai.rds")
cs_a_60 <- readRDS("shuffled_ochiai/cs_a_60_ARD_20000_ochiai.rds")
cs_a_70 <- readRDS("shuffled_ochiai/cs_a_70_ARD_20000_ochiai.rds")
cs_a_80 <- readRDS("shuffled_ochiai/cs_a_80_ARD_20000_ochiai.rds")
mean_n_CS <- shuffling_frame(cs_a_50, cs_a_60, cs_a_70, cs_a_80, 30, a = 20001)
saveRDS(mean_n_CS, "shuffled_ochiai/mean_n_CS_shuffled_posterior.rds")
mean_n_CS <- readRDS("shuffled_ochiai/mean_n_CS_shuffled_posterior.rds")

# SCE 
sce_a_50 <- network_smoothed_values(age_50, shuffled_SCE_20000)
sce_a_60 <- network_smoothed_values(age_60, shuffled_SCE_20000)
sce_a_70 <- network_smoothed_values(age_70, shuffled_SCE_20000)
sce_a_80 <- network_smoothed_values(age_80, shuffled_SCE_20000)
saveRDS(sce_a_50, "shuffled_ochiai/sce_a_50_ARD_20000_ochiai.rds")
saveRDS(sce_a_60, "shuffled_ochiai/sce_a_60_ARD_20000_ochiai.rds")
saveRDS(sce_a_70, "shuffled_ochiai/sce_a_70_ARD_20000_ochiai.rds")
saveRDS(sce_a_80, "shuffled_ochiai/sce_a_80_ARD_20000_ochiai.rds")
sce_a_50 <- readRDS("shuffled_ochiai/sce_a_50_ARD_20000_ochiai.rds")
sce_a_60 <- readRDS("shuffled_ochiai/sce_a_60_ARD_20000_ochiai.rds")
sce_a_70 <- readRDS("shuffled_ochiai/sce_a_70_ARD_20000_ochiai.rds")
sce_a_80 <- readRDS("shuffled_ochiai/sce_a_80_ARD_20000_ochiai.rds")
mean_n_SCE <- shuffling_frame(sce_a_50, sce_a_60, sce_a_70, sce_a_80, 30, a = 20001)
saveRDS(mean_n_SCE, "shuffled_ochiai/mean_n_SCE_shuffled_posterior.rds")
mean_n_SCE <- readRDS("shuffled_ochiai/mean_n_SCE_shuffled_posterior.rds")

# AIC
aic_a_50 <- network_smoothed_values(age_50, shuffled_AIC_20000)
aic_a_60 <- network_smoothed_values(age_60, shuffled_AIC_20000)
aic_a_70 <- network_smoothed_values(age_70, shuffled_AIC_20000)
aic_a_80 <- network_smoothed_values(age_80, shuffled_AIC_20000)
saveRDS(aic_a_50, "shuffled_ochiai/aic_a_50_ARD_20000_ochiai.rds")
saveRDS(aic_a_60, "shuffled_ochiai/aic_a_60_ARD_20000_ochiai.rds")
saveRDS(aic_a_70, "shuffled_ochiai/aic_a_70_ARD_20000_ochiai.rds")
saveRDS(aic_a_80, "shuffled_ochiai/aic_a_80_ARD_20000_ochiai.rds")
aic_a_50 <- readRDS("shuffled_ochiai/aic_a_50_ARD_20000_ochiai.rds")
aic_a_60 <- readRDS("shuffled_ochiai/aic_a_60_ARD_20000_ochiai.rds")
aic_a_70 <- readRDS("shuffled_ochiai/aic_a_70_ARD_20000_ochiai.rds")
aic_a_80 <- readRDS("shuffled_ochiai/aic_a_80_ARD_20000_ochiai.rds")
mean_n_AIC <- shuffling_frame(aic_a_50, aic_a_60, aic_a_70, aic_a_80, 30, a = 20001)
saveRDS(mean_n_AIC, "shuffled_ochiai/mean_n_AIC_shuffled_posterior.rds")
mean_n_AIC <- readRDS("shuffled_ochiai/mean_n_AIC_shuffled_posterior.rds")

# Identifying significant subnetworks based on posterior score
GI_shuffled_50 <- sum(mean_n_GI$`50_to_60` >= m$Disease_Ochiai_GI[1])/20000
GI_shuffled_60 <- sum(mean_n_GI$`60_to_70` >= m$Disease_Ochiai_GI[2])/20000
GI_shuffled_70 <- sum(mean_n_GI$`70_to_80` >= m$Disease_Ochiai_GI[3])/20000
GI_shuffled_80 <- sum(mean_n_GI$`80_to_90` >= m$Disease_Ochiai_GI[4])/20000
p_GI_Post = c(GI_shuffled_50, GI_shuffled_60,GI_shuffled_70,GI_shuffled_80)
p_GI_Post = p.adjust(p_GI_Post, method = "fdr", n = length(p_GI_Post))
paste("The p values for the GI 50-59 year subnetworks based on posterior score are: ", p_GI_Post[1])
paste("The p values for the GI 60-69 year subnetworks based on posterior score are: ", p_GI_Post[2])
paste("The p values for the GI 70-79 year subnetworks based on posterior score are: ", p_GI_Post[3])
paste("The p values for the GI 80-89 year subnetworks based on posterior score are: ", p_GI_Post[4])

TA_shuffled_50 <- sum(mean_n_TA$`50_to_60` >= m$Disease_Ochiai_TA[1])/20000 
TA_shuffled_60 <- sum(mean_n_TA$`60_to_70` >= m$Disease_Ochiai_TA[2])/20000
TA_shuffled_70 <- sum(mean_n_TA$`70_to_80` >= m$Disease_Ochiai_TA[3])/20000
TA_shuffled_80 <- sum(mean_n_TA$`80_to_90` >= m$Disease_Ochiai_TA[4])/20000
p_TA_Post = c(TA_shuffled_50, TA_shuffled_60, TA_shuffled_70, TA_shuffled_80)
p_TA_Post = p.adjust(p_TA_Post, method = "fdr", n = length(p_TA_Post))
paste("The p values for the TA 50-59 year subnetworks based on posterior score are: ", p_TA_Post[1])
paste("The p values for the TA 60-69 year subnetworks based on posterior score are: ", p_TA_Post[2])
paste("The p values for the TA 70-79 year subnetworks based on posterior score are: ", p_TA_Post[3])
paste("The p values for the TA 80-89 year subnetworks based on posterior score are: ", p_TA_Post[4])

EA_shuffled_50 <- sum(mean_n_EA$`50_to_60` >= m$Disease_Ochiai_EA[1])/20000
EA_shuffled_60 <- sum(mean_n_EA$`60_to_70` >= m$Disease_Ochiai_EA[2])/20000
EA_shuffled_70 <- sum(mean_n_EA$`70_to_80` >= m$Disease_Ochiai_EA[3])/20000 
EA_shuffled_80 <- sum(mean_n_EA$`80_to_90` >= m$Disease_Ochiai_EA[4])/20000
p_EA_Post = c(EA_shuffled_50, EA_shuffled_60, EA_shuffled_70, EA_shuffled_80)
p_EA_Post = p.adjust(p_EA_Post, method = "fdr", n = length(p_EA_Post))
paste("The p values for the EA 50-59 year subnetworks based on posterior score are: ", p_EA_Post[1])
paste("The p values for the EA 60-69 year subnetworks based on posterior score are: ", p_EA_Post[2])
paste("The p values for the EA 70-79 year subnetworks based on posterior score are: ", p_EA_Post[3])
paste("The p values for the EA 80-89 year subnetworks based on posterior score are: ", p_EA_Post[4])

LOP_shuffled_50 <- sum(mean_n_LOP$`50_to_60` >= m$Disease_Ochiai_LOP[1])/20000
LOP_shuffled_60 <- sum(mean_n_LOP$`60_to_70` >= m$Disease_Ochiai_LOP[2])/20000
LOP_shuffled_70 <- sum(mean_n_LOP$`70_to_80` >= m$Disease_Ochiai_LOP[3])/20000
LOP_shuffled_80 <- sum(mean_n_LOP$`80_to_90` >= m$Disease_Ochiai_LOP[4])/20000
p_LOP_Post = c(LOP_shuffled_50, LOP_shuffled_60, LOP_shuffled_70, LOP_shuffled_80)
p_LOP_Post = p.adjust(p_LOP_Post, method = "fdr", n = length(p_LOP_Post))
paste("The p values for the LOP 50-59 year subnetworks based on posterior score are: ", p_LOP_Post[1])
paste("The p values for the LOP 60-69 year subnetworks based on posterior score are: ", p_LOP_Post[2])
paste("The p values for the LOP 70-79 year subnetworks based on posterior score are: ", p_LOP_Post[3])
paste("The p values for the LOP 80-89 year subnetworks based on posterior score are: ", p_LOP_Post[4])

DNS_shuffled_50 <- sum(mean_n_DNS$`50_to_60` >= m$Disease_Ochiai_DNS[1])/20000
DNS_shuffled_60 <- sum(mean_n_DNS$`60_to_70` >= m$Disease_Ochiai_DNS[2])/20000
DNS_shuffled_70 <- sum(mean_n_DNS$`70_to_80` >= m$Disease_Ochiai_DNS[3])/20000
DNS_shuffled_80 <- sum(mean_n_DNS$`80_to_90` >= m$Disease_Ochiai_DNS[4])/20000
p_DNS_Post = c(DNS_shuffled_50, DNS_shuffled_60, DNS_shuffled_70, DNS_shuffled_80)
p_DNS_Post = p.adjust(p_DNS_Post, method = "fdr", n = length(p_DNS_Post))
paste("The p values for the DNS 50-59 year subnetworks based on posterior score are: ", p_DNS_Post[1])
paste("The p values for the DNS 60-69 year subnetworks based on posterior score are: ", p_DNS_Post[2])
paste("The p values for the DNS 70-79 year subnetworks based on posterior score are: ", p_DNS_Post[3])
paste("The p values for the DNS 80-89 year subnetworks based on posterior score are: ", p_DNS_Post[4])

MD_shuffled_50 <- sum(mean_n_MD$`50_to_60` >= m$Disease_Ochiai_MD[1])/20000
MD_shuffled_60 <- sum(mean_n_MD$`60_to_70` >= m$Disease_Ochiai_MD[2])/20000 
MD_shuffled_70 <- sum(mean_n_MD$`70_to_80` >= m$Disease_Ochiai_MD[3])/20000 
MD_shuffled_80 <-sum(mean_n_MD$`80_to_90` >= m$Disease_Ochiai_MD[4])/20000 
p_MD_Post = c(MD_shuffled_50, MD_shuffled_60, MD_shuffled_70, MD_shuffled_80)
p_MD_Post = p.adjust(p_MD_Post, method = "fdr", n = length(p_MD_Post))
paste("The p values for the MD 50-59 year subnetworks based on posterior score are: ", p_MD_Post[1])
paste("The p values for the MD 60-69 year subnetworks based on posterior score are: ", p_MD_Post[2])
paste("The p values for the MD 70-79 year subnetworks based on posterior score are: ", p_MD_Post[3])
paste("The p values for the MD 80-89 year subnetworks based on posterior score are: ", p_MD_Post[4])

CS_shuffled_50 <- sum(mean_n_CS$`50_to_60` >= m$Disease_Ochiai_CS[1])/20000
CS_shuffled_60 <- sum(mean_n_CS$`60_to_70` >= m$Disease_Ochiai_CS[2])/20000 
CS_shuffled_70 <- sum(mean_n_CS$`70_to_80` >= m$Disease_Ochiai_CS[3])/20000 
CS_shuffled_80 <- sum(mean_n_CS$`80_to_90` >= m$Disease_Ochiai_CS[4])/20000 
p_CS_Post = c(CS_shuffled_50, CS_shuffled_60, CS_shuffled_70, CS_shuffled_80)
p_CS_Post = p.adjust(p_CS_Post, method = "fdr", n = length(p_CS_Post))
paste("The p values for the CS 50-59 year subnetworks based on posterior score are: ", p_CS_Post[1])
paste("The p values for the CS 60-69 year subnetworks based on posterior score are: ", p_CS_Post[2])
paste("The p values for the CS 70-79 year subnetworks based on posterior score are: ", p_CS_Post[3])
paste("The p values for the CS 80-89 year subnetworks based on posterior score are: ", p_CS_Post[4])

SCE_shuffled_50 <- sum(mean_n_SCE$`50_to_60` >= m$Disease_Ochiai_SCE[1])/20000
SCE_shuffled_60 <- sum(mean_n_SCE$`60_to_70` >= m$Disease_Ochiai_SCE[2])/20000
SCE_shuffled_70 <- sum(mean_n_SCE$`70_to_80` >= m$Disease_Ochiai_SCE[3])/20000
SCE_shuffled_80 <- sum(mean_n_SCE$`80_to_90` >= m$Disease_Ochiai_SCE[4])/20000
p_SCE_Post = c(SCE_shuffled_50,SCE_shuffled_60,SCE_shuffled_70,SCE_shuffled_80)
p_SCE_Post = p.adjust(p_SCE_Post, method = "fdr", n = length(p_SCE_Post))
paste("The p values for the SCE 50-59 year subnetworks based on posterior score are: ", p_SCE_Post[1])
paste("The p values for the SCE 60-69 year subnetworks based on posterior score are: ", p_SCE_Post[2])
paste("The p values for the SCE 70-79 year subnetworks based on posterior score are: ", p_SCE_Post[3])
paste("The p values for the SCE 80-89 year subnetworks based on posterior score are: ", p_SCE_Post[4])

AIC_shuffled_50 <- sum(mean_n_AIC$`50_to_60` >= m$Disease_Ochiai_AIC[1])/20000 
AIC_shuffled_60 <- sum(mean_n_AIC$`60_to_70` >= m$Disease_Ochiai_AIC[2])/20000
AIC_shuffled_70 <- sum(mean_n_AIC$`70_to_80` >= m$Disease_Ochiai_AIC[3])/20000
AIC_shuffled_80 <- sum(mean_n_AIC$`80_to_90` >= m$Disease_Ochiai_AIC[4])/20000
p_AIC_Post = c(AIC_shuffled_50, AIC_shuffled_60, AIC_shuffled_70, AIC_shuffled_80)
p_AIC_Post = p.adjust(p_AIC_Post, method = "fdr", n = length(p_AIC_Post))
paste("The p values for the AIC 50-59 year subnetworks based on posterior score are: ", p_AIC_Post[1])
paste("The p values for the AIC 60-69 year subnetworks based on posterior score are: ", p_AIC_Post[2])
paste("The p values for the AIC 70-79 year subnetworks based on posterior score are: ", p_AIC_Post[3])
paste("The p values for the AIC 80-89 year subnetworks based on posterior score are: ", p_AIC_Post[4])

m <- round(m, digits = 4)
m <- t(m)
m_p_vals_Post <- data.frame(p_GI_Post, p_TA_Post, p_EA_Post, p_LOP_Post, p_DNS_Post, p_MD_Post, p_CS_Post, p_SCE_Post, p_AIC_Post)
m_p_vals_Post <- round(m_p_vals_Post, digits = 4)
m_p_vals_Post <- t(m_p_vals_Post)


# Plotting subgraphs based on the posterior score for 60 - 69 year age category
par(mfrow =c(1,1), mar=c(1,1,1,1), oma = c(0,0,0,0))

# Node names
ochiai_frame_ARDs3 <- ochiai_frame_ARDs3[,c(1,2)]
colnames(a_60)[1] <- c("name_node")
a_60_3 <- merge(a_60, ochiai_frame_ARDs3, by = "name_node")
a_60_3 <- a_60_3[, c(1, 11, 2:10)]

# EA - posterior score
mean_dist_EA_60 <- order_frame2(a_60, g_60, a_60_3, 4, "EA", "Age 60 to 69:", 30, layout2=layout_with_fr)
title(paste("EA: 60 to 69 years, p<0.01"))

# DNS - posterior score
mean_dist_DNS_60 <- order_frame2(a_60, g_60, a_60_3, 6, "DNS", "Age 60 to 69:", 30, layout2=layout_with_graphopt)
title(paste("DNS: 60 to 69 years, p<0.01"))

# MD - posterior score
mean_dist_MD_60 <- order_frame2(a_60, g_60, a_60_3, 7, "MD", "Age 60 to 69:", 30, layout2=layout_with_fr)
title(paste("MD: 60 to 69 years, p<0.05"))

# SCE - posterior score
mean_dist_SCE_60 <- order_frame2(a_60, g_60, a_60_3, 8, "SCE", "Age 60 to 69:", 30, layout2=layout_with_graphopt)
title(paste("SCE: 60 to 69 years, p<0.05"))

# AIC - posterior score
mean_dist_AIC_60 <- order_frame2(a_60, g_60, a_60_3, 10, "AIC", "Age 60 to 69:", 30, layout2=layout_with_graphopt)
title(paste("AIC: 60 to 69 years, p<0.05"))

# Look at newly prioritzed ARDs for each significant AH --> EA, DNS, MD, SCE, AIC
colnames(a_50) <- c("name_node", "GI_Post", "TA_Post", "EA_Post", "LOP_Post", "DNS_Post", "MD_Post", "SCE_Post","CS_Post", "AIC_Post")
new_frame_50 <- merge(ochiai_frame_ARDs, a_50, by="name_node", all.x =TRUE, all.y =TRUE)
new_frame_50 <- new_frame_50[, c("name_node", "Disease_Ochiai_EA", "EA_Post", "Disease_Ochiai_DNS", "DNS_Post",
                                 "Disease_Ochiai_MD", "MD_Post", "Disease_Ochiai_SCE", "SCE_Post", "Disease_Ochiai_AIC", "AIC_Post")]
write.csv(new_frame_50, "shuffled_ochiai/new_frame_50.csv")
colnames(a_60) <- c("name_node", "GI_Post", "TA_Post", "EA_Post", "LOP_Post", "DNS_Post", "MD_Post", "SCE_Post","CS_Post", "AIC_Post")
new_frame_60 <- merge(ochiai_frame_ARDs, a_60, by="name_node", all.x =TRUE, all.y =TRUE)
new_frame_60 <- new_frame_60[, c("name_node", "Disease_Ochiai_EA", "EA_Post", "Disease_Ochiai_DNS", "DNS_Post",
                                 "Disease_Ochiai_MD", "MD_Post", "Disease_Ochiai_SCE", "SCE_Post", "Disease_Ochiai_AIC", "AIC_Post")]
write.csv(new_frame_60, "shuffled_ochiai/new_frame_60.csv")
colnames(a_70) <- c("name_node", "GI_Post", "TA_Post", "EA_Post", "LOP_Post", "DNS_Post", "MD_Post", "SCE_Post","CS_Post", "AIC_Post")
new_frame_70 <- merge(ochiai_frame_ARDs, a_70, by="name_node", all.x =TRUE, all.y =TRUE)
new_frame_70 <- new_frame_70[, c("name_node", "Disease_Ochiai_EA", "EA_Post", "Disease_Ochiai_DNS", "DNS_Post",
                                 "Disease_Ochiai_MD", "MD_Post", "Disease_Ochiai_SCE", "SCE_Post", "Disease_Ochiai_AIC", "AIC_Post")]
write.csv(new_frame_70, "shuffled_ochiai/new_frame_70.csv")
colnames(a_80) <- c("name_node", "GI_Post", "TA_Post", "EA_Post", "LOP_Post", "DNS_Post", "MD_Post", "SCE_Post","CS_Post", "AIC_Post")
new_frame_80 <- merge(ochiai_frame_ARDs, a_80, by="name_node", all.x =TRUE, all.y =TRUE)
new_frame_80 <- new_frame_80[, c("name_node", "Disease_Ochiai_EA", "EA_Post", "Disease_Ochiai_DNS", "DNS_Post",
                                 "Disease_Ochiai_MD", "MD_Post", "Disease_Ochiai_SCE", "SCE_Post", "Disease_Ochiai_AIC", "AIC_Post")]
write.csv(new_frame_80, "shuffled_ochiai/new_frame_80.csv")

# Mapping newly prioritized ARDs to associated genes
ARD_Gene_frame_lower_all <- read.csv("GWAS_catalog/ARD_Gene_frame_all_lower_NCBI_final_50000.csv")
essential_tremor <- ARD_Gene_frame_lower_all[grep("essential_tremor", ARD_Gene_frame_lower_all$name_node),]
bells <- ARD_Gene_frame_lower_all[grep("bells", ARD_Gene_frame_lower_all$name_node),]
