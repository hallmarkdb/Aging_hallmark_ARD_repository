library(stringr)
library(eply)
library(plyr)
library(tidyr)
library(eulerr)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(dplyr)
setwd("/Code_Thesis_HCF/Aging_hallmark_ARD_repository/")

#################################### Heatmaps of Ochiai ############################################# 
process_frame <- function(dataframe_new, cut_off, data_AH) {
  # process the data frame
  dataframe_new <- data.frame(dataframe_new$Disease, dataframe_new$Association)
  dataframe_new <- drop_na(dataframe_new)
  dataframe_new <- ddply(dataframe_new, .(dataframe_new.Disease), summarize, Xc=sum(dataframe_new.Association))
  dataframe_new$Xc <- ifelse(dataframe_new$Xc >cut_off, 1, 0)
  colnames(data_AH) <- c("x","ochiai")
  colnames(dataframe_new) <- c("x", "true")
  new_dataframe_AH <- merge(dataframe_new, data_AH, by="x", all.x = TRUE, all.y =TRUE)
  new_dataframe_AH <- drop_na(new_dataframe_AH)
  new_dataframe_AH$multiple <- (new_dataframe_AH$true*new_dataframe_AH$ochiai)
  return(new_dataframe_AH)
}

# Import the aging hallmark-ARD association scores
data <- read.csv("data/aging_hallmark_subnetworks/counts_frame_0321.csv") 
data <- subset(data, data$Disease_Absolute_Count > 250)
write.csv(data, "data/aging_hallmark_subnetworks/counts_frame_updated_0321.csv")
data_updated <- read.csv("data/aging_hallmark_subnetworks/counts_frame_0321.csv")
colnames(data_updated)[1] <- c("Disease")
data_updated$Disease <- trimws(data_updated$Disease)
Disease_list <- data.frame(data_updated$Disease)
colnames(Disease_list)[1] <- c("Disease")
data_updated_copy <- data_updated

# Import the benchmark sentences
benchmark_sentences_GI <- read.csv("benchmark_sentences/benchmark_sentences_GI.csv")
benchmark_sentences_TA <- read.csv("benchmark_sentences/benchmark_sentences_TA.csv")
benchmark_sentences_EA <- read.csv("benchmark_sentences/benchmark_sentences_EA.csv")
benchmark_sentences_LOP <- read.csv("benchmark_sentences/benchmark_sentences_LOP.csv")
benchmark_sentences_DNS <- read.csv("benchmark_sentences/benchmark_sentences_DNS.csv")
benchmark_sentences_MD <- read.csv("benchmark_sentences/benchmark_sentences_MD.csv")
benchmark_sentences_CS <- read.csv("benchmark_sentences/benchmark_sentences_CS.csv")
benchmark_sentences_SCE <- read.csv("benchmark_sentences/benchmark_sentences_SCE.csv")
benchmark_sentences_AIC <- read.csv("benchmark_sentences/benchmark_sentences_AIC.csv")
benchmark_sentences <- unique(rbind(benchmark_sentences_GI, benchmark_sentences_TA, benchmark_sentences_EA, 
                                    benchmark_sentences_LOP, benchmark_sentences_DNS, benchmark_sentences_MD, 
                                    benchmark_sentences_CS, benchmark_sentences_SCE, benchmark_sentences_AIC))

# Less than 2,500 co-mentioning sentences - GI, TA, LOP, MD, CS, SCE --> one sentence to verify a co-occurrence score
data_GI <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_GI)
new_dataframe_GI <- process_frame(benchmark_sentences_GI, 0, data_GI)
data_TA <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_TA)
new_dataframe_TA <- process_frame(benchmark_sentences_TA, 0, data_TA)
data_LOP <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_LOP)
new_dataframe_LOP <- process_frame(benchmark_sentences_LOP, 0, data_LOP)
data_MD <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_MD)
new_dataframe_MD <- process_frame(benchmark_sentences_MD, 0, data_MD)
data_CS <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_CS)
new_dataframe_CS <- process_frame(benchmark_sentences_CS, 0, data_CS)
data_SCE <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_SCE)
new_dataframe_SCE <- process_frame(benchmark_sentences_SCE, 0, data_SCE)

# More than 2,500 co-mentioning sentences - DNS, EA, AIC --> three sentences to verify a co-occurrence score
data_DNS <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_DNS)
new_dataframe_DNS <- process_frame(benchmark_sentences_DNS, 2, data_DNS)
data_EA <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_EA)
new_dataframe_EA <- process_frame(benchmark_sentences_EA, 2, data_EA)
data_AIC <- data.frame(data_updated$Disease, data_updated$Disease_Ochiai_AIC)
new_dataframe_AIC <- process_frame(benchmark_sentences_AIC, 2, data_AIC)

# Updating the data frame
data_updated_copy <- data_updated_copy[,c(1:10)]
colnames(data_updated_copy)[1] <- c("X")
colnames(new_dataframe_GI) <- c("X", "GI_true", "GI_ochiai", "GI_multiple")
colnames(new_dataframe_TA) <- c("X", "TA_true", "TA_ochiai", "TA_multiple")
colnames(new_dataframe_EA) <- c("X", "EA_true", "EA_ochiai", "EA_multiple")
colnames(new_dataframe_LOP) <- c("X", "LOP_true", "LOP_ochiai", "LOP_multiple")
colnames(new_dataframe_DNS) <- c("X", "DNS_true", "DNS_ochiai", "DNS_multiple")
colnames(new_dataframe_MD) <- c("X", "MD_true", "MD_ochiai", "MD_multiple")
colnames(new_dataframe_CS) <- c("X", "CS_true", "CS_ochiai", "CS_multiple")
colnames(new_dataframe_SCE) <- c("X", "SCE_true", "SCE_ochiai", "SCE_multiple")
colnames(new_dataframe_AIC) <- c("X", "AIC_true", "AIC_ochiai", "AIC_multiple")
new_dataframe <- merge(data, new_dataframe_GI, by = "X", all.x = TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_TA, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_EA, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_LOP, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_DNS, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_MD, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_SCE, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_CS, by = "X", all.x=TRUE)
new_dataframe <- merge(new_dataframe, new_dataframe_AIC, by = "X", all.x=TRUE)
#new_dataframe <- merge(new_dataframe, data_updated_copy, by ="X", all.x=TRUE)
#new_dataframe <- new_dataframe[,c(1,48,49,2:47)]
new_dataframe$GI_multiple <- replace_na(new_dataframe$GI_multiple, 0)
new_dataframe$TA_multiple <- replace_na(new_dataframe$TA_multiple, 0)
new_dataframe$EA_multiple <- replace_na(new_dataframe$EA_multiple, 0)
new_dataframe$LOP_multiple <- replace_na(new_dataframe$LOP_multiple, 0)
new_dataframe$DNS_multiple <- replace_na(new_dataframe$DNS_multiple, 0)
new_dataframe$MD_multiple <- replace_na(new_dataframe$MD_multiple, 0)
new_dataframe$CS_multiple <- replace_na(new_dataframe$CS_multiple, 0)
new_dataframe$SCE_multiple <- replace_na(new_dataframe$SCE_multiple, 0)
new_dataframe$AIC_multiple <- replace_na(new_dataframe$AIC_multiple, 0)
write.csv(new_dataframe, "data/aging_hallmark_subnetworks/ochiai_frame_0321.csv")

# Data frame for heatmap of top 30 ranked ARDs per aging hallmark
GI <- new_dataframe[order(-new_dataframe[,"GI_multiple"]),]
GI <- as.vector(GI[1:30,c(1,23)])
GI$GI_order <- c(1:30)
colnames(GI) <- c("Name", "GI_multiple", "GI_order")
TA <- new_dataframe[order(-new_dataframe[,"TA_multiple"]),]
TA <- as.vector(TA[1:30,c(1,25)])
TA$TA_order <- c(1:30)
colnames(TA) <- c("Name", "TA_multiple", "TA_order")
EA <- new_dataframe[order(-new_dataframe[,"EA_multiple"]),]
EA <- as.vector(EA[1:30,c(1,29)])
EA$EA_order <- c(1:30)
colnames(EA) <- c("Name", "EA_multiple", "EA_order")
LOP <- new_dataframe[order(-new_dataframe[,"LOP_multiple"]),]
LOP <- as.vector(LOP[1:30,c(1,32)])
LOP$LOP_order <- c(1:30)
colnames(LOP) <- c("Name", "LOP_multiple", "LOP_order")
DNS <- new_dataframe[order(-new_dataframe[,"DNS_multiple"]),]
DNS <- as.vector(DNS[1:30,c(1,35)])
DNS$DNS_order <- c(1:30)
colnames(DNS) <- c("Name", "DNS_multiple", "DNS_order")
MD <- new_dataframe[order(-new_dataframe[,"MD_multiple"]),]
MD <- as.vector(MD[1:30,c(1,38)])
MD$MD_order <- c(1:30)
colnames(MD) <- c("Name", "MD_multiple", "MD_order")
SCE <- new_dataframe[order(-new_dataframe[,"SCE_multiple"]),]
SCE <- as.vector(SCE[1:30,c(1,41)])
SCE$SCE_order <- c(1:30)
colnames(SCE) <- c("Name", "SCE_multiple", "SCE_order")
CS <- new_dataframe[order(-new_dataframe[,"CS_multiple"]),]
CS <- as.vector(CS[1:30,c(1,44)])
CS$CS_order <- c(1:30)
colnames(CS) <- c("Name", "CS_multiple", "CS_order")
AIC <- new_dataframe[order(-new_dataframe[,"AIC_multiple"]),]
AIC <- as.vector(AIC[1:30,c(1,47)])
AIC$AIC_order <- c(1:30)
colnames(AIC) <- c("Name", "AIC_multiple", "AIC_order")

frames_merge <- merge(GI, TA, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, EA, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, LOP, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, DNS, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, MD, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, CS, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, SCE, by="Name", all.x = TRUE, all.y = TRUE)
frames_merge <- merge(frames_merge, AIC, by="Name", all.x = TRUE, all.y = TRUE)
frames_AH <- frames_merge[,c(1,3,5,7,9,11,13,15,17,19)]
frames_AH <- frames_AH[order(frames_AH[,2]),]
frames_AH <- frames_AH[order(frames_AH[,3]),]
frames_AH <- frames_AH[order(frames_AH[,4]),]
frames_AH <- frames_AH[order(frames_AH[,5]),]
frames_AH <- frames_AH[order(frames_AH[,6]),]
frames_AH <- frames_AH[order(frames_AH[,7]),]
frames_AH <- frames_AH[order(frames_AH[,8]),]
frames_AH <- frames_AH[order(frames_AH[,9]),]
frames_AH <- frames_AH[order(frames_AH[,10]),]
frames_AH <- frames_AH[order(frames_AH[,1]),]

heatmap_frame <- as.matrix(frames_AH[,c(2:10)])
rownames(heatmap_frame) <- frames_AH$Name
colnames(heatmap_frame) <- c("GI", "TA", "EA", "LOP", "DNS", "MD", "CS", "SCE", "AIC")
col2=colorRampPalette((brewer.pal(9, "YlOrRd")))(325)
pheatmap(heatmap_frame, na_col = 'white', cluster_rows = FALSE, col=rev(col2), row_names_side = c("left"), cluster_cols = FALSE, color_space = "LAB",
         fontsize_row = 10, fontsize_col = 10, heatmap_legend_param = list(title = "rank", rect_gp = gpar(col = "black", lwd = 0.7),
                                                                           legend_height = unit(3, "cm"), title_position = "left"))

# Heatmap of all ARDs associated with aging hallmarks
all_AH_df <- new_dataframe[,c(1, 23, 27, 29, 32, 35, 38, 44, 41, 46)]
all_AH_df <- all_AH_df[order(all_AH_df[,2], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,3], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,4], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,5], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,6], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,7], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,8], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,9], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,10], decreasing = TRUE),]
all_AH_df <- all_AH_df[order(all_AH_df[,1], decreasing = FALSE),]
all_AH_df[all_AH_df == 0] <- NA
heatmap_frame2 <- as.matrix(all_AH_df[,c(2:10)])
rownames(heatmap_frame2) <- all_AH_df$X
heatmap_frame2 <- heatmap_frame2[rowSums(is.na(heatmap_frame2)) != ncol(heatmap_frame2), ]
col3 = colorRampPalette((brewer.pal(9, "YlOrRd")))(325)
colnames(heatmap_frame2) <- c("GI", "TA", "EA", "LOP", "DNS", "MD", " CS", "SCE", "AIC")
pheatmap(log(heatmap_frame2), na_col = 'white', col = col3, cluster_rows = FALSE, row_names_side = c("left"), cluster_cols = FALSE, color_space = "LAB",
         fontsize_row = 10, fontsize_col = 10, heatmap_legend_param = list(title = "rank", rect_gp = gpar(col = "black", lwd = 0.7),
                                                                           legend_height = unit(2, "cm"), title_position = "left"))

# Adding name node and disease abbreviations to ochiai frame --> aging hallmark profiles
ochiai_frame <- read.csv("data/aging_hallmark_subnetworks/ochiai_frame_0321.csv")
abbreviation <- read.csv("data/aging_hallmark_subnetworks/mappings_ARD_abbreviation.csv")
merge_frame <- merge(ochiai_frame, abbreviation, by="X")
merge_frame$X.1 <- NULL
merge_frame <- merge_frame[,c(1,48,49, 2:47)]
write.csv(merge_frame, "data/aging_hallmark_subnetworks/ochiai_frame_ARDs.csv")
