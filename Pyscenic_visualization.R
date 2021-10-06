#############################################################################################
#Title: PyScenic output visualization
#Author: Meltem Omur 
#############################################################################################
#Setup: 
library(Seurat)
library(SCENIC)
library(loomR)
library(AUCell)
library(RcisTarget)
library(SCopeLoomR)
##plotting
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(ComplexHeatmap) 
library(data.table)
library(ggplot2)

####################################################################################################
#### Open Pyscenic loom output file and retreive regulon information for the visualization 		
#### Sox17Ng_aggr_36h_72hDE_72hM_filtered_AUCell.loom is the output of pyscenic. 
####################################################################################################
loom <- open_loom("Sox17Ng_aggr_36h_72hDE_72hM_filtered_AUCell.loom", mode="r") 
exprMat <- get_dgem(loom)
dim(exprMat)
exprMat_log <- log(exprMat+1)
cellInfo <- get_cell_annotation(loom)
saveRDS(cellInfo, file="cellInfo.Rds")
# Read information from loom file:
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat) 
regulonsAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")## regulonsAUC include regulon information (a TF and its target genes). 
saveRDS(regulonsAUC, "regulonsAUC.rds") # regulonsAUC will be used for the visualization.
loom$close_all()

#############################################################################################################################
#### For further analysis, https://github.com/FloWuenne/scFunctions/blob/master/Tutorials/process_SCENIC.md was followed #### 
#############################################################################################################################
install_github("FloWuenne/scFunctions")
library("scFunctions")
library(tidyverse)

regulonAUC <- readRDS("regulonsAUC.rds") ## created using SCopeLoomR package. 
kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)
head(kmeans_thresholds)

###Binarize regulons using thresholds
binary_regulons <- binarize_regulons(regulonAUC,kmeans_thresholds)

#Let's take a look at the first regulon in the binary regulon list.
head(binary_regulons$`AHCTF1(+)`)

###Next, reformat the binary regulons into a big data frame that contains all of the binary regulons to calculate RRS scores.

joined_bin_reg <- binary_regulons %>%
  reduce(left_join,by="cells")
rownames(joined_bin_reg) <- joined_bin_reg$cells
joined_bin_reg <- joined_bin_reg[2:ncol(joined_bin_reg)]
binary_regulons_trans <- as.matrix(t(joined_bin_reg))

#Check that the data table is formatted correctly before proceeding:
binary_regulons_trans[1:4,1:3]
head(binary_regulons_trans)
hh_annot_col <- data.frame(row.names = rownames(cellInfo), cluster= cellInfo$clusters )
write.csv(hh_annot_col,"cluster_annotation.csv")
write.csv(binary_regulons_trans,"binary_matrix_of_regulons.csv")

pp<-pheatmap(
  binary_regulon_matrix,color = colorRampPalette(c("white", "gray", "black"))(100),
  annotation_col = head(hh_annot_col,50), show_rownames = F, show_colnames = F,
  clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean",clustering_method = "ward.D")
ggsave(file="binary_regulon_heatmap_01.pdf", plot=pp, width=10, height=10)

#######		Calculate Regulon Specificity Score (RSS) 	#####################################################################################################
#We now want to use the binary regulon activity together with the cell assignments to see how specific each predicted regulon is for each cell type. 
#Regulon specificity score (RSS) is calculated based on the Jensen-Shannon divergence, a measure of the similarity between two probability distributions. 
#Basically for the calculation of the RSS, we will calculate the Jensen-Shannon divergence between each vector of binary regulon activity overlaps with the 
#assignment of cells to a specific cell type (For further details, check https://github.com/FloWuenne/scFunctions/blob/master/Tutorials/process_SCENIC.md)
#############################################################################################################################################################
metadata_sub <- readRDS("cellInfo.Rds")
head(metadata_sub)

#Now that we are ready to calculate the RSS for all regulons over all cell types.
rrs_df <- calculate_rrs(metadata_sub,
                        binary_regulons = binary_regulons_trans, cell_type_column = "clusters")

#The output is a data frame with a RSS score for each regulon 
head(rrs_df)

#We can visualize the RSS by performing ranking on the RSS scores with
#the most specific regulons ranking the highest per cell type.
plot_rrs_ranking(rrs_df, cell_type = "DE",
                 ggrepel_force =1,
                 ggrepel_point_padding = 0.1,
                 top_genes = 4,
                 plot_extended = FALSE)

#Let's check the distribution of RSS over all cell types.
library(ggridges)
rrs_df_nona <- subset(rrs_df,RSS > 0)
ggplot(rrs_df_nona,aes(RSS,cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 5, alpha = 0.75) +
  geom_vline(xintercept = 0.1) +
  theme(legend.position = "none")

#For highlighting purposes, we are gonna filter with 0.35 to be able to easily visualize the result in a heatmap.
rrs_df_wide <- rrs_df %>%
  spread(cell_type,RSS)
rownames(rrs_df_wide) <- rrs_df_wide$regulon 
rrs_df_wide <- rrs_df_wide[,2:ncol(rrs_df_wide)]

### Subset regulons based on RSS score
rrs_df_wide_specific <- rrs_df_wide[apply(rrs_df_wide,MARGIN = 1 ,FUN =  function(x) any(x > 0.35)),]

###Visualize the regulons in heatmap
top10regulonsmanuels<- read.csv("top10_rrs_df_nano_scores_per_cluster.txt", header = F) #manually created file which contains top 10 regulons for each cluster. 
subs<- subset((rrs_df_wide),rownames(rrs_df_wide) %in% top10regulonsmanuels$V1)
##### Create the heatmap of top 10 regulons for each cell cluster 
ComplexHeatmap::Heatmap(subs,column_order = c("hESC","PS","DE","ME"), 
                        col = colorRampPalette(brewer.pal(9,"Reds"))(25), show_row_names = T,
                        column_names_gp = grid::gpar(fontsize=8),row_names_gp =grid::gpar(fontsize=6),
                        heatmap_width = unit(8, "cm"), heatmap_height = unit(8, "cm"),name="RSS")
