################################################################################################
#Title: Trajectory Analysis 
#Author: Meltem Omur 
################################################################################################
library(Monocle)

esc.data.scenic.seurat <- readRDS("../clusters_assigned.rds") #output of Seurat analysis

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(esc.data.scenic.seurat@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = esc.data.scenic.seurat@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
esc.monocle <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


#run this on the server 
esc.monocle <- estimateSizeFactors(esc.monocle)
esc.monocle <- estimateDispersions(esc.monocle, cores=4)

esc.monocle <- detectGenes(esc.monocle, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(esc.monocle), num_cells_expressed >= 10))

# add a UMI column into phenoData
# use Matrix because data is in the sparse format
pData(esc.monocle)$UMI <- Matrix::colSums(exprs(esc.monocle))

pdf("umigeompont.pdf")
ggplot(pData(esc.monocle), aes(num_genes_expressed, UMI)) + geom_point()
dev.off()

#distribution over UMI counts (without UNKs)
y <- pData(esc.monocle) %>% filter(clusters %in% c("DE", "hESC", "ME", "PS"))
x <- y$UMI
df <- data.frame(x = x)
pdf("histogram.pdf")
ggplot(df, aes(x)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(4000, 20000), linetype = "dotted", color = 'red')
  dev.off()

fData(esc.monocle)$use_for_ordering <- fData(esc.monocle)$num_cells_expressed > 0.05 * ncol(esc.monocle)

# how many genes are used?
table(fData(esc.monocle)$use_for_ordering)

## Constructing Single Cell Trajectories 

esc.monocle <- reduceDimension(esc.monocle, max_components = 2, norm_method = 'log',num_dim = 50, reduction_method = 'tSNE', verbose = T)

esc_expressed_genes <-  row.names(subset(fData(esc.monocle),num_cells_expressed >= 10))

clustering_DEG_genes <-differentialGeneTest(esc.monocle[esc_expressed_genes,],fullModelFormulaStr = '~clusters',cores = 4)

esc_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

esc.monocle <- setOrderingFilter(esc.monocle, ordering_genes = esc_ordering_genes)

esc.monocle <- reduceDimension(esc.monocle, method = 'DDRTree', cores=4)

esc.monocle <- orderCells(esc.monocle, reverse=F)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$clusters)[,"hESC"]
    return(as.numeric(names(T0_counts)[which
          (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
esc.rooted <- orderCells(esc.monocle, root_state = GM_state(esc.monocle))

pdf("pseudotime_cluster_starting_from_beginning.pdf")
plot_cell_trajectory(esc.rooted, color_by = "clusters")
dev.off()

pdf("pseudotime_cluster_color_by_pseudotime_starting_from_beginning.pdf")
plot_cell_trajectory(esc.rooted, color_by = "Pseudotime")
dev.off()

#saveRDS(esc.rooted, "esc_rooted_monocle_object_paper.rds")
monocle_plotting <- readRDS("esc_rooted_monocle_object_paper")


######################### Pseudotime clustering with spesific colors ######################################################

pse_cluster <- plot_cell_trajectory(mono, color_by = "clusters", cell_size = 0.5)+scale_y_reverse() +
  scale_color_manual(breaks = c("DE","hESC","ME","PS"), values= c("#dc050c","#1965b0","#f7f056","#4eb265"))
ggsave(file="pseudotime_plot_by_cluster_colored_by_tsne_colors.svg", plot=pse_cluster, width=4, height=4)

pse_pseu<- plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime", cell_size = 0.5)+scale_y_reverse()+scale_color_gradient(low="blue", high="red")
ggsave(file="pseudotime_plot_by_pseudotime.svg", plot=pse_pseu, width=4, height=4)

########################## Marker expressions throughout trajectories #######################################################
mesoderm_markers<- row.names(subset(fData(monocle_plotting),gene_short_name %in% c("TBX6", "TBXT", "MSX1","FOXH1")))
meso_mono_plot <- plot_genes_branched_pseudotime(monocle_plotting[mesoderm_markers,], 
                                       branch_point = 1,
                                       color_by = "clusters",
                                       ncol = 1, cell_size=0.4, panel_order = c("FOXH1","MSX1","TBX6","TBXT"))
meso_1 <- meso_mono_plot +theme(axis.text.x =  element_text(size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10),
                      legend.text = element_text(size = 8), legend.title = element_text(size = 9))
ggsave(file="mesoderm_markers_pseudotime_4x4.svg", plot=meso_1, width=4, height=4)



	
