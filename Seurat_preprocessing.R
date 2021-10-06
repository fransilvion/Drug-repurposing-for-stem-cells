################################################################################################
#Title: Seurat Analysis
#Author: Meltem Omur 
################################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)

esc.data <- Read10X(data.dir = ".../filtered_feature_bc_matrix/")
esc.data <- CreateSeuratObject(counts = esc.data, project = "esc", min.cells = 3, min.features = 200)
#data<- esc.data

esc.data[["percent.mt"]] <- PercentageFeatureSet(esc.data, pattern = "^MT-")

### Visualize QC metrics as a violin plot
VlnPlot(esc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

###
plot1 <- FeatureScatter(esc.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(esc.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```
###
esc.data <- subset(esc.data, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)


#Normalization. Normalized values are stored in data[["RNA"]]@data.
esc.data <- NormalizeData(esc.data, normalization.method = "LogNormalize", scale.factor = 10000)

#Find highly variable genes
esc.data <- FindVariableFeatures(esc.data, selection.method = "vst", nfeatures = 5000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(esc.data), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(esc.data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = F)
plot1 + plot2

#Scaling tha data
all.genes <- rownames(esc.data)
esc.data <- ScaleData(esc.data, features = all.genes)

##run PCA
esc.data <- RunPCA(esc.data, features = VariableFeatures(object = esc.data), verbose = FALSE)
DimPlot(esc.data, reduction = "pca")


esc.data <- FindNeighbors(esc.data, dims = 1:50)
esc.data <- FindClusters(esc.data, resolution = 0.6)

```{r}
esc.data <- RunUMAP(esc.data, dims = 1:50)
DimPlot(esc.data, reduction = "umap", label=T)

esc.data <- RunTSNE(esc.data, dims = 1:50)
DimPlot(esc.data, reduction="tsne", label=T)

VlnPlot(esc.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, group.by = "seurat_clusters")

esc.data.markers <- FindAllMarkers(esc.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(esc.data.markers,"all_markers_9clusters.csv")
x <- esc.data.markers %>% group_by(cluster)


#cluster 3  - stem cells
DotPlot(esc.data, features = c("FGF2", "CDH1", "DPPA4", "SOX2", "NANOG", "POU5F1")) + RotatedAxis()
FeaturePlot(esc.data,features = c("FGF2", "CDH1", "DPPA4", "SOX2", "NANOG", "POU5F1"), label = T)
VlnPlot(esc.data,features = c("FGF2", "CDH1", "DPPA4", "SOX2", "NANOG", "POU5F1"),pt.size = 0)

#cluster  9 - primitive streak --
DotPlot(esc.data, features = c("TBXT", "DKK4", "FGF4", "EOMES", "WNT3", "MIXL1")) + RotatedAxis()
FeaturePlot(esc.data,features = c("TBXT", "DKK4", "FGF4", "EOMES", "WNT3", "MIXL1"), label = T)
VlnPlot(esc.data,features = c("TBXT", "DKK4", "FGF4", "EOMES", "WNT3", "MIXL1"), pt.size = 0)

cluster9PS_vshESC_cls2 <- FindMarkers(esc.data, ident.1 = c(9), ident.2 = c(3), min.pct = 0.25)

#cluster 4-6 definitive endoderm
DotPlot(esc.data, features = c("EOMES", "SOX17", "GATA6", "GSC", "CXCR4", "FOXA2")) + RotatedAxis()
FeaturePlot(esc.data,features = c("EOMES", "SOX17", "GATA6", "GSC", "CXCR4", "FOXA2", "CER1"), label = T)
VlnPlot(esc.data,features = c("EOMES", "SOX17", "GATA6", "GSC", "CXCR4", "FOXA2", "CER1"), pt.size = 0)

#cluster 0-2- mesoderm
DotPlot(esc.data, features = c("TBXT", "DKK4", "FGF4", "WNT3", "MIXL1", "BMP4")) + RotatedAxis()
FeaturePlot(esc.data, features = c("TBXT", "DKK4", "FGF4", "WNT3", "MIXL1", "BMP4"), label = T)
VlnPlot(esc.data,  features = c("TBXT", "DKK4", "FGF4", "WNT3", "MIXL1", "BMP4"),pt.size = 0)

clusters <- esc.data@meta.data$seurat_clusters
clusters <- as.character(clusters)
clusters[clusters == "9"] <- "PS"
clusters[clusters == "3"] <- "hESC"
clusters[clusters == "4" | clusters == "6"] <- "DE"
clusters[clusters == "0" | clusters == "2"] <- "ME"
clusters[clusters == "1" | clusters == "5" | clusters == "7" | clusters == "8"] <- "UNK"
clusters <- as.factor(clusters)
esc.data@meta.data$clusters <- clusters

DimPlot(esc.data, reduction = "umap", group.by = "clusters")
```

## Filter out bad quality cells

esc.data2 <- subset(esc.data, subset = nCount_RNA > 2500)
```
```{r}
DimPlot(esc.data2, reduction = "umap", group.by = "clusters")
```
```{r}
esc.data2 <- RunPCA(esc.data2, features = VariableFeatures(object = esc.data2), verbose = FALSE)
esc.data2 <- FindNeighbors(esc.data2, dims = 1:50)
esc.data2 <- FindClusters(esc.data2, resolution = 0.6)
```
```{r}
esc.data2 <- RunUMAP(esc.data2, dims = 1:50)
```
```{r}
DimPlot(esc.data2, reduction = "umap", label=T)
```
esc.data.markers2 <- FindAllMarkers(esc.data2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(esc.data.markers2, "all_markers_after_filtering.csv")

#cluster  2-4- stem cells, 
DotPlot(esc.data2, features = c("FGF2", "CDH1", "DPPA4", "SOX2", "NANOG", "POU5F1")) + RotatedAxis()
VlnPlot(esc.data2, features = c("FGF2", "CDH1", "DPPA4", "SOX2", "NANOG", "POU5F1"), pt.size = 0)

#cluster 9- primitive streak
DotPlot(esc.data2, features = c("TBXT", "DKK4", "FGF4", "EOMES", "WNT3", "MIXL1")) + RotatedAxis()
VlnPlot(esc.data2, features = c("TBXT", "DKK4", "FGF4", "EOMES", "WNT3", "MIXL1"), pt.size = 0)

#cluster  3-5-7 definitive endoderm
DotPlot(esc.data2, features = c("EOMES", "SOX17", "GATA6", "GSC", "FOXA2", "mNeonGreen")) + RotatedAxis()
VlnPlot(esc.data2, features = c("EOMES", "SOX17", "GATA6", "GSC", "FOXA2", "mNeonGreen"), pt.size = 0)

#cluster  0-1 - mesoderm
DotPlot(esc.data2, features = c("TBXT", "MSX1", "MIXL1", "TBX6")) + RotatedAxis()
VlnPlot(esc.data2, features = c("TBXT", "MSX1", "MIXL1", "TBX6"), pt.size = 0, ncol = 2)

clusters <- esc.data2@meta.data$seurat_clusters
clusters <- as.character(clusters)
clusters[clusters == "9"] <- "PS"
clusters[clusters == "2"| clusters == "4"] <- "hESC"
clusters[clusters == "3" | clusters == "5" | clusters == "7"] <- "DE"
clusters[clusters == "0" | clusters == "1" ] <- "ME"
clusters[clusters == "6" | clusters == "8"] <- "UNK"
clusters <- as.factor(clusters)
esc.data2@meta.data$clusters <- clusters

DimPlot(esc.data2, reduction = "umap", group.by = "clusters", label = T)
DimPlot(esc.data2, reduction = "tsne", group.by = "clusters", label = T)


#saveRDS(esc.data2, "esc.data2.filtered_clustered.rds")

tflist <- read.table(file = '../../dbTF.tsv', sep = '\t', header = FALSE)
######################################################################################
###    FIND TF MARKERS BETWEEN hESC AND OTHER CELL TYPES 							 
######################################################################################
esc.markers <- FindMarkers(esc.data2, ident.1 = "hESC", min.pct = 0.20, 
                           group.by = "clusters", logfc.threshold = 0.2) 

ps.markers <- FindMarkers(esc.data2, ident.1 = "PS", ident.2 = "hESC", min.pct = 0.20, 
                          group.by = "clusters", logfc.threshold = 0.2)
#write.csv(ps.markers,"PS_markers_vs_hESC.csv")
TFs_only_in_PS_vs_hESC <- ps.markers[rownames(ps.markers) %in% tflist$V1, ]
write.csv(TFs_only_in_PS_vs_hESC,"min0.20_TFs_in_PS_comparedto_hESC.csv")

ms.markers <- FindMarkers(esc.data2, ident.1 = "ME", ident.2 = "hESC", min.pct = 0.20, 
                          group.by = "clusters", logfc.threshold = 0.2)
#write.csv(ms.markers,"MS_markersvshESC.csv")

TFs_only_in_ME_vs_hESC <- ms.markers[rownames(ms.markers) %in% tflist$V1, ]
write.csv(TFs_only_in_ME_vs_hESC,"min0.20_TFs_in_ME_comparedto_hESC.csv")

de.markers <- FindMarkers(esc.data2, ident.1 = "DE", ident.2 = "hESC", min.pct = 0.20, 
                          group.by = "clusters", logfc.threshold = 0.2)
#write.csv(de.markers,"DE_markers_vs_hESC.csv")
TFs_only_in_DE_vs_hESC <- de.markers[rownames(de.markers) %in% tflist$V1, ]

write.csv(TFs_only_in_DE_vs_hESC,"min_0.20_TFs_in_DE_comparedto_hESC.csv")

# unk.markers <- FindMarkers(esc.data2, ident.1 = "UNK", ident.2 = "hESC", min.pct = 0, 
#                            group.by = "clusters", logfc.threshold = 0.2)

ms.ps.markers <- FindMarkers(esc.data2, ident.1 = "ME", ident.2 = "PS", min.pct = 0.20, 
                             group.by = "clusters", logfc.threshold = 0.2)
de.ps.markers <- FindMarkers(esc.data2, ident.1 = "DE", ident.2 = "PS", min.pct = 0.20, 
                             group.by = "clusters", logfc.threshold = 0.2)

de.ms.markers <- FindMarkers(esc.data2, ident.1 = "DE", ident.2 = "ME", min.pct = 0.20, 
                             group.by = "clusters", logfc.threshold = 0.2)

ms.de.markers <- FindMarkers(esc.data2, ident.1 = "ME", ident.2 = "DE", min.pct = 0.20, 
                             group.by = "clusters", logfc.threshold = 0.2)

esc.markers$genes <- rownames(esc.markers)
ps.markers$genes <- rownames(ps.markers)
ms.markers$genes <- rownames(ms.markers)
de.markers$genes <- rownames(de.markers)
#unk.markers$genes <- rownames(unk.markers)
de.ps.markers$genes <- rownames(de.ps.markers)
ms.ps.markers$genes <- rownames(ms.ps.markers)
de.ms.markers$genes <- rownames(de.ms.markers)
ms.de.markers$genes <- rownames(ms.de.markers)

ps.markers <- ps.markers %>% filter(ps.markers$p_val_adj < 0.05)
ms.markers <- ms.markers %>% filter(ms.markers$p_val_adj < 0.05)
de.markers <- de.markers %>% filter(de.markers$p_val_adj < 0.05)
de.ps.markers <- de.ps.markers %>% filter(de.ps.markers$p_val_adj < 0.05)
ms.ps.markers <- ms.ps.markers %>% filter(ms.ps.markers$p_val_adj < 0.05)
de.ms.markers <- de.ms.markers %>% filter(de.ms.markers$p_val_adj < 0.05)
ms.de.markers <- ms.de.markers %>% filter(ms.de.markers$p_val_adj < 0.05)

humanTF <- scan("../../humanTF.txt", character(), quote = "")

```
```{r}
ps.markers <- ps.markers %>% filter(genes %in% humanTF) %>% filter(avg_log2FC  > 0) #18
write.csv(ps.markers, "min020TF_PS_vs_hESC.csv")
ms.markers <- ms.markers %>% filter(genes %in% humanTF) %>% filter(avg_log2FC  > 0) #35
write.csv(ms.markers, "min020TF_ME_vs_hESC.csv")

de.markers <- de.markers %>% filter(genes %in% humanTF) %>% filter(avg_log2FC  > 0) #45
write.csv(de.markers, "min020TF_DE_vs_hESC.csv")


mes.markers <- intersect(ms.markers$genes, ps.markers$genes) #"TBXT"  "MSX1"  "TBX6"  "MIXL1" "FOXH1" "SP5"
write.csv(mes.markers, "intersect_ME_PS.csv")

end.markers <- intersect(de.markers$genes, ps.markers$genes) #"LHX1"   "OTX2"   "GATA6"  "EOMES"  "GSC"    "HHEX"   "PLSCR1" "TGIF1"  "MIXL1"
write.csv(end.markers,"intersect_DE_PS.csv")

common.markers <- intersect(de.markers$genes, ms.markers$genes) #"GTF2I"  "ZNF90"  "SALL4"  "ARID3A" "ZNF292" "NR6A1"  "BPTF"   "MIXL1"
write.csv(common.markers, "intersect_DE_MS.csv")

common.sig.markers <- intersect(mes.markers, end.markers) #"MIXL1"  
write.csv(common.sig.markers, "intersect_ME_PS_and_DE_PS.csv")

mes.sig.markers <- setdiff(mes.markers, common.sig.markers) #"TBXT"  "MSX1"  "TBX6"  "FOXH1" "SP5" 
de.sig.markers <- setdiff(end.markers, common.sig.markers)


############################## 
data <- subset(esc.data2, subset = clusters != "UNK")
                 #   (  0" "1" "2" "3" "4" "5" "7" "9")

new.cluster.ids <- c("ME", "ME","hESC","DE","hESC","DE","DE","PS")

names(new.cluster.ids) <- levels(data)
data <- RenameIdents(data, new.cluster.ids)
DimPlot(data, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(data,"clusters_reassigned.rds")
##############################

### save as loom file for Scenic analysis 
#esc.data2.loom <- as.loom("clusters_reassigned.rds", filename = "Sox17Ng_aggr_36h_72hDE_72hM_filtered_R.loom", verbose = FALSE)
#esc.data2.loom$close_all()