################################################################################################
#Title: Gene Set Enrichment Analysis 
#Author: Meltem Omur 
################################################################################################
#Setup: 
library("msigdbr")
library("escape")
library("PairedData")
library("Seurat")

#clusters_reassigned.rds is the output of Seurat analysis.
sce <-readRDS("..../clusters_reassigned.rds") 
exprMatrix <- as.matrix(sce@assays$RNA@data)

GS.hallmark <- getGeneSets(library = "H")

ES.seurat <- enrichIt(obj = exprMatrix, gene.sets = GS.hallmark, groups = 1000, cores = 3)

# Add gene set enrichments to Seurat object-metadata. 
sce <- Seurat::AddMetaData(sce, ES.seurat)
#saveRDS(sce, "sce_with_halmarks.rds")

violin <- multi_dittoPlot(sce, vars = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_MTORC1_SIGNALING","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_MYC_TARGETS_V1","HALLMARK_HYPOXIA","HALLMARK_PI3K_AKT_MTOR_SIGNALING"), 
                          group.by = "clusters", plots = c("boxplot"),  nrow=3,
                          ylab = "Enrichment Scores", color.panel = c("#dc050c","#1965b0","#f7f056","#4eb265")+stat_pvalue_manual(stat.test,label = "p.val", step.increase = 0.1,y.position = 35)
,
                          theme = theme_classic()+theme(plot.title = element_text(size = 8),axis.text.x =  element_text(size = 8), axis.text.y = element_text(size = 8))+stat_compare_means())
ggsave(file="boxplot_colored_with_tsne_all_hallmark_paths.svg", plot=violin, width=10, height=10)


# Subset weight data before treatment
before <- subset(sce,  group == "DE", drop = TRUE)
# subset weight data after treatment
after <- subset(my_data,  group == "after", weight,
                drop = TRUE)
# Plot paired data
pd <- paired(before, after)
plot(pd, type = "profile") + theme_bw()

WNT_BETA_<- FetchData(object = sce, vars = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_MTORC1_SIGNALING","HALLMARK_TGF_BETA_SIGNALING","HALLMARK_MYC_TARGETS_V1","HALLMARK_HYPOXIA","HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                                             "clusters"))

x<- WNT_BETA_[WNT_BETA_$clusters == 'hESC',]$HALLMARK_WNT_BETA_CATENIN_SIGNALING #mean(x) 0.09278243
y<- WNT_BETA_[WNT_BETA_$clusters == 'DE',]$HALLMARK_WNT_BETA_CATENIN_SIGNALING #mean(y) 0.1085079
z<- WNT_BETA_[WNT_BETA_$clusters == 'ME',]$HALLMARK_WNT_BETA_CATENIN_SIGNALING #mean(z) 0.2001037
res_x_y <- wilcox.test(x,y) # p-value < 2.2e-16
res_x_z <- wilcox.test(x,z) # p-value < 2.2e-16


a<- WNT_BETA_[WNT_BETA_$clusters == 'hESC',]$HALLMARK_MTORC1_SIGNALING
b<- WNT_BETA_[WNT_BETA_$clusters == 'DE',]$HALLMARK_MTORC1_SIGNALING
c<- WNT_BETA_[WNT_BETA_$clusters == 'ME',]$HALLMARK_MTORC1_SIGNALING
res_ab <- wilcox.test(a,b) # hESC vs DE: p-value < 2.2e-16
res_ac <- wilcox.test(a,c) #  hESC vs ME: p-value < 2.2e-16

t<- WNT_BETA_[WNT_BETA_$clusters == 'hESC',]$HALLMARK_TGF_BETA_SIGNALING
g<- WNT_BETA_[WNT_BETA_$clusters == 'DE',]$HALLMARK_TGF_BETA_SIGNALING
f<- WNT_BETA_[WNT_BETA_$clusters == 'ME',]$HALLMARK_TGF_BETA_SIGNALING
res_tg <- wilcox.test(t,g) # hESC vs DE: p-value < 2.2e-16 
res_tf <- wilcox.test(t,f) #  hESC vs ME: p-value < 2.2e-16

m<- WNT_BETA_[WNT_BETA_$clusters == 'hESC',]$HALLMARK_MYC_TARGETS_V1
y<- WNT_BETA_[WNT_BETA_$clusters == 'DE',]$HALLMARK_MYC_TARGETS_V1
c<- WNT_BETA_[WNT_BETA_$clusters == 'ME',]$HALLMARK_MYC_TARGETS_V1
res_my<- wilcox.test(m,y) # hESC vs DE: p-value < 2.2e-16 
res_mc <- wilcox.test(m,c) #  hESC vs ME: p-value < 2.2e-16


hesc<- WNT_BETA_[WNT_BETA_$clusters == 'hESC',]$HALLMARK_HYPOXIA
de<- WNT_BETA_[WNT_BETA_$clusters == 'DE',]$HALLMARK_HYPOXIA
me<- WNT_BETA_[WNT_BETA_$clusters == 'ME',]$HALLMARK_HYPOXIA
res_ed<- wilcox.test(hesc,de) # hESC vs DE: p-value < 2.2e-16 
res_mc <- wilcox.test(hesc,me) #  hESC vs ME: p-value < 2.2e-16


hesc2<- WNT_BETA_[WNT_BETA_$clusters == 'hESC',]$HALLMARK_PI3K_AKT_MTOR_SIGNALING
de2<- WNT_BETA_[WNT_BETA_$clusters == 'DE',]$HALLMARK_PI3K_AKT_MTOR_SIGNALING
me2<- WNT_BETA_[WNT_BETA_$clusters == 'ME',]$HALLMARK_PI3K_AKT_MTOR_SIGNALING
res_ed2<- wilcox.test(hesc2,de2) # hESC vs DE: p-value 0.3437
res_mc2 <- wilcox.test(hesc2,me2) #  hESC vs ME: p-value < 2.2e-16
