#loading required packages
library(limma) 
library(edgeR)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tibble)
library(pheatmap)
library(sva)
library(fgsea)
library(stats)
library(gplots)
library(RSvgDevice)
library(stringr)
library(statmod)
library(pcaMethods)
library(fastICA)
library(RSvgDevice)
library(cairoDevice)
library(ggplot2)
library(tximportData)
library(tximport)
library(tximeta)
library(GenomicFeatures)
library(data.table)
library(msigdbr)

#Loading the datasets and preparing the data input for ssCMAP/LINCS
v96v0_75748 <- readRDS("v96v0_75748.rds")
v96v0_109658 <- readRDS("v96v0_109658.rds")

#getting up-reg genes
v96v0_75748_pos <- v96v0_75748 %>% filter(logFC > 2) %>%
  filter(adj.P.Val < 0.05)
v96v0_109658_pos <- v96v0_109658 %>% filter(logFC > 2) %>%
  filter(adj.P.Val < 0.05)

#getting down-reg genes
v96v0_75748_neg <- v96v0_75748 %>% filter(logFC < -2) %>%
  filter(adj.P.Val < 0.05)
v96v0_109658_neg <- v96v0_109658 %>% filter(logFC < -2) %>%
  filter(adj.P.Val < 0.05)

#common up- and down-reg genes
inter_pos <- intersect(rownames(v96v0_75748_pos), rownames(v96v0_109658_pos))
inter_neg <- intersect(rownames(v96v0_75748_neg), rownames(v96v0_109658_neg))
length(inter_pos)
length(inter_neg)

com_pos_75748 <- v96v0_75748_pos[inter_pos,]
com_pos_109658 <- v96v0_109658_pos[inter_pos,]

com_neg_75748 <- v96v0_75748_neg[inter_neg,]
com_neg_109658 <- v96v0_109658_neg[inter_neg,]

com_pos <- merge(x = com_pos_75748, y = com_pos_109658, by = "gene", all = TRUE)
com_neg <- merge(x = com_neg_75748, y = com_neg_109658, by = "gene", all = TRUE)

all_com <- rbind(com_pos, com_neg)
colnames(all_com) <- c("gene","logFC.75748","AveExpr.75748","t.75748",
                       "P.Value.75748","adj.P.Val.75748","B.75748","logFC.109658",
                       "AveExpr.109658","t.109658","P.Value.109658","adj.P.Val.109658",
                       "B.109658")
write.table(all_com, file="Supplementary_table_2.tsv", sep="\t", row.names = FALSE)

######################################################################################
v96v0_75748$genes <- rownames(v96v0_75748)
v96v0_109658$genes <- rownames(v96v0_109658)

common_up <- intersect(v96v0_75748Up$genes, v96v0_109658Up$genes)
common_down <- intersect(v96v0_75748Down$genes, v96v0_109658Down$genes)

lfc109658 <- v96v0_109658$logFC
names(lfc109658) <- v96v0_109658$genes
lfc109658UP <- lfc109658[common_up]
lfc109658DOWN <- lfc109658[common_down]

lfc75748 <- v96v0_75748$logFC
names(lfc75748) <- v96v0_75748$genes
lfc75748UP <- lfc75748[common_up]
lfc75748DOWN <- lfc75748[common_down]

dfUP <- data.frame(lfc75748UP, lfc109658UP)
dfDOWN <- data.frame(lfc75748DOWN, lfc109658DOWN)

#for selecting the genes
dfUP <- apply(dfUP, 1,function(x) { max(x) })
dfDOWN <- apply(dfDOWN, 1,function(x) { min(x) })

res <- c(dfUP, dfDOWN)

write(names(res[1:150]), "list_of_up_reg_genes_75748+109658_NEW_180221.txt")
write(names(tail(res, n =150)), "list_of_dowm_reg_genes_75748+109658_NEW_180221.txt")

######################################################################################
#converting to gene id
library(hgu133a.db)
library(annotate)

x <- hgu133aSYMBOL

mapped_probes <- mappedkeys(x)
# Convert to a dataframe
genesym.probeid <- as.data.frame(x[mapped_probes])

head(genesym.probeid)

positive_probes <- genesym.probeid %>% filter(symbol %in% names(res[1:150]))
positive_probes$regulation <- 1
negative_probes <- genesym.probeid %>% filter(symbol %in% names(tail(res, n =150)))
negative_probes$regulation <- -1

probes <- c(positive_probes$probe_id, negative_probes$probe_id)
regulation <- c(positive_probes$regulation, negative_probes$regulation)

for_sscMAP <- data.frame(AffyProbeSetID = probes, regulation = regulation)

write.table(for_sscMAP,"for_sscMAP_75748_109658_180221.sig",sep="\t",row.names=FALSE, quote=F)