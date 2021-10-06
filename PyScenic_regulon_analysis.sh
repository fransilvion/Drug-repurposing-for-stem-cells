#############################################################################################
#Title: PyScenic regulon analysis
#Author: Meltem Omur 
#############################################################################################
#Setup: 
#conda create -n pyscenic python=3.7
#conda activate pyscenic
#pip install loompy==2.0.17
#pip install pyscenic
#pip install scanpy
#conda install numpy pandas matplotlib seaborn
#conda install -c anaconda cytoolz  
#conda install seaborn scikit-learn statsmodels numba pytables
#conda install -c conda-forge python-igraph louvain
##############################################################################################

# transcription factors:
wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt

# motif to TF annotation database:
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

# genome ranking database:
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

##############			Network inference (GRNBoost2)     ####################################################################
#The first step in the SCENIC pipeline is to perform GRN inference.The primary algorithm used for this is GRNBoost2. 
##############################################################################################################################
pyscenic grn \
--num_workers 4 \
--output adj.tsv \
--method grnboost2 \
Sox17Ng_aggr_36h_72hDE_72hM_filtered_R.loom \
hs_hgnc_tfs.txt


##### Candidate regulon generation and regulon prediction using motifs ####################################################################
pyscenic ctx \
adj.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname Sox17Ng_aggr_36h_72hDE_72hM_filtered_R.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 6 \
--mask_dropouts


#### Cellular enrichment -- AUCell Step  ####################################################################################################
#After the creation of the refined regulons, cellular enrichment can be assayed with AUCell. 
#############################################################################################################################################
 
pyscenic aucell \
Sox17Ng_aggr_36h_72hDE_72hM_filtered_R.loom \
reg.csv \
--output Sox17Ng_aggr_36h_72hDE_72hM_filtered_AUCell.loom \
--num_workers 6