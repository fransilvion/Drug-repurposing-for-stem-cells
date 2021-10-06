# Description 

This repo contains R and bash scripts to reproduce the results of our paper ***In silico* discovery of small molecules for efficient stem cell differentiation into definitive endoderm** (Novakovsky & Sasaki et al.)

1. GSE75748_analysis.Rmd - differential gene expression analysis of the public dataset ([GSE75748](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75748));
2. GSE109658_analysis.Rmd - differential gene expression analysis of the public dataset ([GSE109658](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109658));
3. The following scripts include analysis of the single-cell data:
   1. Seurat_processing.R;
   2. GSEA_analysis_singlecell.R;
   3. Monocole_trajectory_analysis.R;
   4. PyScenic_regulon_analysis.sh;
   5. Pyscenic_visualization.R
4. Data_processing_for_ssCMAP_LINCS.R - get datasets for "ssCMAP/LINCS" approach;
5. cTRAP_analysis.Rmd - includes LINCS part of the "ssCMAP/LINCS" approach and full "Correlation" approach;
6. Pathway_TF_approach.Rmd - includes full "Pathway/TF" approach;
7. RNA_seq_analysis_of_different_treatments.Rmd - differential gene expression analysis of the identified compounds (BRD/JNJ);
