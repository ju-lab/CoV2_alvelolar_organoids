# Three-dimensional human alveolar stem cell culture models reveal infection response to SARS-CoV-2

Youk, J., Kim, T., Evans, K.V., Jeong, Y.-I., Hur, Y., Hong, S.P., Kim, J.H., Yi, K., Kim, S.Y., Na, K.J., Bleazard, T., Kim, H.M., Fellows, M., Mahbubani, K.T., Saeb-Parsy, K., Kim, S.Y., Kim, Y.T., Koh, G.Y., Choi, B.-S., Ju, Y.S., Lee, J.-H., Three-dimensional human alveolar stem cell culture models reveal infection response to SARS-CoV-2, Cell Stem Cell (2020), doi: https://doi.org/10.1016/j.stem.2020.10.004.

Thank you for your interest in this paper.
This page provides R scipts to reproduce representative figures for bulk RNA sequencing and single cell RNA sequenicng analysis in this project.

There are three ways to start the repoduction.
- You can download raw fastq files for bulk RNA sequencing and bam rils for single cell RNA sequencing. They have been deposited in the European Genome-Phenome Archive (EGA) with accession ID EGAS00001004508.
- You can download RSEM output files for bulk RNA sequenicng and filtered_feature_bc_matrix files for single cell RNA sequencing. They have been uploaded in https://www.synapse.org/#!Synapse:syn22146555/files/.
- You can depict figures from Seurat/SingleCellExperiment objects for single cell RNA sequencing. They have been uploaded in https://www.synapse.org/#!Synapse:syn22146555/files/.


## 1. Bulk RNA sequencing
Detailed instruciton will be updated soon.

## 2. Single Cell RNA sequencing
There are two files related to scRNAseq analysis; 02_scRNAseq.R and 03_Figure_scRNA.R
If you try to generate Seurat/SingleCellExperiment Objects from feature matrices, please follow "02_scRNAseq.R" file.
If you already generated or downloaded Seurat/SingleCellExperiment Objects, please follow "03_Figure_scRNA.R" file. 




If you have any questions, please send me an e-mail to jhyouk@gmail.com.

Jeonghwan Youk
(last update: Dec. 3rd, 2020)
