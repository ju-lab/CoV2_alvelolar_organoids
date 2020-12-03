################################################
## Three-Dimensional Human Alveolar Stem Cell Culture Models Reveal Infection Response to SARS-CoV-2
## Cell Stem Cell, October 21, 2020, DOI:https://doi.org/10.1016/j.stem.2020.10.004
## 
## Scripts for Figure generation of single cell RNA sequencing
## If you want to follow the script to generate R data files (.RDS), please see 02_scRNAseq.R in https://github.com/ju-lab/SARS-CoV-2_alveolar_organoids
##
## Scripts have been posted in https://github.com/ju-lab/SARS-CoV-2_alveolar_organoids
## Feature_bc_matrix data and Seurat/SingleCellExperiment objects have been uploaded in https://www.synapse.org/#!Synapse:syn22146555/files/
## Bam files for the five single cell data have been deposited in the European Genome-Phenome Archive (EGA) with accession ID EGAS00001004508 
## 3rd Dec. 2020. by Jeonghwan Youk, Su Yeon Kim
################################################

##-------------------------------------------
## Directory setting
##-------------------------------------------
#working_dir <- "/home/users/jhyouk/13_SARS_COV_2/github_script"
working_dir <- "/PATH_To_workding_directory"
setwd(working_dir)

##-------------------------------------------
## Load libraries and source files
##-------------------------------------------
require(RColorBrewer)
require(Seurat)

source("func_SC3_CoreFunctions.R")
source("func_SC3_ShinyFunctions.R")



##-------------------------------------------
## Load final seurat object and final single cell experiment object
## You can download final seurat object and final single cell experiment object in https://www.synapse.org/#!Synapse:syn22146555/files/
## Or, if you want to generate two above objects from filtered_feature_bc_matrix, please see 02_scRNAseq.R in https://github.com/ju-lab/SARS-CoV-2_alveolar_organoids
## Or, if you want to download bam files, please download the files from the European Genome-Phenome Archive (EGA) with accession ID EGAS00001004508 
##-------------------------------------------
sroF <- readRDS("SeuratObj_RefitV3.rmlowSF.rmBlood_UMAPpc30_clust0.1_FINAL_cleaned.RDS")
sceF <- readRDS("SCE_rmlowSF.rmBlood_g5k_SC3_k9_trained_FINAL.RDS")                    

## rename clusters for brief name 
sceF$SC3_clusters <- factor(sceF$sc3_9_clusters)
sceF$Seurat_clusters <- sceF$ClusterrefitV3Fac
sceF$Experiment <- sceF$Condition
sceF$Level_SARSCoV2 <- sceF$`Level-SARS-CoV-2`


##-------------------------------------------
## colors for annotations
##-------------------------------------------

cols.cond <- c("olivedrab2", "goldenrod2", "orangered1", "red3", "violetred3")
cols.clusterV3 <- c(brewer.pal(12, "Paired")[c(2, 4)],
                    "darkorange1",
                    brewer.pal(12, "Paired")[c(12)],
                    brewer.pal(12, "Paired")[c(10)],
                    "magenta", "navyblue")
cols.sc3_clusters <- c(brewer.pal(8, "Accent")[1:7], brewer.pal(12, "Paired")[c(10, 6)])

cols.Experiment <- cols.cond
names(cols.Experiment) <- levels(sceF$Experiment)

cols.Seurat_clusters <- cols.clusterV3
names(cols.Seurat_clusters) <- levels(sceF$Seurat_clusters)

cols.SC3_clusters <- c(brewer.pal(8, "Accent")[1:7], brewer.pal(12, "Paired")[c(10, 6)])
names(cols.SC3_clusters) <- levels(sceF$SC3_clusters)

cols.levelSARS <- c("#CCCCCC","#A088BF","#7765AD","brown3")
cols.Level_SARSCoV2 <- cols.levelSARS
names(cols.Level_SARSCoV2) <- levels(sceF$Level_SARSCoV2)


annot_colors <- list(Experiment=cols.Experiment,
                     Seurat_clusters=cols.Seurat_clusters,
                     SC3_clusters=cols.SC3_clusters,
                     Level_SARSCoV2=cols.Level_SARSCoV2)






##==================================================
## Make figures
##==================================================

##--------------------------------------------------
## Figure 6A
##--------------------------------------------------
DimPlot(sroF, reduction = "umap", 
        group.by = "Condition", cols=cols.cond)

##--------------------------------------------------
## Figure 6B
##--------------------------------------------------
DimPlot(sroF, reduction = "umap", 
        group.by = "ClusterrefitV3Fac", cols=cols.clusterV3)

##--------------------------------------------------
## Figure 6C
##--------------------------------------------------
MyBoxplot  <- function (dat, xvar, yvar, f.col, myCols, 
                        mycex=1.0, plotLg=FALSE, showaxis=TRUE,...) {
  
  if(class(dat[,xvar])!="factor") dat[,xvar] <- factor(dat[,xvar])
  if(class(dat[,f.col])!="factor") dat[,f.col] <- factor(dat[,f.col])
  
  ## discard those without counts
  b <- boxplot(dat[,yvar] ~ dat[,xvar], xaxt='n', 
               ...,
               outline=FALSE)
  abline (v=1:length(levels(dat[,xvar])), lty=3, col=gray(0.9))
  if(showaxis) axis(1, 1:length(levels(dat[,xvar])), levels(dat[,xvar]), las=2, cex.axis=1)
  points ( jitter(as.numeric(dat[,xvar]), 1),  dat[,yvar],
           col=myCols[as.numeric(dat[,f.col])], 
           pch=19, cex=mycex)
  
  ct <- table(dat[,f.col])
  if(plotLg) legend("topleft", levels(dat[,f.col]),
                    col=myCols, pch=19, cex=0.9,
                    bg="transparent", bty='n')
  b
}
CreateBoxplot <- function (gdat, gname, f.group, cols.group, ymax) {
  
  b <- MyBoxplot (gdat, xvar=f.group, yvar=gname, 
                  ylim=c(0, ymax),
                  f.col=f.group, myCols=cols.group, mycex=0.4, 
                  plotLg=FALSE, showaxis=FALSE, 
                  xlab="", ylab="Expression level",
                  main=gname)
  text(1:length(b$names), par('usr')[3], labels=b$names, srt=60,
       adj=c(1.1, 1.0), xpd = TRUE, cex=1.0)
  
  grmed <- b$stats[3,] 
  points(grmed, type='o', pch=16, cex=2, lwd=2, col="dodgerblue3")
  points(grmed, pch=1, cex=2.1, lwd=1.5, col="black", lty=2)
  
  grmean <- tapply(gdat[,gname], gdat[,f.group], mean)
  points(grmean, type='o', pch=18, cex=2, lwd=2, col="gold2", lty=2)
  points(grmean, pch=1, cex=2.1, lwd=2, col="black", lty=2)
  
  ct <- table(gdat[,f.group]); ct
  ct0 <- tapply(gdat[,gname], gdat[,f.group], function (v) sum(v==0, na.rm=T))
  
}

gdat <- sroF@meta.data
gname <- "SARS-CoV-2"
f.group <- "Condition"
cols.group <- cols.cond


CreateBoxplot (gdat, gname, f.group, cols.group, ymax=20) 
abline(h=c(1, 3, 9), lty=2, col=gray(0.8))


##--------------------------------------------------
## Figure 6D
##--------------------------------------------------

DimPlot(sroF, reduction = "umap", 
        group.by = "Level-SARS-CoV-2", cols=cols.levelSARS,
        pt.size = 0.85)


##--------------------------------------------------
## Figure 6E
##--------------------------------------------------
#sroF@meta.data$'log2-Total-UMI' <- log2(sroF@meta.data$`Total-UMI`)
FeaturePlot(sroF, reduction = "umap", feature="log2-Total-UMI")
FeaturePlot(sroF, reduction = "umap", feature="SFTPB")
FeaturePlot(sroF, reduction = "umap", feature="NKX2-1")
FeaturePlot(sroF, reduction = "umap", feature="ACE2")
FeaturePlot(sroF, reduction = "umap", feature="TMPRSS2")

summary(sroF@meta.data$`Total-UMI`[sroF@meta.data$ClusterrefitV3Fac=="rcl3_60h_m1.0_highCoV2.mix"])
summary(sroF@meta.data$`Total-UMI`[sroF@meta.data$ClusterrefitV3Fac!="rcl3_60h_m1.0_highCoV2.mix"])

tapply(sroF@assays$RNA@scale.data["SFTPB",],sroF@meta.data$ClusterrefitV3Fac,mean)
tapply(sroF@assays$RNA@scale.data["NKX2-1",],sroF@meta.data$ClusterrefitV3Fac,mean)

table(sroF@assays$RNA@scale.data["ACE2",]>0)[2]/(table(sroF@assays$RNA@scale.data["ACE2",]>0)[2]+table(sroF@assays$RNA@scale.data["ACE2",]>0)[1])
table(sroF@assays$RNA@scale.data["TMPRSS2",]>0)[2]/(table(sroF@assays$RNA@scale.data["TMPRSS2",]>0)[2]+table(sroF@assays$RNA@scale.data["TMPRSS2",]>0)[1])

##--------------------------------------------------
## Figure 7A
## This figure can be different whenever SC3 clustering was newly performed.
## However, key features on each SC3 cluster are usually consistent, such as AQP5, HSPA1A, and IFI27.
##--------------------------------------------------

sc3_plot_markers.revised  <- function(object, k, auroc, p.val, show_pdata, annot_colors) {
  
  if (is.null(metadata(object)$sc3$consensus)) {
    warning(paste0("Please run sc3_consensus() first!"))
    return(object)
  }
  hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
  dataset <- get_processed_dataset(object)
  if (!is.null(metadata(object)$sc3$svm_train_inds)) {
    dataset <- dataset[, metadata(object)$sc3$svm_train_inds]
  }
  
  add_ann_col <- FALSE
  ann <- NULL
  if (!is.null(show_pdata)) {
    ann <- make_col_ann_for_heatmaps(object, show_pdata)
    if (!is.null(ann)) {
      add_ann_col <- TRUE
      # make same names for the annotation table
      rownames(ann) <- colnames(dataset)
    }
  }
  
  # get all marker genes
  markers <- organise_marker_genes(object, k, p.val, auroc)
  
  if(!is.null(markers)) {
    # get top 10 marker genes of each cluster
    markers <- markers_for_heatmap(markers)
    
    row.ann <- data.frame(Cluster = factor(markers[, 1], levels = unique(markers[, 1])))
    rownames(row.ann) <- markers$feature_symbol
    
    do.call(pheatmap::pheatmap, 
            c(list(dataset[markers$feature_symbol, , drop = FALSE], 
                   show_colnames = FALSE, 
                   cluster_rows = FALSE, 
                   cluster_cols = hc, 
                   cutree_cols = k, 
                   annotation_row = row.ann, 
                   annotation_names_row = FALSE, 
                   gaps_row = which(diff(markers[, 1]) != 0), 
                   cellheight = NA),  ## set cellheight as NA
              list(annotation_col = ann)[add_ann_col],
              list(annotation_colors=annot_colors)[add_ann_col])) ## !! added
  } else {
    message("No markers have been found, try to lower significance thresholds!")
  }
}


sc3_plot_markers.revised(sceF, k=9, 
                         auroc = 0.8, p.val = 0.01,                # In the paper, Figure 7B was illustrated with auroc = 0.85
                         show_pdata = c("Level_SARSCoV2", 
                                        "Experiment", 
                                        "Seurat_clusters", 
                                        "SC3_clusters"),
                         annot_colors)


##--------------------------------------------------
## Figure 7B
##--------------------------------------------------
FeaturePlot(sroF, reduction = "umap", feature="AQP5")
FeaturePlot(sroF, reduction = "umap", feature="HSPA1A")
FeaturePlot(sroF, reduction = "umap", feature="IFI27")
FeaturePlot(sroF, reduction = "umap", feature="PPP1R15A")

summary(sroF@assays$RNA@scale.data["AQP5",][sroF@meta.data$ClusterrefitV3Fac=="rcl0_0h_m0.0"])
summary(sroF@assays$RNA@scale.data["AQP5",][sroF@meta.data$ClusterrefitV3Fac!="rcl0_0h_m0.0"])

summary(sroF@assays$RNA@scale.data["HSPA1A",][sroF@meta.data$ClusterrefitV3Fac=="rcl4_4h_m0.1"])
summary(sroF@assays$RNA@scale.data["HSPA1A",][sroF@meta.data$ClusterrefitV3Fac!="rcl4_4h_m0.1"])

summary(sroF@assays$RNA@scale.data["IFI27",][sroF@meta.data$ClusterrefitV3Fac=="rcl2_60h_m0.1"])
summary(sroF@assays$RNA@scale.data["IFI27",][sroF@meta.data$ClusterrefitV3Fac=="rcl1_60h_m1.0_lowCoV2"])
summary(sroF@assays$RNA@scale.data["IFI27",][sroF@meta.data$ClusterrefitV3Fac!="rcl1_60h_m1.0_lowCoV2" & sroF@meta.data$ClusterrefitV3Fac!="rcl2_60h_m0.1" ])

summary(sroF@assays$RNA@scale.data["PPP1R15A",][sroF@meta.data$ClusterrefitV3Fac=="rcl3_60h_m1.0_highCoV2.mix"])
summary(sroF@assays$RNA@scale.data["PPP1R15A",][sroF@meta.data$ClusterrefitV3Fac!="rcl3_60h_m1.0_highCoV2.mix"])



##--------------------------------------------------
## Figure S6A
##--------------------------------------------------

sc3_plot_consensus.revised <- function(object, k, show_pdata, annot_colors) {
  if (is.null(metadata(object)$sc3$consensus)) {
    warning(paste0("Please run sc3_consensus() first!"))
    return(object)
  }
  hc <- metadata(object)$sc3$consensus[[as.character(k)]]$hc
  consensus <- metadata(object)$sc3$consensus[[as.character(k)]]$consensus
  
  add_ann_col <- FALSE
  ann <- NULL
  if (!is.null(show_pdata)) {
    ann <- make_col_ann_for_heatmaps(object, show_pdata)
    if (!is.null(ann)) {
      add_ann_col <- TRUE
      # make same names for the annotation table
      rownames(ann) <- colnames(consensus)
    }
  }
  do.call(pheatmap::pheatmap, 
          c(list(consensus, cluster_rows = hc, cluster_cols = hc, cutree_rows = k, 
                 cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), 
            list(annotation_col = ann)[add_ann_col],
            list(annotation_colors=annot_colors)[add_ann_col])) ## !! added
}


sc3.consensus.plot <- sc3_plot_consensus.revised(sceF, k=9, 
                                                 show_pdata = c("Level_SARSCoV2", 
                                                                "Experiment", 
                                                                "Seurat_clusters", 
                                                                "SC3_clusters"),
                                                 annot_colors)

sc3.consensus.plot

##--------------------------------------------------
## Figure S6B
##--------------------------------------------------
FeaturePlot(sroF, reduction = "umap", feature="SOX2")
FeaturePlot(sroF, reduction = "umap", feature="KRT5")
FeaturePlot(sroF, reduction = "umap", feature="TP63")

table(sroF@assays$RNA@scale.data["SOX2",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[2]/(table(sroF@assays$RNA@scale.data["SOX2",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[2]+table(sroF@assays$RNA@scale.data["SOX2",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[1])
table(sroF@assays$RNA@scale.data["SOX2",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[2]/(table(sroF@assays$RNA@scale.data["SOX2",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[2]+table(sroF@assays$RNA@scale.data["SOX2",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[1])

table(sroF@assays$RNA@scale.data["KRT5",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[2]/(table(sroF@assays$RNA@scale.data["KRT5",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[2]+table(sroF@assays$RNA@scale.data["KRT5",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[1])
table(sroF@assays$RNA@scale.data["KRT5",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[2]/(table(sroF@assays$RNA@scale.data["KRT5",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[2]+table(sroF@assays$RNA@scale.data["KRT5",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[1])

table(sroF@assays$RNA@scale.data["TP63",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[2]/(table(sroF@assays$RNA@scale.data["TP63",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[2]+table(sroF@assays$RNA@scale.data["TP63",][sroF@meta.data$ClusterrefitV3Fac=="rcl6_airway-like"]>0)[1])
table(sroF@assays$RNA@scale.data["TP63",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[2]/(table(sroF@assays$RNA@scale.data["TP63",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[2]+table(sroF@assays$RNA@scale.data["TP63",][sroF@meta.data$ClusterrefitV3Fac!="rcl6_airway-like"]>0)[1])


##--------------------------------------------------
## Figure S6F
##--------------------------------------------------

FeaturePlot(sroF, reduction = "umap", feature="MKI67")
table(sroF@assays$RNA@scale.data["MKI67",]>0)[2]/(table(sroF@assays$RNA@scale.data["MKI67",]>0)[2]+table(sroF@assays$RNA@scale.data["MKI67",]>0)[1])
