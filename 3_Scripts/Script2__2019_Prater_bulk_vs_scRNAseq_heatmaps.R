#!/usr/local/bin/Rscript
rm(list=ls())

##############################################################################################################
#   Project: CTR_gjb2_0001_placenta_1st_vs_2nd_trimester   
#                 
#   Malwina Prater (mn367@cam.ac.uk), 2019                     
#   Centre for Trophobast Research, University of Cambridge
#   Script: Make heatmaps for bulk vs scRNA-seq datasets for selected enriched terms in preflow vs postflow 
#   
##############################################################################################################

message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(ggrepel)
  library(biomaRt)
  library(DESeq2)
  library(ggdendro)
  library(ggalt)
  library(dplyr)
  library(pheatmap)
  library(enrichR)
  library(reshape2)
  library(pathview)
  library(gage)
  library(gageData)
  library(DOSE)
  library(clusterProfiler)
  library(enrichR)
  library(ComplexHeatmap)
  library(circlize)
  library(Seurat)
})


Project       <- "CTR_gjb2_0001_STAR_DESeq2_shrinkage___bulk_vs_scRNA_seq"
col_1st <- "firebrick2"
col_2nd <- "steelblue3"

Base.dir      <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR"
Res.dir  <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq"
setwd(Res.dir)


message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")
ensembl    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)          
#ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006' ), mart = ensembl)    
ensEMBL2id_go <- read.csv("ensEMBL2id_go.csv")


message("+-------------------------------------------------------------------------------")
message("+                       Load bulk rna seq data                           ")
message("+-------------------------------------------------------------------------------")

rld_df <- read.table(file.path("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/Methylation_data/Methylation_data_1st_2nd_trimester/CTR_gjb2_0001_STAR_DESeq2_shrinkage_DESeq2-rlog-transformed-count.txt"), sep = "\t", header = TRUE) 
rld_df$ensembl_gene_id <- rownames(rld_df)
rld_df_ann <- unique(merge(rld_df, ensEMBL2id[,c(1,2)], by.x = "ensembl_gene_id"))
rownames(rld_df_ann) <- rld_df_ann$ensembl_gene_id 
head(rld_df_ann)

#    rld means:
rld_df_ann$FT <- rowMeans(rld_df_ann[,c(2:9)])
rld_df_ann$ST <- rowMeans(rld_df_ann[,c(10:15)])
rld_df_means <- rld_df_ann[,-c(1:15)]
head(rld_df_means)

#    MeanCentred - means:
rld_df_meanCentered <- rld_df_ann[,c(2:15)] - rowMeans(rld_df_ann[,c(2:15)])   
rld_df_meanCentered$FT <- rowMeans(rld_df_meanCentered[,c(1:8)])
rld_df_meanCentered$ST <- rowMeans(rld_df_meanCentered[,c(9:14)])
rld_df_meanCentered <- rld_df_meanCentered[,-c(1:14)]
rld_df_meanCentered$ensembl_gene_id <- rownames(rld_df_meanCentered)
head(rld_df_meanCentered)

rld_df_meanCentered_ann <- unique(merge(rld_df_meanCentered, ensEMBL2id[,c(1,2)], by.x = "ensembl_gene_id"))
rownames(rld_df_meanCentered_ann) <- rld_df_meanCentered_ann$ensembl_gene_id
rld_df_meanCentered_ann <- rld_df_meanCentered_ann[,c(4,2,3)]
head(rld_df_meanCentered_ann)


resSig.ann <- readRDS("Preflow_vs_postflow__resSig.ann.rds")
DEGs_pre_postflow <- resSig.ann


message("+-------------------------------------------------------------------------------")
message("+                       Load sc rna seq data                           ")
message("+-------------------------------------------------------------------------------")
GSE89497_scaledata <- as.data.frame(readRDS("GSE89497_matrix.su__scale.data_FINAL.rds")) #  dims: 24910 1478
GSE89497_metadata <- readRDS("GSE89497_matrix.su__meta.data.rds")
dim(GSE89497_metadata)
dim(GSE89497_scaledata)
GSE89497_metadata <- GSE89497_metadata[match(  colnames(GSE89497_scaledata),rownames(GSE89497_metadata)),]
rownames(GSE89497_metadata) == colnames(GSE89497_scaledata)

GSE89497_metadata$Split <- GSE89497_metadata$Celltype
GSE89497_metadata$Split <- gsub("HE8W_CTB" , "CTB_8W" , GSE89497_metadata$Split)
GSE89497_metadata$Split <- gsub("HE8W_EVT" , "EVT_8W" , GSE89497_metadata$Split)
GSE89497_metadata$Split <- gsub("HE8W_STR" , "STR_8W" , GSE89497_metadata$Split)
GSE89497_metadata$Split <- gsub("HE8W_STB" , "STB_8W" , GSE89497_metadata$Split)
GSE89497_metadata$Split <- gsub("HE24W_EVT" , "EVT_24W" , GSE89497_metadata$Split)

GSE89497_mat_means <- as.data.frame(GSE89497_scaledata)
dim(GSE89497_mat_means) # [1] 21021 1420
GSE89497_mat_means$STB_8W <- rowMeans(GSE89497_mat_means[ , colnames(GSE89497_mat_means) %in% rownames(GSE89497_metadata[GSE89497_metadata$Split == "STB_8W",])  ]  )  
GSE89497_mat_means$CTB_8W <- rowMeans(GSE89497_mat_means[ , colnames(GSE89497_mat_means) %in% rownames(GSE89497_metadata[GSE89497_metadata$Split == "CTB_8W",])  ]  )  
GSE89497_mat_means$ETV_8W <- rowMeans(GSE89497_mat_means[ , colnames(GSE89497_mat_means) %in% rownames(GSE89497_metadata[GSE89497_metadata$Split == "EVT_8W",])  ]  )  
GSE89497_mat_means$EVT_24W <- rowMeans(GSE89497_mat_means[ , colnames(GSE89497_mat_means) %in% rownames(GSE89497_metadata[GSE89497_metadata$Split == "EVT_24W",])  ]  )  
GSE89497_mat_means$STR_8W <- rowMeans(GSE89497_mat_means[ , colnames(GSE89497_mat_means) %in% rownames(GSE89497_metadata[GSE89497_metadata$Split == "STR_8W",])  ]  )  
dim(GSE89497_mat_means) 
GSE89497_mat_means <- GSE89497_mat_means[,c(1472:1476)]


message("+-------------------------------------------------------------------------------")
message("+                       Load gene lists to plot on heatmaps                     ")
message("+-------------------------------------------------------------------------------")

genes_Transport <- readRDS("All_Transport_GO_genes.rds")
genes_Lipid_Transport <- readRDS("Lipid_Transport_GO_genes.rds")
genes_ECM <- readRDS("ECM_genes.rds")


message("+-------------------------------------------------------------------------------")
message("+                 set up log2 fold change cut off for heatmaps                  ")
message("+-------------------------------------------------------------------------------")

significance  <- 0.05
l2fc          <- 1


message("+-------------------------------------------------------------------------------")
message("+                     prepare gene lists for heatmaps                           ")
message("+-------------------------------------------------------------------------------")

genes_Transport       <- unique(DEGs_pre_postflow[DEGs_pre_postflow$external_gene_name %in% genes_Transport,])
genes_Transport <- rbind(genes_Transport, resSig.ann[resSig.ann$external_gene_name == "HBZ" | resSig.ann$external_gene_name == "HBE1" | resSig.ann$external_gene_name == "HBB" | resSig.ann$external_gene_name == "HBG2", ])
genes_Lipid_Transport <- unique(DEGs_pre_postflow[DEGs_pre_postflow$external_gene_name %in% genes_Lipid_Transport,])
genes_ECM             <- unique(DEGs_pre_postflow[DEGs_pre_postflow$external_gene_name %in% genes_ECM,])

genes_Transport       <- subset(genes_Transport, genes_Transport$padj < significance & abs(genes_Transport$log2FoldChange) > l2fc)
genes_Lipid_Transport <- subset(genes_Lipid_Transport, genes_Lipid_Transport$padj < significance & abs(genes_Lipid_Transport$log2FoldChange) > l2fc)
genes_ECM             <- subset(genes_ECM, genes_ECM$padj < significance & abs(genes_ECM$log2FoldChange) > l2fc)



# transport with SPLIT :::
ensEMBL2id_go_transport <- ensEMBL2id_go[grepl("transport", ensEMBL2id_go$name_1006),]
genes_Transport_selected <- genes_Transport
genes_Transport_selected$Split <- ""
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "metal ion", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Metal ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "ion transmembrane transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("APO", genes_Transport_selected$external_gene_name), "Apolipoproteins",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("ABC", genes_Transport_selected$external_gene_name), "Lipid transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("SLC", genes_Transport_selected$external_gene_name), "SLC transporters",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "glucose", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Glucose transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "protein transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Protein transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("cation", genes_Transport_selected$description), "Cation transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("sodium", genes_Transport_selected$description), "Cation transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("transferrin", genes_Transport_selected$description), "Iron transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("potassium", genes_Transport_selected$description), "Cation transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("calcium", genes_Transport_selected$description), "Calcium transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "oxygen transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Oxygen transport",  genes_Transport_selected$Split)
genes_Transport_selected <- subset(genes_Transport_selected, genes_Transport_selected$Split != "")



message("+-------------------------------------------------------------------------------")
message("+                     makeHeatmap function without split                        ")
message("+-------------------------------------------------------------------------------")

makeHeatmap <- function(selected_genes, name_of_heatmap, gene_number = 50, by_cellType = FALSE){
  selected_genes <- unique(DEGs_pre_postflow[DEGs_pre_postflow$external_gene_name %in% selected_genes,])
  selected_genes <- selected_genes[order(abs(selected_genes$log2FoldChange), decreasing = TRUE),] 
  selected_genes <- subset(selected_genes, abs(selected_genes$log2FoldChange) > l2fc)
  selected_genes <- unique(selected_genes$external_gene_name)
  if (length(selected_genes) > gene_number )  {
    selected_genes <- selected_genes[1:gene_number]
  }
  rld_mat <- rld_df_meanCentered_ann[rld_df_meanCentered_ann$external_gene_name %in% selected_genes,]
  rownames(rld_mat) <- rld_mat$external_gene_name
  rld_mat <- rld_mat[,-1]
  rld_mat2 <- rld_mat
  scRNAseq_mat <- GSE89497_mat_means[match( rownames(rld_mat2), rownames(GSE89497_mat_means)),]
  
  scRNAseq_mat2 <- na.omit(scRNAseq_mat)
  rld_mat2 <-  rld_mat2[rownames(rld_mat2) %in% rownames(scRNAseq_mat2),]
  scRNAseq_mat2 <- scRNAseq_mat2[match( rownames(rld_mat2), rownames(scRNAseq_mat2)),]
  head(rld_mat2)
  head(scRNAseq_mat2)
  colnames(scRNAseq_mat2) <- c("STB_8w", "CTB_8w", "ETV_8w", "ETV_24w", "STR_8w")

  if (by_cellType == TRUE){
    # this is for ordering instead of clustering - by bulk and scrna seq!!!
    rld_mat2$scRNA_celltype <- colnames(scRNAseq_mat2)[apply(scRNAseq_mat2,1,which.max)]
    rld_mat2$bulkRNA        <- colnames(rld_mat2)[apply(rld_mat2[,c(1,2)],1,which.max)]
    rld_mat2$scRNA_celltype <- factor(rld_mat2$scRNA_celltype, levels= c("STB_8w", "CTB_8w", "ETV_8w", "ETV_24w", "STR_8w"))
    rld_mat2 <- rld_mat2[order(rld_mat2$scRNA_celltype),]
    rld_mat2 <- rld_mat2[order(rld_mat2$bulkRNA),]
    scRNAseq_mat2 <- GSE89497_mat_means[match( rownames(rld_mat2), rownames(GSE89497_mat_means)),]
    head(rld_mat2)
    head(scRNAseq_mat2)
  }

  f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
  f2 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB")
  lgd1 = Legend(col_fun = f1, title = "Bulk RNA-seq", at = c(-2,-1,0,1,2)  )
  lgd2 = Legend(col_fun = f2, title = "sc RNA-seq",  at = c(-2,-1,0,1,2)  )
  
  pd = packLegend(lgd1, lgd2, direction = "vertical")
  
  
  if (by_cellType == TRUE){
    ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = TRUE, heatmap_legend_param = list(title = "Bulk RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , width = unit(4, "cm"), row_names_side ="left")
  }else{
    ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "Bulk RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(4, "cm"), row_names_side ="left")
  }  
  
  ht2 = Heatmap(scRNAseq_mat2,   col = f2, name = "sc RNA-seq", row_title = "", column_title = "sc RNA-seq", show_row_names = FALSE, heatmap_legend_param = list( title = "sc RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(6, "cm"))
  ht_list = ht1 + ht2
  
  heatmap_height <- nrow(rld_mat)*0.2+1
  #  pdf(paste("Fig_xx", Project, "ComplexHeatmap",  name_of_heatmap, "_na.rm.pdf", sep="_"), onefile=FALSE, width=8, height=heatmap_height) 
  pdf(paste("Fig_xx", Project, "ComplexHeatmap",  name_of_heatmap, "_na.rm_LEGEND_FIXED.pdf", sep="_"), onefile=FALSE, width=8, height=heatmap_height) 
  par(bg=NA)
  draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
  
  print(nrow(rld_mat2))
  print(rownames(rld_mat2))
  return(ht_list)
}


genes_ht_ER <- c("ERP44","ANK1","HSPA5","HERPUD1","KIF2A","FCGR2B","LMAN1","EDEM2","DNAJB11","SRPX","DNAJC3","TUSC3","AREG","MAN1C1","SDF2L1","DNAJB9","THBS1","LMAN1L","PDIA6","MANF","ANK2","HYOU1","KIF5A","TGFA","SPTA1","HSP90B1","PDIA3","HSPA6","GOLT1A","CALR","MCFD2","SSR4","GAS6","PROS1","PDIA2","LRRK2","SERPINA1","ERO1A","MAN1A2","SCAMP5","FICD") 
genes_ht_secret <- c("SLC22A16","ABCC8","CD22","IGF1","FCGR2B","EDN1","TLR8","TNFSF13B","RGCC","PIK3CG","C5","ADTRP","SLC16A10","PLA2G4A","INHBA","INHA","CCR7","WNK4","GAD1","C5AR2","CGA","TLR2","FGF7","FCGR2A","ADRA2A","MCOLN2","SYN2","TNFRSF14","APOA2","FCGR3B","SLC30A8","TLR1","LEP","MCTP1","NLRP10","MAFA","C5AR1","SCAMP5","CR1","ANG","CARD17") 
genes_ht_oxy <- c("TBXAS1","ADCY2","EDN1","ICAM1","RGCC","TRPA1","CCL2","NCF2","AKR1D1","CCR7","KLF2","CBFA2T3","HBZ","IL10","TLR2","CH25H","CYP1A1","CYP1A2","PKLR","KCNMB1","SESN3","GUCY1A2","KCNMA1","NCF1","GDNF","TLR6","LEP","APOLD1","PENK","MAFA","CYP2R1","MAPT","P2RX2","DPYD","APOD","HBG2","SLC30A10","HIF3A","ANG","COX4I2","SOD3") 
genes_ht_WNT <-c("DKK3","CDK14","LRP6","TNRC6C","PORCN","FGF9","WNT2","WNT5B","TBX18","IGFBP2","SDC1","RECK","WNT1","APOE","ZBED3","PDE6B","WNT10A","SCEL","TLR2","PLCB2","LEF1","NKD1","ARRB2","IGFBP4","NKD2","DCDC2","DAAM2","RSPO3","DRD2","PTPRO","WNT3A","PSMA8","APCDD1","DKK2","EDA","PRKAA2","ZEB2","CTNND2","WNT10B","CDH2","SOX7","LYPD6","LRRK2","NRARP","PSMB8","PSMB9") 
genes_ht_cellcycl <-c("MLXIPL","FHL1","FOXC1","DUSP13","PTPRC","ABCB1","CDKL1","FAM83D","RGCC","TUBB4A","LILRB1","RASSF4","CCL2","HYAL1","PLAGL1","TEX11","TEX14","INHBA","ESX1","INHA","IL10","CGREF1","CYP1A1","MAPK4","PRKACB","HMGA2","PPP1R1C","PPP2R2B","PYHIN1","GEM","GFI1B","FAM107A","AFAP1L2","WNT10B","LEP","RNF212","GAS6","EVI2B","PSMB8","AIF1","SPDYC","GMNC","PSMB9","TUBB3","SMIM22","TGFB1") 
genes_ht_TFs <-c("MLXIPL","FOXC1","HES2","GATA1","MYEF2","KLF1","TFEC","LHX6","ZBTB16","TBX18","HAND1","CDX1","TFAP2E","KLF7","TCF21","KLF9","HOXC13","SPDEF","TCF15","KLF2","MYOD1","CASZ1","EHF","HEY2","STAT4","RORC","SCML4","ZNF157","NFIB","HMGA2","NPAS3","EBF1","GFI1B","SOX14","EMX2","SOX7","BNC2","MAF","MSC","TSHZ2","MAFA","MYT1","TCF4","SOX18","SP5","FOXD1") 
genes_hormone <-c("RORC","SPX","CGA","LEP","CGB2","INHA","SLC16A10","CGB5","NR2F1","CGB8","FBN1","PENK","CGB7","CGB3","IGF1","KL","LHB","INSL4","EDN1","METRN","CGB1","RARB","HAMP","OR51E2","INHBA","PTH2R","GPER1","SLC16A2","NR6A1","CSH1","ERFE","GUCA2A","STC1","CSHL1","TTR","RLN1","UTS2B","CRH","NPPB","RXRG","PGR","EPO","NR4A2","EDN2","AMHR2","FSHR","EGFR")
glycolytic_genes <- c("HK3",  "GCK",  "PFKM", "BPGM" ,"PFKP" ,"HK2",  "PKLR")

ht_ER       <- makeHeatmap(genes_ht_ER,      name_of_heatmap= "genes_ER",        by_cellType = TRUE)
ht_secret   <- makeHeatmap(genes_ht_secret,  name_of_heatmap= "genes_secretory", by_cellType = TRUE)
ht_oxy      <- makeHeatmap(genes_ht_oxy,     name_of_heatmap= "genes_oxy",       by_cellType = TRUE)

ht_hormone  <- makeHeatmap(genes_hormone,    name_of_heatmap= "genes_hormone",     by_cellType = TRUE)

ht_WNT      <- makeHeatmap(genes_ht_WNT,         name_of_heatmap= "genes_WNT",         by_cellType = TRUE)
ht_cellcycl <- makeHeatmap(genes_ht_cellcycl,    name_of_heatmap= "genes_cell_cycle",  by_cellType = TRUE)
ht_TFs      <- makeHeatmap(genes_ht_TFs,         name_of_heatmap= "genes_TFs",         by_cellType = TRUE)

ht_glycol   <- makeHeatmap(glycolytic_genes,                     name_of_heatmap= "glycolytic_genes")

ht_lipid_transport  <- makeHeatmap(genes_Lipid_Transport$external_gene_name,    name_of_heatmap= "genes_hormone",     by_cellType = TRUE)





message("+-------------------------------------------------------------------------------")
message("+                     makeHeatmap TRANSPORT genes                               ")
message("+-------------------------------------------------------------------------------")

name_of_heatmap <- "transport_genes"
l2fcx <- 1.0
gene_number <- 60

genes_Transport_selected$Split2 <- genes_Transport_selected$Split

genes_Transport_selected <- genes_Transport_selected[order(abs(genes_Transport_selected$log2FoldChange), decreasing = TRUE),] 
selected_genes <- subset(genes_Transport_selected, abs(genes_Transport_selected$log2FoldChange) > l2fcx)
selected_genes_tmp <- selected_genes[grepl("Ion", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 2)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]
selected_genes_tmp <- selected_genes[grepl("Cation", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 2)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]
selected_genes_tmp <- selected_genes[grepl("SLC", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 1.9)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]
selected_genes_tmp <- selected_genes[grepl("Lipid", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 1.5)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]
selected_genes_tmp <- selected_genes[grepl("Apolipoproteins", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 1.0)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]
selected_genes_tmp <- selected_genes[grepl("Glucose", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 1.1)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]
selected_genes_tmp <- selected_genes[grepl("Calcium transport", selected_genes$Split), ]
selected_genes_tmp <- subset(selected_genes_tmp, abs(selected_genes_tmp$log2FoldChange) < 1.4)
selected_genes <- selected_genes[ ! selected_genes$external_gene_name %in% selected_genes_tmp$external_gene_name, ]

genes_Transport_selected$Split2 <- gsub("Cation transport", "Calcium & ion transport", genes_Transport_selected$Split2)
genes_Transport_selected$Split2 <- gsub("Calcium transport", "Calcium & ion transport", genes_Transport_selected$Split2)
genes_Transport_selected$Split2 <- gsub("Ion transport", "Calcium & ion transport", genes_Transport_selected$Split2)

selected_genes <- unique(selected_genes$external_gene_name)
if (length(selected_genes) > gene_number )  {
  selected_genes <- selected_genes[1:gene_number]
}
rld_mat <- rld_df_meanCentered_ann[rld_df_meanCentered_ann$external_gene_name %in% selected_genes,]
rownames(rld_mat) <- rld_mat$external_gene_name
rld_mat <- rld_mat[,-1]

scRNAseq_mat <- GSE89497_mat_means[match( rownames(rld_mat), rownames(GSE89497_mat_means)),]
head(rld_mat)
head(scRNAseq_mat)

scRNAseq_mat2 <- na.omit(scRNAseq_mat)
rld_mat2 <-  rld_mat[rownames(rld_mat) %in% rownames(scRNAseq_mat2),]
scRNAseq_mat2 <- scRNAseq_mat2[match( rownames(rld_mat2), rownames(scRNAseq_mat2)),]
head(rld_mat2)
head(scRNAseq_mat2)

colnames(scRNAseq_mat2) <- c("STB_8w", "CTB_8w", "ETV_8w", "ETV_24w", "STR_8w")
rld_mat2$scRNA_celltype <- colnames(scRNAseq_mat2)[apply(scRNAseq_mat2,1,which.max)]
rld_mat2$bulkRNA        <- colnames(rld_mat2)[apply(rld_mat2[,c(1,2)],1,which.max)]
rld_mat2$scRNA_celltype <- factor(rld_mat2$scRNA_celltype, levels= c("STB_8w", "CTB_8w", "ETV_8w", "ETV_24w", "STR_8w"))
rld_mat2 <- rld_mat2[order(rld_mat2$scRNA_celltype),]
rld_mat2 <- rld_mat2[order(rld_mat2$bulkRNA),]
scRNAseq_mat2 <- GSE89497_mat_means[match( rownames(rld_mat2), rownames(GSE89497_mat_means)),]
head(rld_mat2)
head(scRNAseq_mat2)

Split1 <- genes_Transport_selected[genes_Transport_selected$external_gene_name %in% rownames(rld_mat),]
rld_mat2$Split <- Split1$Split2[match(rownames(rld_mat2), Split1$external_gene_name )]
rld_mat2$Split <- factor(rld_mat2$Split, levels = c("SLC transporters","Glucose transport","Protein transport","Lipid transport","Apolipoproteins","Iron transport","Calcium & ion transport", "Oxygen transport"))

f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") #floralwhite
f2 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB")
lgd1 = Legend(col_fun = f1, title = "Bulk RNA-seq", at = c(-2,-1,0,1,2)  )
lgd2 = Legend(col_fun = f2, title = "sc RNA-seq",  at = c(-2,-1,0,1,2)  )
pd = packLegend(lgd1, lgd2, direction = "vertical")

ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "Bulk RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE, width = unit(4, "cm"), row_names_side ="left", split = rld_mat2$Split, row_title_rot = 0, show_row_dend= TRUE, row_names_gp = gpar(fontsize = 12))
ht2 = Heatmap(scRNAseq_mat2,   col = f2, name = "sc RNA-seq", row_title = "", column_title = "sc RNA-seq", show_row_names = FALSE, heatmap_legend_param = list( title = "sc RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(5, "cm"))
ht_list = ht1 + ht2

heatmap_height <- nrow(rld_mat)*0.18
pdf(paste("Fig_xx", Project, "ComplexHeatmap",  name_of_heatmap, "_na.rm2.pdf", sep="_"), onefile=FALSE, width=10, height=heatmap_height) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()



message("+-------------------------------------------------------------------------------")
message("+                         makeHeatmap ECM genes                                 ")
message("+-------------------------------------------------------------------------------")

ensEMBL2id_go_ECM <- unique(ensEMBL2id_go[ensEMBL2id_go$external_gene_name %in% genes_ECM$external_gene_name, ])
tmp_tbl <- as.data.frame(table(ensEMBL2id_go_ECM$name_1006))
tmp_tbl <- tmp_tbl[order(tmp_tbl$Freq, decreasing = T),]
tmp_tbl <- subset(tmp_tbl, tmp_tbl$Freq >5)
genes_ECM_selected <- genes_ECM
genes_ECM_selected$Split <- ""
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "extracellular matrix structural constituent", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "ECM structural constituent",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "extracellular matrix disassembly", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "ECM disassembly",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "integrin-mediated signaling pathway", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Integrin-mediated signaling",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "cell migration", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Cell migration",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "leukocyte migration", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Leukocyte migration",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "cell-matrix adhesion", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Cell-matrix adhesion",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "metallo", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Metallopeptidase",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "oxidoreductase activity", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Oxidoreductase",  genes_ECM_selected$Split)
genes_ECM_TFs <- genes_ECM_selected[genes_ECM_selected$external_gene_name %in% genes_TFs,]$external_gene_name
genes_ECM_selected <- subset(genes_ECM_selected, genes_ECM_selected$Split != "")


name_of_heatmap <- "ECM_genes"
l2fcx <- 1.1
gene_number <- 50

genes_ECM_selected <- genes_ECM_selected[order(abs(genes_ECM_selected$log2FoldChange), decreasing = TRUE),] 
selected_genes <- subset(genes_ECM_selected, abs(genes_ECM_selected$log2FoldChange) > l2fcx)
selected_genes <- unique(selected_genes$external_gene_name)

rld_mat <- rld_df_meanCentered_ann[rld_df_meanCentered_ann$external_gene_name %in% selected_genes,]
rownames(rld_mat) <- rld_mat$external_gene_name
rld_mat <- rld_mat[,-1]

scRNAseq_mat <- GSE89497_mat_means[match( rownames(rld_mat), rownames(GSE89497_mat_means)),]
head(rld_mat)
head(scRNAseq_mat)

Split1 <- genes_ECM_selected[genes_ECM_selected$external_gene_name %in% rownames(rld_mat),]
rld_mat$Split <- Split1$Split[match(rownames(rld_mat), Split1$external_gene_name )]

scRNAseq_mat2 <- GSE89497_mat_means[match( rownames(rld_mat), rownames(GSE89497_mat_means)),]
scRNAseq_mat2 <- na.omit(scRNAseq_mat2)
rld_mat2 <-  rld_mat[rownames(rld_mat) %in% rownames(scRNAseq_mat2),]
scRNAseq_mat2 <- scRNAseq_mat2[match( rownames(rld_mat2), rownames(scRNAseq_mat2)),]
head(rld_mat2)
head(scRNAseq_mat2)

rld_mat2$Split <- Split1$Split[match(rownames(rld_mat2), Split1$external_gene_name )]
rownames(rld_mat2)[rownames(rld_mat2) == "PECAM1"] <- "CD31"

f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "LAB") #floralwhite
f2 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "LAB")
lgd1 = Legend(col_fun = f1, title = "Bulk RNA-seq", at = c(min(rld_mat2[,c(1,2)]),0, max(rld_mat2[,c(1,2)]))  )
lgd2 = Legend(col_fun = f2, title = "sc RNA-seq", at = c(min(scRNAseq_mat2, na.rm = TRUE),0, max(scRNAseq_mat2, na.rm = TRUE)) , border = TRUE )
pd = packLegend(lgd1, lgd2, direction = "vertical")

ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "Bulk RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(4, "cm"), row_names_side ="left", split = rld_mat2$Split, row_title_rot = 0, show_row_dend= TRUE, row_names_gp = gpar(fontsize = 12))
ht2 = Heatmap(scRNAseq_mat2,   col = f2, name = "sc RNA-seq", row_title = "", column_title = "sc RNA-seq", show_row_names = FALSE, heatmap_legend_param = list( title = "sc RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(5, "cm"))
ht_list = ht1 + ht2

heatmap_height <- nrow(rld_mat)*0.20
pdf(paste("Fig_xx", Project, "ComplexHeatmap",  name_of_heatmap, "_na.rm2.pdf", sep="_"), onefile=FALSE, width=10, height=heatmap_height) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


