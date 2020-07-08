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


Project       <- "CTR_gjb2_0001_STAR_DESeq2_shrinkage___bulk_vs_scRNA_seq_BlockSex"
col_1st <- "firebrick2"
col_2nd <- "steelblue3"

Res.dir      <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex"
setwd(Res.dir)


message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")
ensembl    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)          
ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'go_id', 'name_1006' ), mart = ensembl)    
ensEMBL2id <- unique(ensEMBL2id_go[,-c(1,6,7)])


message("+-------------------------------------------------------------------------------")
message("+                       Load bulk rna seq data                           ")
message("+-------------------------------------------------------------------------------")

rld_df <- as.data.frame(assay(rld))
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


rld_df_meanCentered_ann <- readRDS("rld_df_meanCentered_ann.Rds")
GSE89497_mat_means <- readRDS("GSE89497_mat_means.Rds")

resSig.ann <- read.csv("CTR_gjb2_0001_STAR_DESeq2_shrinkage_BlockSex_deseq2_DEGs_padj0.05.csv")
DEGs_pre_postflow <- resSig.ann

resdata_ann <- read.csv("CTR_gjb2_0001_STAR_DESeq2_shrinkage_BlockSex_deseq2_resdata_ann.csv")

message("+-------------------------------------------------------------------------------")
message("+                       Load sc rna seq data                           ")
message("+-------------------------------------------------------------------------------")

GSE89497_scaledata <- as.data.frame(readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/GSE89497_matrix.su__scale.data_FINAL.rds")) 
GSE89497_metadata <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/GSE89497_matrix.su__meta.data.rds")
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
dim(GSE89497_mat_means) 
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

genes_Transport <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex/rds_BlockSex/All_Transport_GO_genes.rds")
genes_Lipid_Transport <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex/rds_BlockSex/Lipid_Transport_GO_genes.rds")
genes_ECM <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex/rds_BlockSex/ECM_genes_BlockSex.rds")


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


# transport with SPLIT :::
ensEMBL2id_go_transport <- ensEMBL2id_go[grepl("transport", ensEMBL2id_go$name_1006),]
genes_Transport_selected <- genes_Transport
genes_Transport_selected <- rbind(genes_Transport_selected, resSig.ann[resSig.ann$external_gene_name == "LTF",])
genes_Transport_selected <- rbind(genes_Transport_selected, resSig.ann[resSig.ann$external_gene_name == "HEPH",])

genes_Transport_selected$Split <- ""
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "ion transmembrane transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "ion transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("ABC", genes_Transport_selected$external_gene_name), "Lipid transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "lipid transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Lipid transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("APO", genes_Transport_selected$external_gene_name), "Apolipoproteins",  genes_Transport_selected$Split)

genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "glucose", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Glucose transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "protein transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Protein transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "amino acid transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Amino acid transport",  genes_Transport_selected$Split)

genes_Transport_selected$Split <- ifelse( grepl("cation", genes_Transport_selected$description), "Cation transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("sodium", genes_Transport_selected$description), "Cation transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("potassium", genes_Transport_selected$description), "Cation transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("calcium", genes_Transport_selected$description), "Calcium transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "oxygen transport", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Oxygen transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("transferrin", genes_Transport_selected$description), "Iron & metal ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( genes_Transport_selected$external_gene_name %in% ensEMBL2id_go_transport[grepl( "metal ion", ensEMBL2id_go_transport$name_1006),]$external_gene_name, "Iron & metal ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("SLC30A2", genes_Transport_selected$external_gene_name), "Iron & metal ion transport",  genes_Transport_selected$Split)
genes_Transport_selected$Split <- ifelse( grepl("SLC30A10", genes_Transport_selected$external_gene_name), "Iron & metal ion transport",  genes_Transport_selected$Split)

genes_Transport_selected <- subset(genes_Transport_selected, genes_Transport_selected$Split != "")
genes_Transport_selected <- unique(genes_Transport_selected)


message("+-------------------------------------------------------------------------------")
message("+                     makeHeatmap function without split                        ")
message("+-------------------------------------------------------------------------------")

makeHeatmap <- function(selected_genes, name_of_heatmap, gene_number = 50, by_cellType = FALSE){
  selected_genes <- unique(DEGs_pre_postflow[DEGs_pre_postflow$external_gene_name %in% selected_genes,])
  selected_genes <- selected_genes[order(abs(selected_genes$log2FoldChange), decreasing = TRUE),] 
  selected_genes <- unique(selected_genes$external_gene_name)
  if (length(selected_genes) > gene_number )  {
    selected_genes <- selected_genes[1:gene_number]
  }
  rld_mat <- rld_df_meanCentered_ann[rld_df_meanCentered_ann$external_gene_name %in% selected_genes,]
  rownames(rld_mat) <- rld_mat$external_gene_name
  rld_mat <- rld_mat[,-1]
  rld_mat2 <- rld_mat
  scRNAseq_mat <- GSE89497_mat_means[match( rownames(rld_mat2), rownames(GSE89497_mat_means)),]
  
  scRNAseq_mat2 <- scRNAseq_mat
  rownames(scRNAseq_mat2) <- rownames(rld_mat2)
 
  head(rld_mat2)
  head(scRNAseq_mat2)
  colnames(scRNAseq_mat2) <- c("STB_8W", "CTB_8W", "EVT_8W", "EVT_24W", "STR_8W")
  rownames(scRNAseq_mat2) <- rownames(rld_mat2)
  scRNAseq_mat2[is.na(scRNAseq_mat2)] <- -10
  
  if (by_cellType == TRUE){
    # this is for ordering instead of clustering - by bulk and scrna seq!!!
    rld_mat2$scRNA_celltype <- colnames(scRNAseq_mat2)[apply(scRNAseq_mat2,1,which.max)]
    rld_mat2$bulkRNA        <- colnames(rld_mat2)[apply(rld_mat2[,c(1,2)],1,which.max)]
    rld_mat2$scRNA_celltype <- factor(rld_mat2$scRNA_celltype, levels= c("STB_8W", "CTB_8W", "EVT_8W","STR_8W", "EVT_24W" ))
    rld_mat2 <- rld_mat2[order(rld_mat2$scRNA_celltype),]
    rld_mat2 <- rld_mat2[order(rld_mat2$bulkRNA),]
    scRNAseq_mat2 <- scRNAseq_mat2[match( rownames(rld_mat2), rownames(scRNAseq_mat2)),]
    head(rld_mat2)
    head(scRNAseq_mat2)
  }

  f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
  f2 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB")
  lgd1 = Legend(col_fun = f1, title = "Expression (rld)", at = c(-2,-1,0,1,2)  )
  lgd2 = Legend(col_fun = f2, title = "sc RNA-seq",  at = c(-2,-1,0,1,2)  )
  pd = packLegend(lgd1, lgd2, direction = "vertical")

  scRNAseq_mat2 <- scRNAseq_mat2[,c(1,2,3,5,4)]
  
  scRNAseq_mat2[scRNAseq_mat2== -10] <- NA
  
  
  if (by_cellType == TRUE){
    ht1 = Heatmap(as.matrix(rld_mat2[,c(1,2)]),  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = TRUE, heatmap_legend_param = list(title = "Expression (rld)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = FALSE , width = unit(3, "cm"), row_names_side ="left")
  }else{
    ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "Expression (rld)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(3, "cm"), row_names_side ="left")
  }  
  
  ht2 = Heatmap(as.matrix(scRNAseq_mat2),   col = f2, name = "sc RNA-seq", row_title = "", column_title = "sc RNA-seq", show_row_names = FALSE, heatmap_legend_param = list( title = "sc RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(4, "cm"), column_split = c(rep("08 w",4), "24 w"), show_heatmap_legend = FALSE) 
  ht_list = ht1 + ht2
  
  heatmap_height <- nrow(rld_mat)*0.2+1
  pdf(paste("Fig_xx_", "ComplexHeatmap_",  name_of_heatmap, ".pdf", sep=""), onefile=FALSE, width=8, height=heatmap_height) 
  par(bg=NA)
  draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),  column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
  dev.off()
  
  print(nrow(rld_mat2))
  print(rownames(rld_mat2))
  return(ht_list)
}


genes_ht_ER <- c("ERP44","ANK1","HSPA5","HERPUD1","KIF2A","FCGR2B","DNAJB11","DNAJC3","TUSC3","AREG","MAN1C1","SDF2L1","DNAJB9","PDIA6","MANF","HYOU1","KIF5A","TGFA","SPTA1","HSP90B1","PDIA3","HSPA6","GOLT1A","CALR","MCFD2","SSR4","GAS6","PROS1","PDIA2","LRRK2","SERPINA1","ERO1A","MAN1A2","SCAMP5","FICD","XBP1","LMAN1", "ERN1", "XBP1","EIF2AK3", "EDEM2", "ANK2" ) 
genes_ht_secret <- c("SLC22A16","ABCC8","CD22","IGF1","FCGR2B","EDN1","TLR8","TNFSF13B","RGCC","PIK3CG","C5","ADTRP","SLC16A10","PLA2G4A","INHBA","INHA","CCR7","WNK4","GAD1","C5AR2","CGA","TLR2","FGF7","FCGR2A","ADRA2A","MCOLN2","SYN2","TNFRSF14","APOA2","FCGR3B","SLC30A8","TLR1","LEP","MCTP1","NLRP10","MAFA","C5AR1","SCAMP5","CR1","ANG","CARD17") 
genes_ht_oxy <- c("TBXAS1","ADCY2","EDN1","ICAM1","RGCC","TRPA1","CCL2","NCF2","CCR7","KLF2","CBFA2T3","HBZ","IL10","TLR2","CH25H","CYP1A1","CYP1A2","PKLR","KCNMB1","GUCY1A2","KCNMA1","NCF1","GDNF","TLR6","LEP","APOLD1","PENK","MAFA","CYP2R1","MAPT","P2RX2","DPYD","HBG2","SLC30A10","HIF3A","ANG","COX4I2","SOD3","ACADL", "LOX", "CAT") 
genes_ht_WNT <-c("DKK3","LRP6","TNRC6C","PORCN","FGF10","WNT2","WNT5B","TBX18","IGFBP2","SDC1","WNT4","PDE6B","WNT10A","TLR2","PLCB2","LEF1","NKD1","ARRB2","IGFBP4","NKD2","DCDC2","RSPO3","DRD2","PTPRO","WNT3A","LYPD6","APCDD1","DKK2","EDA","PRKAA2","ZEB2","CTNND2","WNT10B","CDH2","SOX7","LRRK2","PSMB8","PSMB9", "WNT9A", "RYR2", "CCND1",  "RECK","SCEL","APOE","NR4A2","FOLR1")         
genes_ht_cellcycl <-c("CDK14","MLXIPL","FHL1","FOXC1","DUSP13","PTPRC","ABCB1","CDKL1","FAM83D","RGCC","TUBB4A","LILRB1","RASSF4","CCL2","HYAL1","PLAGL1","INHBA","ESX1","INHA","IL10","CGREF1","CYP1A1","MAPK4","PRKACB","HMGA2","PPP1R1C","PPP2R2B","PYHIN1","GEM","GFI1B","FAM107A","AFAP1L2","WNT10B","LEP","RNF212","GAS6","EVI2B","PSMB8","AIF1","SPDYC","GMNC","PSMB9","TUBB3","SMIM22","TGFB1", "GKN1") 
genes_ht_TFs <-c("MLXIPL","FOXC1","HES2","GATA1","MYEF2","KLF1","TFEC","LHX6","ZBTB16","TBX18","HAND1","TFAP2E","KLF7","TCF21","KLF9","HOXC13","SPDEF","TCF15","KLF2","MYOD1","CASZ1","EHF","HEYL","STAT4","RORC","SCML4","ZNF157","NFIB","HMGA2","NPAS3","EBF1","GFI1B","SOX14","EMX2","SOX7","BNC2","MSC","TSHZ2","MAFA","MYT1","SOX18","FOXD1","TCF4","PBX4","MEF2C","GPER1") 
genes_hormone <-c("RORC","SPX","CGA","LEP","CGB2","INHA","CGB5","NR2F1","CGB8","FBN1","PENK","CGB7","CGB3","IGF1","KL","LHB","INSL4","EDN1","METRN","CGB1","RARB","INHBA","PTH2R","GPER1","CSH1","ERFE","GUCA2A","STC1","CSHL1","TTR","RLN1","UTS2B","CRH","NPPB","RXRG","PGR", "NR6A1","NR4A2","EDN2","EGFR", "HCRTR2", "CRYM","RLN2", "NR1D1", "CH25H") 
glycolytic_genes <- c("HK3",  "GCK",  "PFKM", "PGM1" ,"PFKP" ,"HK2",  "PKLR", "CBFA2T3", "IGF1", "MLXIPL", "PFKFB2", "GPD1", "ACADL")  
PSG_genes <- c("Psg16", "Ceacam3", "Psg29","Psg-ps1","Ceacam5","Ceacam14", "Gm5155", "Ceacam11", "Ceacam13", "Ceacam12","Psg18", "Psg28", "Psg26","Psg25","Psg27","Psg23","Psg21","Psg20", "Psg22","Sycp1-ps1", "Psg19","Psg17" )
CTS_genes <- c("Tpbpb", "Tpbpa", "Ctsj","Ctsq","Ctsr","Cts6","Cts8","Cts7","Ctsm","Cts3","Fbp2","Fbp1","Cdc14b","Gm16907","Ctsl")
PRL_genes <- c("Prl","Prl3d1","Prl3d2","Prl3d3","Prl3c1","Prl3b1","Prl3a1","Prl6a1","Prl8a2","Prl2b1","Prl8a6","Prl8a8","Prl8a9","Prl7b1","Prl7a2","Prl7d1","Prl7c1","Prl2a1","Prl2c1","Prl4a1","Prl5a1")


ht_ER       <- makeHeatmap(genes_ht_ER,      name_of_heatmap= "genes_ER",        by_cellType = FALSE) # 41
ht_oxy      <- makeHeatmap(genes_ht_oxy,     name_of_heatmap= "genes_oxy",       by_cellType = FALSE) # 41
ht_hormone  <- makeHeatmap(genes_hormone,    name_of_heatmap= "genes_hormone",     by_cellType = FALSE) # 41
ht_WNT      <- makeHeatmap(genes_ht_WNT,         name_of_heatmap= "genes_WNTxx",         by_cellType = FALSE) # 46
ht_cellcycl <- makeHeatmap(genes_ht_cellcycl,    name_of_heatmap= "genes_cell_cyclex",  by_cellType = FALSE) # 46
ht_TFs      <- makeHeatmap(genes_ht_TFs,         name_of_heatmap= "genes_TFs",         by_cellType = FALSE) # 46
ht_glycol   <- makeHeatmap(glycolytic_genes,                     name_of_heatmap= "glycolytic_genes") # 13
genes_Lipid_Transport <- subset(genes_Lipid_Transport, abs(genes_Lipid_Transport$log2FoldChange) > 1)
ht_lipid_transport  <- makeHeatmap(genes_Lipid_Transport$external_gene_name,    name_of_heatmap= "genes_lipidTransport",     by_cellType = FALSE) # 50



message("+-------------------------------------------------------------------------------")
message("+                     makeHeatmap TRANSPORT genes                               ")
message("+-------------------------------------------------------------------------------")

name_of_heatmap <- "transport_genes"
l2fcx <- 1.0
gene_number <- 64

genes_Transport_selected <- genes_Transport_selected[abs(genes_Transport_selected$log2FoldChange) > l2fcx,]

genes_Transport_selected$Split2 <- genes_Transport_selected$Split
genes_Transport_selected[genes_Transport_selected$external_gene_name == "MCOLN2",]$Split2 <-  "Ion transport"
genes_Transport_selected[genes_Transport_selected$external_gene_name == "HEPH",]$Split2 <-  "Iron & metal ion transport"
genes_Transport_selected$Split2 <- gsub("Cation transport", "Calcium & ion transport", genes_Transport_selected$Split2)
genes_Transport_selected$Split2 <- gsub("Calcium transport", "Calcium & ion transport", genes_Transport_selected$Split2)
genes_Transport_selected$Split2 <- gsub("Ion transport", "Calcium & ion transport", genes_Transport_selected$Split2)


genes_Transport_selected <- genes_Transport_selected[order(abs(genes_Transport_selected$log2FoldChange), decreasing = TRUE),]
unique(genes_Transport_selected$Split2)
final_genes1 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Oxygen transport")
final_genes2 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Calcium & ion transport")[c(1:15),]
final_genes3 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Lipid transport")[c(1:6,8,9),]
final_genes4 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Iron & metal ion transport")[c(1:6),]
final_genes5 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Protein transport")[c(1:4),]
final_genes6 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Apolipoproteins")[c(1:5),]
final_genes7 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Amino acid transport")[c(1:5),]
final_genes8 <- subset(genes_Transport_selected, genes_Transport_selected$Split2 == "Glucose transport")[c(1:4),]
final_genes <- unique(c(as.character(final_genes1$external_gene_name), as.character(final_genes2$external_gene_name), as.character(final_genes3$external_gene_name), as.character(final_genes4$external_gene_name), as.character(final_genes5$external_gene_name), as.character(final_genes6$external_gene_name), as.character(final_genes7$external_gene_name), as.character(final_genes8$external_gene_name)))


selected_genes <- genes_Transport_selected[genes_Transport_selected$external_gene_name %in% final_genes,]
selected_genes <- selected_genes[order(abs(selected_genes$log2FoldChange), decreasing = TRUE),] 

selected_genes <- unique(selected_genes$external_gene_name)

rld_mat <- rld_df_meanCentered_ann[rld_df_meanCentered_ann$external_gene_name %in% selected_genes,]
rownames(rld_mat) <- rld_mat$external_gene_name
rld_mat <- rld_mat[,-1]

scRNAseq_mat <- GSE89497_mat_means[match( rownames(rld_mat), rownames(GSE89497_mat_means)),]
head(rld_mat)
head(scRNAseq_mat)

scRNAseq_mat2 <- scRNAseq_mat
rld_mat2 <- rld_mat
rownames(scRNAseq_mat2) <- rownames(rld_mat2)


colnames(scRNAseq_mat2) <- c("STB_8w", "CTB_8w", "EVT_8w", "EVT_24w", "STR_8w")
head(scRNAseq_mat2)
rld_mat2$scRNA_celltype <- colnames(scRNAseq_mat2)[apply(scRNAseq_mat2,1,which.max)]
rld_mat2$bulkRNA        <- colnames(rld_mat2)[apply(rld_mat2[,c(1,2)],1,which.max)]
rld_mat2$scRNA_celltype <- factor(rld_mat2$scRNA_celltype, levels= c("STB_8w", "CTB_8w", "EVT_8w", "EVT_24w", "STR_8w"))
rld_mat2 <- rld_mat2[order(rld_mat2$scRNA_celltype),]
rld_mat2 <- rld_mat2[order(rld_mat2$bulkRNA),]
scRNAseq_mat2 <- GSE89497_mat_means[match( rownames(rld_mat2), rownames(GSE89497_mat_means)),]
head(rld_mat2)
head(scRNAseq_mat2)
colnames(scRNAseq_mat2) <- c("STB_8w", "CTB_8w", "EVT_8w", "EVT_24w", "STR_8w")

Split1 <- genes_Transport_selected[genes_Transport_selected$external_gene_name %in% rownames(rld_mat),]
rld_mat2$Split <- Split1$Split2[match(rownames(rld_mat2), Split1$external_gene_name )]
rld_mat2$Split <- factor(rld_mat2$Split, levels = c("Glucose transport","Protein transport","Amino acid transport" ,"Lipid transport", "Apolipoproteins","Iron & metal ion transport","Calcium & ion transport", "Oxygen transport"))

scRNAseq_mat2 <- scRNAseq_mat2[,c(1,2,3,5,4)]


f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
f2 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB")
lgd1 = Legend(col_fun = f1, title = "Bulk RNA-seq", at = c(-2,-1,0,1,2)  )
lgd2 = Legend(col_fun = f2, title = "sc RNA-seq",  at = c(-2,-1,0,1,2)  )
pd = packLegend(lgd1, lgd2, direction = "vertical")

ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "Expression (rld)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE, width = unit(3, "cm"), row_names_side ="left", split = rld_mat2$Split, row_title_rot = 0, show_row_dend= TRUE, cluster_row_slices = FALSE, row_names_gp = gpar(fontsize = 12))
ht2 = Heatmap(scRNAseq_mat2,   col = f2, name = "sc RNA-seq", row_title = "", column_title = "sc RNA-seq", show_row_names = FALSE, heatmap_legend_param = list( title = "sc RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(4, "cm"),show_heatmap_legend = FALSE,column_split = c(rep("08 w",4), "24 w"))
ht_list = ht1 + ht2

heatmap_height <- nrow(rld_mat)*0.18
pdf(paste("Fig_xx_", "ComplexHeatmap_",  name_of_heatmap, "3.pdf", sep=""), onefile=FALSE, width=9, height=heatmap_height) 
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

genes_ECM_selected <- rbind(genes_ECM, resSig.ann[resSig.ann$external_gene_name == "MSC",])
genes_ECM_selected <- rbind(genes_ECM_selected, resSig.ann[resSig.ann$external_gene_name == "ANG",])
genes_ECM_selected <- rbind(genes_ECM_selected, resSig.ann[resSig.ann$external_gene_name == "ENG",])
genes_ECM_selected$Split <- ""
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "extracellular matrix structural constituent", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "ECM structural",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "extracellular matrix disassembly", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "ECM disassembly",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "integrin-mediated signaling pathway", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Cell-matrix adhesion",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "cell migration", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Cell migration",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go[grepl( "leukocyte migration", ensEMBL2id_go$name_1006),]$external_gene_name, "Immunomodulation",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go[grepl( "cellular response to leukemia inhibitory factor", ensEMBL2id_go$name_1006),]$external_gene_name, "Immunomodulation",  genes_ECM_selected$Split)

genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "cell-matrix adhesion", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Cell-matrix adhesion",  genes_ECM_selected$Split)

genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "metallo", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Proteases",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go_ECM[grepl( "oxidoreductase activity", ensEMBL2id_go_ECM$name_1006),]$external_gene_name, "Oxidoreductase",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse( genes_ECM_selected$external_gene_name %in% ensEMBL2id_go[grepl( "vasculature development", ensEMBL2id_go$name_1006),]$external_gene_name, "Vasculature",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "SMOC2", "Vasculature",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "PECAM1", "Vasculature",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "VWF", "Vasculature",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "ICAM1", "Vasculature",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "ICAM2", "Immunomodulation",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "ANG", "Vasculature",  genes_ECM_selected$Split)
genes_ECM_selected$Split <- ifelse(genes_ECM_selected$external_gene_name == "ENG", "Vasculature",  genes_ECM_selected$Split)

genes_ECM_selected <- subset(genes_ECM_selected, genes_ECM_selected$Split != "")

MMP_genes <- genes_ECM_selected$external_gene_name[grepl( "MMP", genes_ECM_selected$external_gene_name)]


name_of_heatmap <- "ECM_genes"
l2fcx <- 1.055
gene_number <- 55

genes_ECM_selected <- genes_ECM_selected[order(abs(genes_ECM_selected$log2FoldChange), decreasing = TRUE),] 
genes_ECM_selected <- subset(genes_ECM_selected, abs(genes_ECM_selected$log2FoldChange) > l2fcx | genes_ECM_selected$external_gene_name %in% MMP_genes)
selected_genes <- unique(genes_ECM_selected$external_gene_name)



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

scRNAseq_mat2 <- scRNAseq_mat2[,c(1,2,3,5,4)]


f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
f2 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB")
lgd1 = Legend(col_fun = f1, title = "Expression (rld)", at = c(-2,-1,0,1,2)  )
lgd2 = Legend(col_fun = f2, title = "sc RNA-seq",  at = c(-2,-1,0,1,2)  )
pd = packLegend(lgd1, direction = "vertical")


ht1 = Heatmap(rld_mat2[,c(1,2)],  col = f1, name = "Bulk RNA-seq",  row_title = "", column_title = "Bulk RNA-seq", show_row_names = T, heatmap_legend_param = list(title = "Expression (rld)", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, cluster_rows = TRUE, width = unit(3, "cm"), row_names_side ="left", split = rld_mat2$Split, row_title_rot = 0, show_row_dend= TRUE, cluster_row_slices = FALSE, row_names_gp = gpar(fontsize = 12))
ht2 = Heatmap(scRNAseq_mat2,   col = f2, name = "sc RNA-seq", row_title = "", column_title = "sc RNA-seq", show_row_names = FALSE, heatmap_legend_param = list( title = "sc RNA-seq", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, width = unit(4, "cm"),show_heatmap_legend = FALSE,column_split = c(rep("08 w",4), "24 w"))
ht_list = ht1 + ht2


heatmap_height <- nrow(rld_mat)*0.20
pdf(paste("Fig_xx_", "ComplexHeatmap_",  name_of_heatmap, "2.pdf", sep=""), onefile=FALSE, width=10, height=heatmap_height) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()


