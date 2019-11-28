#!/usr/local/bin/Rscript
rm(list=ls())

##############################################################################################################
#   Project: CTR_gjb2_0001_placenta_1st_vs_2nd_trimester   
#                 
#   Malwina Prater (mn367@cam.ac.uk), 2019                     
#   Centre for Trophobast Research, University of Cambridge
#   Script: DMR analysis
#   
##############################################################################################################

Project       <- "CTR_gjb2_0001_Methylation_1500promoter"

significance  <- 0.05
l2fc          <- 1
col_1st <- "firebrick2"
col_2nd <- "steelblue3"

Base.dir      <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/Methylation_data"
setwd(Base.dir)
list.files(Base.dir)


library(ggplot2)
library(cowplot)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(karyoploteR)

message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version="GRCh37")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)          
#ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description', 'go_id', 'name_1006' ), mart = ensembl)
ensEMBL2id_go <- read.csv("ensEMBL2id_go.csv")


message("--------------------------------------------------------------------------------")
message("+    make plot of methylation vs l2fc expression to identify interesting genes  ")
message("+-------------------------------------------------------------------------------")

message("--------------------------------------------------------------------------------")
message("+                 load in tables                        ")
message("+-------------------------------------------------------------------------------")

resSig.ann <- read.table(file.path(Base.dir, "CTR_gjb2_0001_STAR_DESeq2_shrinkage_deseq2_DEGs_padj0.05.csv"), sep = ",", header = TRUE) 
resSig.ann <- resSig.ann[,-1]
DMRs_closest <- read.table(file.path(Base.dir, "EPIC_DifferentialMethylationFiles/Bedtools_closest_wo_DMRs_to_genes_2kb_promoter.txt"), sep = "\t", header = F) 
colnames(DMRs_closest) <- c("chr", "start", "end", "DMR", "chr_2", "start_2", "end_2", "ensembl_gene_id", "distance") 
bedtools_intersect_wao_table <- read.table(file.path(Base.dir, "EPIC_DifferentialMethylationFiles/Bedtools_intersect_wao_DMRs_to_genes_2kb_promoter.txt"), sep = "\t", header = F) 
results.df <- read.csv("CTR_gjb2_0001_STAR_DESeq2_shrinkage_res_all.csv")
colnames(bedtools_intersect_wao_table) <- c("chr", "start", "end", "DMR", "chr2", "start2", "end2", "ensembl_gene_id", "overlap")
rld_df <- read.table(file.path(Base.dir, "Methylation_data_1st_2nd_trimester/CTR_gjb2_0001_STAR_DESeq2_shrinkage_DESeq2-rlog-transformed-count.txt"), sep = "\t", header = TRUE) 

message("--------------------------------------------------------------------------------")
message("+                top DE genes with DMRs                       ")
message("+-------------------------------------------------------------------------------")

DMRs_closest_DEGs <- subset(DMRs_closest, abs(DMRs_closest$distance) < 500) # here are both promoters and very close to promoters (2kb UPSTREAM of gene)
DMRs_closest_DEGs <- DMRs_closest_DEGs[DMRs_closest_DEGs$ensembl_gene_id %in% resSig.ann$ensembl_gene_id,] 
DMRs_to_DEGs <- unique(merge(DMRs_closest_DEGs[,c(4,8)], resSig.ann[,c(1,3,7,8,9,10)], by.x = "ensembl_gene_id"))
DMRs_to_DEGs_0.6 <- subset(DMRs_to_DEGs, abs(DMRs_to_DEGs$log2FoldChange) > 0.6)

df_expr_up_meth_down <- subset(df_l2fc_expression_fc_methylation1500, df_l2fc_expression_fc_methylation1500$l2fc_expression > 0 & df_l2fc_expression_fc_methylation1500$delta_methylation_oxBS < 0 )
df_expr_down_meth_up <- subset(df_l2fc_expression_fc_methylation1500, df_l2fc_expression_fc_methylation1500$l2fc_expression < 0 & df_l2fc_expression_fc_methylation1500$delta_methylation_oxBS > 0 )

df_meth_expr <- rbind(df_expr_down_meth_up, df_expr_up_meth_down)
df_meth_expr <- df_meth_expr[df_meth_expr$ensembl_gene_id %in% DMRs_to_DEGs_0.6$ensembl_gene_id,]

METH_DIFF <- 0.1 
DE_L2FC   <- 1  

df_meth_expr <- df_meth_expr[abs(df_meth_expr$delta_methylation_oxBS) > METH_DIFF & abs(df_meth_expr$l2fc_expression) > DE_L2FC,] 


message("--------------------------------------------------------------------------------")
message("+              make scatterplot                       ")
message("+-------------------------------------------------------------------------------")

bedtools_intersect_ann <- unique(merge(bedtools_intersect_wao_table, results.df[,c(1,2,3,7)], by.x = "ensembl_gene_id", by.y = "X", all.x = T))

genes_with_DMRs <- unique(bedtools_intersect_ann[bedtools_intersect_ann$overlap != ".",]$ensembl_gene_id)
length(unique(bedtools_intersect_ann[bedtools_intersect_ann$overlap != ".",]$ensembl_gene_id)) # 316 unique ens genes with DMRs
length(unique(genes_with_DMRs[genes_with_DMRs %in% resSig.ann$ensembl_gene_id])) # 126 DEGs with DMR
length(unique(genes_with_DMRs[genes_with_DMRs %in% resSig.ann[abs(resSig.ann$log2FoldChange) > 1,]$ensembl_gene_id])) # 62 DEGs with DMRs l2fc>1
length(unique(bedtools_intersect_ann[bedtools_intersect_ann$overlap != ".",]$DMR)) # 264 DMRs within gene promoters



df_expression_methylation1500 <- read.table(file.path(Base.dir, "CTR_gjb2_0001_Methylationdf_expression_methylation1500.csv"), sep = ",", header = TRUE) 
df_expression_methylation1500 <- df_expression_methylation1500[,-1]
head(df_expression_methylation1500)

df_expression_methylation1500$methylation_FC_oxBS <- df_expression_methylation1500$oxBS_2nd-df_expression_methylation1500$oxBS_1st
df_expression_methylation1500$methylation_FC_BS <- df_expression_methylation1500$BS_2nd-df_expression_methylation1500$BS_1st

df_l2fc_expression_fc_methylation1500 <- df_expression_methylation1500[,c(1,6,7,8,9,12,13)]
colnames(df_l2fc_expression_fc_methylation1500)[2] <- "l2fc_expression"
colnames(df_l2fc_expression_fc_methylation1500)[6] <- "delta_methylation_oxBS"
colnames(df_l2fc_expression_fc_methylation1500)[7] <- "delta_methylation_BS"
head(df_l2fc_expression_fc_methylation1500)


selected_genes <- as.character(df_meth_expr$external_gene_name)

plt_expr_meth <- ggplot(data = df_l2fc_expression_fc_methylation1500, aes(x=delta_methylation_oxBS, y= l2fc_expression)) + 
  geom_point(size=1, alpha=0.3,  col="grey") +
  geom_point(data = subset(df_l2fc_expression_fc_methylation1500, df_l2fc_expression_fc_methylation1500$l2fc_expression > 1 & df_l2fc_expression_fc_methylation1500$delta_methylation_oxBS < -0.1 ), size=1, alpha=0.5,  col="steelblue3") +
  geom_point(data = subset(df_l2fc_expression_fc_methylation1500, df_l2fc_expression_fc_methylation1500$l2fc_expression < -1 & df_l2fc_expression_fc_methylation1500$delta_methylation_oxBS > 0.1 ), size=1, alpha=0.5,  col="firebrick2") +
  ylab("log2 Fold Change expression (ST/FT)") + xlab("Methylation change (ST-FT) ") +
  geom_hline(yintercept=0) + geom_vline(xintercept=0) +
  geom_text_repel(aes(label=ifelse(df_l2fc_expression_fc_methylation1500$external_gene_name %in% selected_genes , as.character(external_gene_name),'')) ) + xlim(-0.5, 0.5)



message("+-------------------------------------------------------------------------------")
message("+                       Load gene lists to plot on heatmaps                     ")
message("+-------------------------------------------------------------------------------")

genes_TFs  <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/TF_genes_tf_checkpoint_db.rds")
genes_oxy <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/oxy_genes_ego.rds")
genes_secretory <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/secretory_genes_ego.rds")
genes_ER  <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/ER_genes_stress_kegg_go.rds")
genes_hormone <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/hormone_activity_GO_genes.rds")
genes_Prolif_diff<- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/Prolif_differentiation_GO_genes.rds")
genes_Transport <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/All_Transport_GO_genes.rds")
genes_Lipid_Transport <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/Lipid_Transport_GO_genes.rds")
genes_ECM <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/ECM_genes.rds")
genes_WNT <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/WNT_signaling_resSig.rds")
genes_cell_cycle <- readRDS("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/bulk_vs_scRNA_seq/cell_cycle_genes_resSig.rds")


message("--------------------------------------------------------------------------------")
message("+                prepare table useful for SPLIT in complex heatmaps             ")
message("+-------------------------------------------------------------------------------")

df_meth_expr$Split <- ""
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_TFs, "TF", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_secretory, "Secretion", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_ER, "ER", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_Transport, "Transport", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_ECM, "ECM", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_Prolif_diff, "Proliferation & differentiation", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_cell_cycle$external_gene_name, "Cell cycle", df_meth_expr$Split ) # none!
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_hormone, "Hormone", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% genes_WNT$external_gene_name, "WNT", df_meth_expr$Split ) # none!±
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% c("ADAMS19", "STAMBPL1", "CNDP1", "MMP3"), "Metalopeptidase", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% c("CLDN10"), "Cell-cell adhesion", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% c("EVA1C"), "Glycosaminoglycan binding", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% c("ASTN1"), "Cell-cell adhesion", df_meth_expr$Split )
df_meth_expr$Split <- ifelse(df_meth_expr$external_gene_name %in% c("SH3RF3"), "Protein ubiquitination", df_meth_expr$Split )

df_meth_expr_for_heatmap <- subset(df_meth_expr, abs(df_meth_expr$l2fc_expression) > 1.0)
df_meth_expr_for_heatmap <- subset(df_meth_expr_for_heatmap, abs(df_meth_expr_for_heatmap$delta_methylation_oxBS) > 0.1)



message("--------------------------------------------------------------------------------")
message("+                 complex heatmaps - expression vs methylation                  ")
message("+-------------------------------------------------------------------------------")

rld_df$ensembl_gene_id <- rownames(rld_df)
rld_df_ann <- unique(merge(rld_df, ensEMBL2id[,c(1,2)], by = "ensembl_gene_id", all.x = TRUE))
rownames(rld_df_ann) <- rld_df_ann$ensembl_gene_id 
head(rld_df_ann)

#    rld means:
rld_df_ann$FT <- rowMeans(rld_df_ann[,c(2:9)])
rld_df_ann$ST <- rowMeans(rld_df_ann[,c(10:15)])
rld_df_means <- rld_df_ann[,-c(1:15)]

#    MeanCentred - means:
rld_df_meanCentered <- rld_df_ann[,c(2:15)] - rowMeans(rld_df_ann[,c(2:15)])   
rld_df_meanCentered$FT <- rowMeans(rld_df_meanCentered[,c(2:9)])
rld_df_meanCentered$ST <- rowMeans(rld_df_meanCentered[,c(10:15)])
rld_df_meanCentered <- rld_df_meanCentered[,-c(1:14)]
rld_df_meanCentered$ensembl_gene_id <- rownames(rld_df_meanCentered)

rld_df_meanCentered_ann <- unique(merge(rld_df_meanCentered, ensEMBL2id[,c(1,2)], by.x = "ensembl_gene_id"))
rownames(rld_df_meanCentered_ann) <- rld_df_meanCentered_ann$ensembl_gene_id
rld_df_meanCentered_ann <- rld_df_meanCentered_ann[,c(4,2,3)]

head(df_expression_methylation1500)
methylation_mat <- unique(df_expression_methylation1500[,c(1:3,8)])
colnames(methylation_mat)[2:3] <- c("FT", "ST")
head(methylation_mat)

rld_mat <- rld_df_meanCentered_ann[rld_df_meanCentered_ann$external_gene_name %in% selected_genes,]
rownames(rld_mat) <- rld_mat$external_gene_name
rld_mat <- rld_mat[,-1]

meth_mat <- methylation_mat[methylation_mat$external_gene_name %in% rownames(rld_mat) , ] 
meth_mat <- meth_mat[match( rownames(rld_mat), meth_mat$external_gene_name),]

rownames(meth_mat) <- rownames(rld_mat)
meth_mat <- (meth_mat[,c(2,3)])
meth_mat$delta <- meth_mat$ST - meth_mat$FT

# add info about the Split- by GO term etc:::
Split1 <- df_meth_expr_for_heatmap[df_meth_expr_for_heatmap$external_gene_name %in% rownames(rld_mat),]
rld_mat$Split <- Split1$Split[match(rownames(rld_mat), Split1$external_gene_name )]
unique(rld_mat$Split)
rld_mat$Split <- factor(rld_mat$Split, levels = c("Oxidative stress" , "Proliferation & differentiation" ,"Transport" ,  "WNT",   "Cell-cell adhesion" ,  "Glycosaminoglycan binding" , "Protein ubiquitination" , ""  ))

f1 = colorRamp2( c(-2,-0.8,0,0.8,2), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
f2 = colorRamp2( c(0,0.25,0.5,0.75,1), c("darkgreen", "green3", "grey95", "slateblue", "darkorchid4"), space = "RGB") 
lgd1 = Legend(col_fun = f2, title = "Expression", at =  c(-2,-1,0,1,2))
lgd2 = Legend(col_fun = f2, title = "Methylation", at = c(0, 0.5, 1)) 
pd = packLegend(lgd1, lgd2, direction = "vertical")

ht1 = Heatmap(rld_mat[,c(1,2)],  col = f1, name = "Expression",  row_title = "", column_title = "Expression", show_row_names = T, heatmap_legend_param = list(title = "Expression", legend_height = unit(3, "cm"), title_position = "topleft"), cluster_columns = FALSE, row_names_side ="left",  split = (rld_mat$Split), row_title_rot = 0, show_row_dend = F) 
ht2 = Heatmap(meth_mat[,c(1,2)],   col = f2, name = "Methylation", row_title = "", column_title = "Methylation",  heatmap_legend_param = list( at = c( 0, 0.5, 1), labels = c("0", "0.5", "1"), title = "Methylation", legend_height = unit(3, "cm"), title_position = "topleft"), split = meth_mat$meth_reg, cluster_columns = FALSE, show_row_names = FALSE)
ht_list = ht1 + ht2
ht_list

pdf(paste("Fig_6x", Project,  "__FT_ST_gene_expression_methylation"  , "__selected_GO_groups__", "ComplexHeatmap", ".pdf", sep="_"), onefile=FALSE, width=9, height=5) 
par(bg=NA)
draw(ht_list, row_title = " ", row_title_gp = gpar(col = "red"),
     column_title = " ", column_title_side = "bottom", gap = unit(1, "cm"))
dev.off()






message("+-------------------------------------------------------------------------------")
message("+                                karyoplotes                                    ")
message("+-------------------------------------------------------------------------------")

resSig_meth.ann <- resSig.ann[resSig.ann$ensembl_gene_id %in% df_meth_expr$ensembl_gene_id,]

UP_l2fc1 <- subset(resSig.ann, resSig.ann$log2FoldChange > 1)
DOWN_l2fc1 <- subset(resSig.ann, resSig.ann$log2FoldChange < -1)

txdb_ens <- makeTxDbFromUCSC(genome="hg19", tablename="ensGene")
columns(txdb_ens)
all.genes <- genes(txdb_ens, columns=c("CDSNAME" , "TXNAME", "GENEID"))

gene.UP.symbols <- UP_l2fc1[UP_l2fc1$ensembl_gene_id %in% df_meth_expr$ensembl_gene_id,]$external_gene_name
genes.UP        <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
                                   filters = 'external_gene_name', values =gene.UP.symbols, mart = ensembl))
seqlevelsStyle(genes.UP) <- "UCSC"

gene.DN.symbols <- DOWN_l2fc1[DOWN_l2fc1$ensembl_gene_id %in% df_meth_expr$ensembl_gene_id,]$external_gene_name
genes.DN        <- toGRanges(getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'external_gene_name'),
                                   filters = 'external_gene_name', values =gene.DN.symbols, mart = ensembl))
seqlevelsStyle(genes.DN) <- "UCSC"

kp <- plotKaryotype(genome="hg19", plot.type=2)
kpDataBackground(kp, color = "#FFFFFF00", data.panel=2)
kp <- kpPlotDensity(kp, all.genes)
kpPlotRegions(kp, data=genes.UP, data.panel=2, col="red",  r0=-1.0, r1=0.5, lwd=2)
kpPlotRegions(kp, data=genes.DN, data.panel=2, col="blue", r0=-1.0, r1=0.5, lwd=2)


pdf(paste("Fig.6.c", "__DEG_l2fc1_SigGenes_meth0.1_karyoplot.pdf", sep=""), width=8, height=10, onefile=FALSE)
par(bg=NA)
kp <- plotKaryotype(genome="hg19", plot.type=2)
kpDataBackground(kp, color = "#FFFFFF00", data.panel=2)
kp <- kpPlotDensity(kp, all.genes)
kpPlotRegions(kp, data=genes.UP, data.panel=2, col="red",  r0=-1.0, r1=0.5, lwd=2)
kpPlotRegions(kp, data=genes.DN, data.panel=2, col="blue", r0=-1.0, r1=0.5, lwd=2)
dev.off()


sessionInfo()
