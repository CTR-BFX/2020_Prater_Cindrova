#!/usr/local/bin/Rscript
rm(list=ls())

##############################################################################################################
#   Project: CTR_gjb2_0001_placenta_1st_vs_2nd_trimester   
#                 
#   Malwina Prater (mn367@cam.ac.uk), 2020                     
#   Centre for Trophobast Research, University of Cambridge
#   Script: DESeq2 (with shrinkage) Analysis of STAR bam files 
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
})

Project       <- "CTR_gjb2_0001_STAR_DESeq2_shrinkage_BlockSex"
significance  <- 0.05
l2fc          <- 1
col_1st       <- "firebrick2"
col_2nd       <- "steelblue3"

Base.dir      <- "/Users/malwinaprater/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/"
setwd(Base.dir)
Res.dir      <- "/Users/malwinaprater/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex"

message("+-------------------------------------------------------------------------------")
message("+                       Prepare sample table                                    ")
message("+-------------------------------------------------------------------------------")

sampleTable  <- read.csv("sampleTable.csv")
sample_id    <- sampleTable$sample
HTSeq.dir    <- paste(Base.dir,"/HTSEQ-COUNTS_CTR_gjb2_0001", sep="")
list.files(HTSeq.dir)
sampleTable$Sex <- c("male", "female", "female", "male","male",  "female", "female", "female", "male","male",   "female", "male","male","male")

message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)          
ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description', 'go_id', 'name_1006' ), mart = ensembl)    


message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ Sex + condition)
dds <- DESeq(ddsHTSeq)

setwd(Res.dir)

message("+-------------------------------------------------------------------------------")
message("+  Get the RESULTS    ")
message("+-------------------------------------------------------------------------------")

resultsNames(dds)

res <- lfcShrink(dds, coef="condition_Second_trimester_vs_First_trimester", type="normal")
res <- res[order(res$padj),]

message("+-------------------------------------------------------------------------------")
message("+ DESeq2 Results Tables")
message("+-------------------------------------------------------------------------------")

contrast_levels <- levels(dds$condition)
mcols(res, use.names = T)

results.df <- as.data.frame(res)
results.df$ensembl_gene_id <- rownames(results.df)
results.df <- merge(results.df, ensEMBL2id, by="ensembl_gene_id")
results.df$description <- gsub("..Source.*", "", results.df$description)
results.df <- results.df[order(results.df$padj),]

resSig.ann <- subset(results.df, padj < significance)
resSig.ann <- resSig.ann[order(resSig.ann$padj),]

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'ensembl_gene_id'
resdata_ann <- merge(resdata, ensEMBL2id[,c(1:3)], by.x= "ensembl_gene_id")

counts_norm <- counts(dds, normalized=TRUE)


message("+-------------------------------------------------------------------------------")
message("+                           add fold change values                              ")
message("+-------------------------------------------------------------------------------")

resSig.ann$FoldChange <-  2^resSig.ann$log2FoldChange
resSig.ann$FC <- 2^(abs(resSig.ann$log2FoldChange)) 
#change the sign of FC according to log2FC
resSig.ann$FC<-ifelse(resSig.ann$log2FoldChange<0,resSig.ann$FC*(-1),resSig.ann$FC*1) 
resSig.ann <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > 1)




message("+             Sex gened in our DEGs                          ")

sex_genes_Gong <- read.table("sex_genes_Gong2018.txt")


RESULTS_1 <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > 1 & resSig.ann$padj < 0.05)
length(unique(RESULTS_1$external_gene_name[RESULTS_1$external_gene_name %in% sex_genes_Gong$V2])) # 0
length(unique(RESULTS_1$external_gene_name)) # 3242
length(unique(sex_genes_Gong$V2)) # 101



message("+-------------------------------------------------------------------------------")
message("+ Now creating gene plots.")
message("+-------------------------------------------------------------------------------")

elementTextSize <- 8

makeGeneCountPlot <- function(gene2plot,outdir) {
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(dds, gene=gene2plot, intgroup=c("condition"), normalized=TRUE, returnData=TRUE)
  t2$samples    <- rownames(t2)
  pdf(paste(outdir, Project, "-DGE_", gene2plot, "_individual.pdf", sep=""),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ ggplot(t2, aes(x=samples, y=count, fill=condition)) + geom_bar(stat="identity", alpha=.75) + 
      scale_fill_manual(values = c("green4","skyblue2",  "blue4", "#5f022f", "violetred", "olivedrab", "cornflowerblue", "blue")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ggtitle(paste("CTR_mz205_0001 ::: ", gene2plot, " ::: ", genename2plot, sep="")) })
  dev.off()
  print(paste("Created plot for", gene2plot), sep=" ")
}

makeGeneCountPlot('ENSG00000087245')
makeGeneCountPlot('ENSG00000118113')
makeGeneCountPlot('ENSG00000113916')
makeGeneCountPlot('ENSG00000122565')



message("+-------------------------------------------------------------------------------")
message("+ Run transformations")
message("+-------------------------------------------------------------------------------")

rld <- DESeq2::rlogTransformation(dds, blind=TRUE)


message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

elementTextSize <- 14
pca = prcomp(t(assay(rld)))
rv = rowVars(assay(rld))

pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

scores <- data.frame(pca$x, sampleTable)
sampleTable$Sample_short <- sampleTable$Sample
sampleTable$Sample_short <- gsub("First_trimester","FT",sampleTable$Sample_short)
sampleTable$Sample_short <- gsub("Second_trimester","ST",sampleTable$Sample_short)

Supp_Fig_1_A <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sampleTable$condition)), shape = (factor(sampleTable$Sex)) )) + 
  geom_point(size = 5 ) + 
  xlab(pc1lab) + ylab(pc2lab) + 
  geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) + 
  scale_shape_manual(name="Sex", values = c(1,2)) +
  scale_colour_manual(name="Stage", values = c(col_1st, col_2nd)) +
  scale_fill_manual(name="Stage", values = c(col_1st, col_2nd)) +
  theme(text = element_text(size=elementTextSize)) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))+
  coord_fixed(ratio = 1, xlim = c(-60,60), ylim = c(-60,60), expand = TRUE, clip = "on")

Supp_Fig_1_A_labs <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sampleTable$condition)) )) + 
  geom_point(size = 5 ) + #geom_text_repel(aes(label=sampleTable$Sample_short), col = "black") +
  xlab(pc1lab) + ylab(pc2lab) + 
  geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) + 
  scale_fill_manual(name="Stage", values = c(col_1st, col_2nd)) +
  scale_colour_manual(name="Stage", values = c(col_1st, col_2nd)) +
  theme(text = element_text(size=elementTextSize)) + 
  coord_fixed(ratio = 1, xlim = c(-60,60), ylim = c(-60,60), expand = TRUE, clip = "on")




message("+-------------------------------------------------------------------------------")
message("+                                PCA explained                                  ")
message("+-------------------------------------------------------------------------------")

ensanno <- ensEMBL2id[,c(1:2)]
ensanno <- ensanno[!duplicated(ensanno),]
loadings                         <- as.data.frame(pca$rotation)
loadings$ensembl_gene_id         <- rownames(loadings)
loadings                         <- merge(loadings, ensanno, by="ensembl_gene_id")

pca.1         <-  loadings[ order(loadings$PC1,decreasing=TRUE), ]
pca.1.25      <-  pca.1[c(1:25),]
pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

pca.2         <-  loadings[ order(loadings$PC2,decreasing=TRUE), ]
pca.2.25      <-  pca.2[c(1:25),]
pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + geom_point(size = 5 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

Supp_Fig_1_B <- plot_grid(pca.1.25.plot, pca.2.25.plot, labels=c(" ", " "), ncol = 1, nrow = 2)



message("+-------------------------------------------------------------------------------")
message("+       Hierarchical clustering        ")
message("+-------------------------------------------------------------------------------")

sample_distances <- dist(t(assay(rld)))
Supp_Fig_1_C_hclust <- ggdendrogram(hclust(sample_distances), rotate = FALSE, segments = TRUE)
Supp_Fig_1_B_C  <-  plot_grid(Supp_Fig_1_B, Supp_Fig_1_C_hclust, labels=c("B", "C"), ncol = 2, nrow = 1, scale = c(0.95, 0.95), rel_widths = c(1, 1))



message("+-------------------------------------------------------------------------------")
message("+                                volcano plot                                   ")
message("+-------------------------------------------------------------------------------")

data <- data.frame(gene = results.df$ensembl_gene_id,
                   symbol = results.df$external_gene_name,
                   padj = results.df$padj,
                   pvalue = -log10(results.df$padj), 
                   lfc = results.df$log2FoldChange)
data <- na.omit(data)

data <- data %>%  dplyr::mutate(color = ifelse(data$lfc > 1 & data$pvalue > 1.3, yes = "Second_trimester", no = ifelse(data$lfc < -1 & data$pvalue > 1.3, yes = "First_trimester", no = "none")))
# pvalue = 2 for padj = 0.01
# pvalue = 1.3 for padj = 0.05

colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
  scale_x_continuous( limits = c(-9, 9)) + theme(legend.position = "none") + 
  #xlab(expression(log[2]("First_trimester" / "Second_trimester"))) + 
  xlab(expression( title = "log[2] NormCounts (First_trimester/Second_trimester)")) + 
  ylab(expression(-log[10]("adjusted p-value"))) +   ggtitle(label = "" , subtitle = "") +  
  geom_vline(xintercept = -1, colour = "black",linetype="dotted") +  
  geom_vline(xintercept = 1, colour = "black",linetype="dotted") + 
  geom_hline(yintercept = 1.3, colour = "black",linetype="dotted") + 
  annotate(geom = "text", label = "First trimester", x = -4, y = 55, size = 7, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 4.5, y = 55, size = 7, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) +  
  theme(text = element_text(size=elementTextSize)) 

Fig_1_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 3 & padj < 0.00000001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"))

pdf(paste("Fig_1b_volcano", Project,  "corrLabs.pdf", sep="_"), onefile=FALSE, width=7, height=7) 
par(bg=NA)
Fig_1_volcano
dev.off()


message("+-------------------------------------------------------------------------------")
message("+                             Starting GO/KEGG                                "  )
message("+-------------------------------------------------------------------------------")

RESULTS_1 <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > l2fc & resSig.ann$padj < significance)
resdata_simplified <- resdata_ann
resdata_simplified$mean_1st <- rowMeans(resdata_simplified[,c(8:15)])
resdata_simplified$mean_2nd <- rowMeans(resdata_simplified[,c(16:21)])
resdata_simplified <- resdata_simplified[,c(1,2,3,7,22:25)]

message("+-------------------------------------------------------------------------------")
message("+                                  enrichKegg                                   ")
message("+-------------------------------------------------------------------------------")

Kegg_genes <- na.omit(RESULTS_1)
head(Kegg_genes)
colnames(Kegg_genes)[colnames(Kegg_genes)=="entrezgene_id"] <- "Gene_ID" 

foldchanges = Kegg_genes$log2FoldChange
names(foldchanges) = Kegg_genes$Gene_ID
foldchanges <- sort(foldchanges, decreasing = T)
head(foldchanges)

kk_down <- enrichKEGG(names(foldchanges[foldchanges<0]), organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg", universe = unique(as.character(resdata_ann$entrezgene_id)))
head(summary(kk_down))
kk2_down <- setReadable(kk_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kk_res_down <- as.data.frame(kk2_down)
kk_res_down$direction <- "down"

kk_up <- enrichKEGG(names(foldchanges[foldchanges>0]), organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg", universe = unique(as.character(resdata_ann$entrezgene_id)))
head(summary(kk_up))
kk2_up <- setReadable(kk_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kk_res_up <- as.data.frame(kk2_up)
kk_res_up$direction <- "up"


kk_results <- rbind(kk_res_up, kk_res_down)

selected_kegg <- c("hsa04514" ,"hsa04060", "hsa04022",  "hsa04015", "hsa04141",  "hsa04921", "hsa04062", "hsa04350", "hsa04270", "hsa04010","hsa04911" ,"hsa03320" ,"hsa04918" ,"hsa00510",  "hsa04923", "hsa02010", "hsa04310", "hsa04141", "hsa04151", "hsa04072","hsa04020",  "hsa04512", "hsa04915", "hsa00140","hsa04928","hsa04612", "hsa04141",  "hsa05320", "hsa04371")   # "hsa04913","hsa04929",

enrichKegg_selected <- kk_results[kk_results$ID %in% selected_kegg,]
enrichKegg_selected$Description <- gsub("endoplasmic reticulum", "ER"  , enrichKegg_selected$Description)

list_up <- list()
list_down <- list()
for (i in 1:nrow(enrichKegg_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$external_gene_name %in% unlist(strsplit(enrichKegg_selected[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <-tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <-tmp_down
}
enrichKegg_selected$genes_UP <- as.numeric(list_up)
enrichKegg_selected$genes_DOWN <- as.numeric(list_down)
enrichKegg_selected$genes_DOWN <- -enrichKegg_selected$genes_DOWN
enrichKegg_molten <- melt(enrichKegg_selected[,c(1,2,6,7,11:12)], id.vars=c("ID","Description","p.adjust", "qvalue") )
enrichKegg_molten$Description <- gsub( "Parathyroid hormone synthesis, secretion and action", "Parathyroid hormone synth., secr., action", enrichKegg_molten$Description)

p_kegg_mlt <- ggplot(enrichKegg_molten, aes(x=reorder(Description, -qvalue), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue))) +
  coord_flip() +  xlab("KEGG Pathways") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +  ylab("Gene count") +
  ylim(-max(abs(enrichKegg_molten$value)), max(abs(enrichKegg_molten$value)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))

p_kegg_mlt







enrichr_biocart <- biocarta_enrichr_up
enrichr_biocart$Description <- gsub(" Homo sapiens h.*", ""  , enrichr_biocart$Term)
enrichr_biocart$Description <- gsub("Activation of Csk by cAMP-dependent Protein Kinase Inhibits Signaling through the T Cell Receptor", "Activ. of Csk by cAMP-dep. Prot. Kinase Inhibits Sig. through T Cell Receptor"  , enrichr_biocart$Description)

list_up <- list()
list_down <- list()
for (i in 1:nrow(enrichr_biocart)){
  df_tmp <- resdata_simplified[resdata_simplified$external_gene_name %in% unlist(strsplit(enrichr_biocart[i, "Genes"] , split=";")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <-tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <-tmp_down
}
enrichr_biocart$genes_UP <- as.numeric(list_up)
enrichr_biocart$genes_DOWN <- as.numeric(list_down)
enrichr_biocart$genes_DOWN <- -enrichr_biocart$genes_DOWN


p_biocarta <- ggplot(enrichr_biocart, aes(x=reorder(Description, -Adjusted.P.value), y=as.numeric(genes_UP))) + geom_bar(stat="identity", aes(alpha = Combined.Score) ,fill = col_2nd) +
  coord_flip() +  xlab("Biocarta 2016") +
  ylab("Gene count") +
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))

p_biocarta

pdf(paste(Project, "barplot_biocarta2016",  ".pdf", sep="_"), onefile=FALSE, width=10, height=4) 
par(bg=NA)
p_biocarta
dev.off()


message("+-------------------------------------------------------------------------------")
message("+                              Cluster profiler                                 ")
message("+-------------------------------------------------------------------------------")

bkcg_genes <- as.character(unique(results.df$entrezgene_id))
geneListx <- unique(RESULTS_1[,c(9,3)])
geneListx <- geneListx[!is.na(geneListx$entrezgene),]
geneList <- geneListx$log2FoldChange
names(geneList) <- geneListx$entrezgene
geneList <- sort(geneList, decreasing = T )
geneList <- geneList[unique(names(geneList))]
gene <- names(geneList)[abs(geneList) > 1]

# GO 
gsego_cc <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "CC", nPerm= 1000, minGSSize  = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsego_bp <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "BP", nPerm= 1000, minGSSize  = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsego_mf <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "MF", nPerm= 1000, minGSSize  = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsecc_df <- as.data.frame(gsego_cc)
gsebp_df <- as.data.frame(gsego_bp)
gsemf_df <- as.data.frame(gsego_mf)

# GO over-representation test
ego_bp <- enrichGO(gene = gene, universe = bkcg_genes, OrgDb = org.Hs.eg.db, ont = "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
ego_mf <- enrichGO(gene = gene, universe = bkcg_genes, OrgDb = org.Hs.eg.db, ont = "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
ego_cc <- enrichGO(gene = gene, universe = bkcg_genes, OrgDb = org.Hs.eg.db, ont = "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)

ego_bp_df <- as.data.frame(ego_bp)
ego_cc_df <- as.data.frame(ego_cc) 
ego_mf_df <- as.data.frame(ego_mf)
ego_bp_df$category <- "BP"
ego_mf_df$category <- "MF"
ego_cc_df$category <- "CC"

ego <- rbind(ego_bp_df, ego_cc_df, ego_mf_df )
head(ego)


message("+-------------------------------------------------------------------------------")
message("+                              Plot gse   results                               ")
message("+-------------------------------------------------------------------------------")

gsebp_selected_terms <-  c("ER to Golgi vesicle-mediated transport", "response to endoplasmic reticulum stress", "cytokine-mediated signaling pathway", "adaptive immune response", "regulation of vasculature development", "cellular response to oxygen-containing compound",  "G protein-coupled receptor signaling pathway", "positive regulation of ERK1 and ERK2 cascade",  "carboxylic acid transport", "activation of MAPK activity", "cellular response to lipid", "calcium ion transport",  "positive regulation of peptide secretion", "regulation of transmembrane transport", "Golgi vesicle transport",  "regulation of GTPase activity", "positive regulation of cell adhesion", "regulation of cytosolic calcium ion concentration", "response to oxidative stress", "immune response-activating signal transduction",  "extracellular matrix organization", "IRE1-mediated unfolded protein response", "endoplasmic reticulum unfolded protein response", "regulation of complement activation", "negative regulation of response to endoplasmic reticulum stress", "cell redox homeostasis", 'ribosome biogenesis', "regulation of G protein-coupled receptor signaling pathway","protein folding in endoplasmic reticulum", 'regulation of protein secretion', "lipid transport","reactive oxygen species metabolic process") 

ego_BP_selected <- ego_bp_df
ego_BP_selected <- ego_bp_df[ego_bp_df$Description %in% gsebp_selected_terms,]
ego_BP_selected$Description <- gsub("G protein-coupled receptor", "GPCR" , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("regulation", "reg." , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("oxygen", "O2" , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("nucleotide", "nt" , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("GPCR signaling pathway, coupled to cyclic nt second messenger", "GPCR signal. pathway*" , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("signaling pathway", "signaling" , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("growth factor stimulus", "growth factor" , ego_BP_selected$Description)
ego_BP_selected$Description <- gsub("endoplasmic reticulum", "ER" , ego_BP_selected$Description)

list_up <- list()
list_down <- list()
for (i in 1:nrow(ego_BP_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$external_gene_name %in% unlist(strsplit(ego_BP_selected[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <- tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <- tmp_down
}
ego_BP_selected$genes_UP <- as.numeric(list_up)
ego_BP_selected$genes_DOWN <- -as.numeric(list_down)



ego_BP_molten <- melt(ego_BP_selected[,c(2,6,7,8,9:11)], id.vars=c("Description","p.adjust","qvalue","geneID","Count") )


p_bp_mlt <- ggplot(ego_BP_molten, aes(x=reorder(Description, -qvalue), y=value,fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue)))+
  coord_flip()+
  xlab("GO: Biological Processes") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +
  ylab("Gene count") +
  ylim(-max(abs(ego_BP_molten$value)), max(abs(ego_BP_molten$value)))+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

#pdf(paste(Project, "ego_BP_selected", "barplot20_padj0.05_l2fc1_MLT", ".pdf", sep="_"), onefile=FALSE, width=20, height=60) 
#par(bg=NA)
#p_bp_mlt 
#dev.off()





gsemf_selected_terms <-  c("molecular transducer activity", "signaling receptor activity", "transmembrane signaling receptor activity" , "calcium ion binding" , "extracellular matrix structural constituent", "receptor ligand activity","substrate-specific channel activity" , "receptor regulator activity" , "passive transmembrane transporter activity", "cation channel activity","G protein-coupled receptor activity", "cytokine activity", "G protein-coupled receptor binding", "apolipoprotein binding","calcium-release channel activity","metal ion transmembrane transporter activity","growth factor activity","glycosaminoglycan binding")


gse_MF_selected <- ego_mf_df[ego_mf_df$Description %in% gsemf_selected_terms,]
gse_MF_selected$Description <- gsub("transmembrane", "transmem." , gse_MF_selected$Description)
gse_MF_selected$Description <- gsub("G protein-coupled receptor", "GPCR" , gse_MF_selected$Description)

list_up <- list()
list_down <- list()
for (i in 1:nrow(gse_MF_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$external_gene_name %in% unlist(strsplit(gse_MF_selected[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <- tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <- tmp_down
}
gse_MF_selected$genes_UP <- as.numeric(list_up)
gse_MF_selected$genes_DOWN <- -as.numeric(list_down)


gse_MF_selected <- melt(gse_MF_selected[,c(2,6,7,8,9:11)], id.vars=c("Description","p.adjust","qvalue","geneID","Count") )

p_MF_mlt <- ggplot(gse_MF_selected, aes(x=reorder(Description, -qvalue), y=value,fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue)))+
  coord_flip()+
  xlab("GO: Molecular Function") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +
  ylab("Gene count") +
  ylim(-max(abs(gse_MF_selected$value)), max(abs(gse_MF_selected$value)))+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 





gsecc_selected_terms <- c("secretory vesicle", "cell surface" ,"endocytic vesicle membrane" , "transmembrane transporter complex", "external side of plasma membrane", "ion channel complex" ,
                          "secretory granule membrane", "transporter complex", "collagen-containing extracellular matrix" ,"plasma membrane protein complex", "ribosomal subunit", "catalytic step 2 spliceosome" ,"endoplasmic reticulum chaperone complex", "endoplasmic reticulum lumen" ,"MHC class II protein complex", "ER to Golgi transport vesicle membrane","cation channel complex")     

gse_CC_selected <- ego_cc_df[ego_cc_df$Description %in% gsecc_selected_terms,]
gse_CC_selected$Description <- gsub("extracellular matrix", "ECM" , gse_CC_selected$Description)
gse_CC_selected$Description <- gsub("endoplasmic reticulum", "ER" , gse_CC_selected$Description)


list_up <- list()
list_down <- list()
for (i in 1:nrow(gse_CC_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$external_gene_name %in% unlist(strsplit(gse_CC_selected[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <- tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <- tmp_down
}
gse_CC_selected$genes_UP <- as.numeric(list_up)
gse_CC_selected$genes_DOWN <- -as.numeric(list_down)

gse_CC_selected <- melt(gse_CC_selected[,c(2,6,7,8,9:11)], id.vars=c("Description","p.adjust","qvalue","geneID","Count") )


p_CC_mlt <- ggplot(gse_CC_selected, aes(x=reorder(Description, -qvalue), y=value,fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue)))+
  coord_flip()+
  xlab("GO: Cellular Component") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +
  ylab("Gene count") +
  ylim(-max(abs(gse_CC_selected$value)), max(abs(gse_CC_selected$value)))+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



plot_GO_Kegg_2sided <- plot_grid(p_kegg_mlt, p_bp_mlt, p_MF_mlt, p_CC_mlt,  labels=c("", "", "", ""), align="hv", ncol = 2, nrow = 2, rel_heights = c(1, 0.7))

pdf(paste(Project, "plot_GO_Kegg_sided", ".pdf", sep="_"), width=15, height=9, onefile=FALSE)
par(bg=NA)
plot_GO_Kegg_2sided 
dev.off()





message("+-------------------------------------------------------------------------------")
message("+                             Transport genes                                 "  )
message("+-------------------------------------------------------------------------------")

ensEMBL2id_go_transport <-  ensEMBL2id_go[grepl("transport", ensEMBL2id_go$name_1006),]
ensEMBL2id_go_transport$calcium_transport <- grepl("calcium", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport$transmembrane_transport <- grepl("transmembrane transport", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport$lipid_transport <- grepl("lipid", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport$hormone_transport <- grepl("hormone", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport$secretion_transport <- grepl("secretion", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport$drug_transport <- grepl("drug", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport$oxygen <- grepl("oxygen", ensEMBL2id_go_transport$name_1006)
ensEMBL2id_go_transport_1 <- subset(ensEMBL2id_go_transport, ensEMBL2id_go_transport$lipid_transport == TRUE | ensEMBL2id_go_transport$calcium_transport == TRUE |ensEMBL2id_go_transport$transmembrane_transport == TRUE |ensEMBL2id_go_transport$hormone_transport == TRUE |ensEMBL2id_go_transport$secretion_transport == TRUE |ensEMBL2id_go_transport$drug_transport == TRUE |ensEMBL2id_go_transport$oxygen == TRUE )

results_transport <-  na.omit(unique(results.df[results.df$ensembl_gene_id %in% ensEMBL2id_go_transport_1$ensembl_gene_id, ])) 
results_transport <- results_transport[order(results_transport$log2FoldChange),]

data <- data.frame(gene = results_transport$ensembl_gene_id,
                   symbol = results_transport$external_gene_name,
                   padj = results_transport$padj,
                   pvalue = -log10(results_transport$padj), 
                   lfc = results_transport$log2FoldChange)
data <- na.omit(data)

data <- data %>%  dplyr::mutate(color = ifelse(data$lfc > 1 & data$pvalue > 1.3, yes = "Second_trimester", no = ifelse(data$lfc < -1 & data$pvalue > 1.3, yes = "First_trimester", no = "none")))

colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
  scale_x_continuous( limits = c(-5.5, 5.5)) +  theme(legend.position = "none") + 
  ggtitle(label = "" , subtitle = "") +  
  xlab(expression(log[2]("First_trimester" / "Second_trimester"))) +  
  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_vline(xintercept = -1, colour = "black",linetype="dotted") +  
  geom_vline(xintercept = 1, colour = "black",linetype="dotted") + 
  geom_vline(xintercept = 0, colour = "black") +  geom_hline(yintercept = 1.3, colour = "black") + 
  annotate(geom = "text", label = "First trimester", x = -4, y = 33, size = 6, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 3, y = 33, size = 6, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) + 
  theme(text = element_text(size=elementTextSize)) 

Fig_4_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 2 & padj < 0.00000001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.iter = 10, force = TRUE)

pdf(paste("Fig_4_", "volcano_plot_TRANSPORT", ".pdf", sep=""), width=7, height=7, onefile=FALSE)
par(bg=NA)
Fig_4_volcano
dev.off()



ensEMBL2id_go_lipid_transport <- subset(ensEMBL2id_go_transport, ensEMBL2id_go_transport$lipid_transport == TRUE)
results_lipid_transport <-  na.omit(unique(results.df[results.df$ensembl_gene_id %in% ensEMBL2id_go_lipid_transport$ensembl_gene_id, ])) 
results_lipid_transport <- results_lipid_transport[order(results_lipid_transport$log2FoldChange),]
saveRDS(results_lipid_transport$external_gene_name, "Lipid_Transport_GO_genes.rds")



message("+-------------------------------------------------------------------------------")
message("+              Proliferation and differentiation  genes                       "  )
message("+-------------------------------------------------------------------------------")

ensEMBL2id_go_proliferation  <- ensEMBL2id_go
ensEMBL2id_go_proliferation$prolif <- grepl("proliferation", ensEMBL2id_go_proliferation$name_1006)
ensEMBL2id_go_proliferation$diff <- grepl("differentiation", ensEMBL2id_go_proliferation$name_1006)
ensEMBL2id_go_diff <- subset(ensEMBL2id_go_proliferation,  ensEMBL2id_go_proliferation$diff == TRUE)
ensEMBL2id_go_prolif <- subset(ensEMBL2id_go_proliferation, ensEMBL2id_go_proliferation$prolif == TRUE)
ensEMBL2id_go_prolif_diff <- subset(ensEMBL2id_go_proliferation, ensEMBL2id_go_proliferation$prolif == TRUE | ensEMBL2id_go_proliferation$diff == TRUE)

results_prolif_diff <-  unique(results.df[results.df$ensembl_gene_id %in% ensEMBL2id_go_prolif_diff$ensembl_gene_id, ]) 
results_prolif_diff <- na.omit(results_prolif_diff)
results_prolif_diff <- results_prolif_diff[order(results_prolif_diff$log2FoldChange),]
saveRDS(results_prolif_diff$external_gene_name, "Prolif_differentiation_GO_genes.rds")

data <- data.frame(gene = results_prolif_diff$ensembl_gene_id,
                   symbol = results_prolif_diff$external_gene_name,
                   padj = results_prolif_diff$padj,
                   pvalue = -log10(results_prolif_diff$padj), 
                   lfc = results_prolif_diff$log2FoldChange)
data <- na.omit(data)

data <- data %>%  dplyr::mutate(color = ifelse(data$lfc > 1 & data$pvalue > 1.3, yes = "Second_trimester", no = ifelse(data$lfc < -1 & data$pvalue > 1.3, yes = "First_trimester", no = "none")))

colored <- ggplot(data, aes(x = lfc, y = pvalue)) + 
  geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
  scale_x_continuous( limits = c(-5, 5)) +  theme(legend.position = "none") +  
  ggtitle(label = "" , subtitle = "") +  
  xlab(expression(log[2]("First_trimester" / "Second_trimester"))) + ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_vline(xintercept = 0, colour = "black") +  geom_hline(yintercept = 1.3, colour = "black") + 
  annotate(geom = "text", label = "First trimester", x = -3, y = 30, size = 6, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 3, y = 30, size = 6, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) +  
  theme(text = element_text(size=elementTextSize)) 

Fig_5_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 2.2 & padj < 0.00000001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.iter = 10, force = T)



message("+-------------------------------------------------------------------------------")
message("+                 Genes associated with hormonal activity                       ")
message("+-------------------------------------------------------------------------------")

ensEMBL2id_go_hormone <- ensEMBL2id_go
ensEMBL2id_go_hormone$hormone <- grepl("hormone", ensEMBL2id_go_hormone$name_1006)
ensEMBL2id_go_hormone <- subset(ensEMBL2id_go_hormone, ensEMBL2id_go_hormone$hormone == TRUE)
ensEMBL2id_go_hormone$hormone <- grepl("activity", ensEMBL2id_go_hormone$name_1006)
ensEMBL2id_go_hormone <- subset(ensEMBL2id_go_hormone, ensEMBL2id_go_hormone$hormone == TRUE)

hormone_genes_table <- unique(results.df[results.df$ensembl_gene_id %in% ensEMBL2id_go_hormone$ensembl_gene_id,])
hormone_genes_table <- hormone_genes_table[!duplicated(hormone_genes_table$external_gene_name) | !duplicated(hormone_genes_table$ensembl_gene_id),]
hormone_genes_table <- na.omit(hormone_genes_table)
hormone_genes_table <- hormone_genes_table[order(hormone_genes_table$log2FoldChange),]
saveRDS(hormone_genes_table$external_gene_name, "hormone_activity_GO_genes.rds")

data <- data.frame(gene = hormone_genes_table$ensembl_gene_id,
                   symbol = hormone_genes_table$external_gene_name,
                   padj = hormone_genes_table$padj,
                   pvalue = -log10(hormone_genes_table$padj), 
                   lfc = hormone_genes_table$log2FoldChange)
data <- na.omit(data)

data <- data %>%  dplyr::mutate(color = ifelse(data$lfc > 1 & data$pvalue > 1.3, yes = "Second_trimester", no = ifelse(data$lfc < -1 & data$pvalue > 1.3, yes = "First_trimester", no = "none")))

colored <- ggplot(data, aes(x = lfc, y = pvalue)) +   geom_point(aes(color = factor(color)), size = 1.75, alpha = 0.8, na.rm = T) + 
  scale_x_continuous( limits = c(-4, 4)) +  theme(legend.position = "none") +  ggtitle(label = "" , subtitle = "") +  
  xlab(expression(log[2]("First_trimester" / "Second_trimester"))) +  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_vline(xintercept = 0, colour = "black") + geom_hline(yintercept = 1.3, colour = "black") + 
  annotate(geom = "text", label = "First trimester", x = -3, y = 15, size = 6, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 3, y = 15, size = 6, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) +  
  theme(text = element_text(size=elementTextSize)) 

Fig_3_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 1 & padj < 0.001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.iter = 10)



message("+-------------------------------------------------------------------------------")
message("+              Fig.2a      ER Genes                                             ")
message("+-------------------------------------------------------------------------------")

kk_genes <- strsplit(kk_results[kk_results$ID == "hsa04141", ]$geneID, split = "/")
ensEMBL2id_go_ER     <- ensEMBL2id_go
ensEMBL2id_go_ER$ER1 <- grepl("endoplasmic reticulum", ensEMBL2id_go_ER$name_1006)
ensEMBL2id_go_ER$ER2 <- grepl("endoplasmic reticulum", ensEMBL2id_go_ER$description)
ensEMBL2id_go_ER     <- subset(ensEMBL2id_go_ER, ensEMBL2id_go_ER$ER1 == TRUE | ensEMBL2id_go_ER$ER2 == TRUE)
ensEMBL2id_go_ER$ER1 <- grepl("endoplasmic reticulum stress", ensEMBL2id_go_ER$name_1006)
ensEMBL2id_go_ER$ER2 <- grepl("transport", ensEMBL2id_go_ER$name_1006)
ensEMBL2id_go_ER$ER3 <- grepl("endoplasmic reticulum to Golgi vesicle-mediated", ensEMBL2id_go_ER$name_1006)
ensEMBL2id_go_ER     <- subset(ensEMBL2id_go_ER, ensEMBL2id_go_ER$ER1 == TRUE |  ensEMBL2id_go_ER$ER3 == TRUE  | ensEMBL2id_go_ER$ER2 == TRUE )
ER_assoc_GO_genes    <- unique(ensEMBL2id_go_ER$external_gene_name)
ER_genes <- unique(unlist((append(as.vector(kk_genes), as.vector(ER_assoc_GO_genes)))))
selected_genes <- resSig.ann[resSig.ann$external_gene_name %in% ER_genes,]
selected_genes <- subset(selected_genes, abs(selected_genes$log2FoldChange) > 0.6)


message("+-------------------------------------------------------------------------------")
message("+          Fig. 2b     secretory Genes                                          ")
message("+-------------------------------------------------------------------------------")

ego_terms_go_secretion <- ensEMBL2id_go
ego_terms_go_secretion$secretion <- grepl("secret", ego_terms_go_secretion$name_1006)
ego_terms_go_secretion <- subset(ego_terms_go_secretion, ego_terms_go_secretion$secretion == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% ego_terms_go_secretion$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)


message("+-------------------------------------------------------------------------------")
message("+             Fig. 2c  oxy genes                                                ")
message("+-------------------------------------------------------------------------------")

ego_terms_go_Oxy <- ego
ego_terms_go_Oxy$Oxygen <- grepl("oxygen", ego_terms_go_Oxy$Description)
ego_terms_go_Oxy$Oxygen2 <- grepl("oxid", ego_terms_go_Oxy$Description)
ego_terms_go_Oxy <- subset(ego_terms_go_Oxy, ego_terms_go_Oxy$Oxygen == TRUE | ego_terms_go_Oxy$Oxygen2 == TRUE )
oxy_genes <- as.character(unlist(ego_terms_go_Oxy$geneID))
oxy_genes    <- strsplit(oxy_genes, ("/"))
oxy_genes <- unique(unlist(oxy_genes))
selected_genes <- resSig.ann[resSig.ann$external_gene_name %in% oxy_genes,]
selected_genes <- subset(selected_genes, abs(selected_genes$log2FoldChange) > 0.6)


message("+-------------------------------------------------------------------------------")
message("+             Fig. TFs :                                            ")
message("+-------------------------------------------------------------------------------")

TF_db <- read.csv("/Users/malwina/Documents/CTR-Groups/Graham_Burton/2019_11_07__tfcheckpoint_database.csv")
TF_genes <- resSig.ann[resSig.ann$external_gene_name %in% TF_db$Gene_Name,]
TF_genes <- subset(TF_genes, abs(TF_genes$log2FoldChange) > 0.6)
as.character(TF_genes$external_gene_name)


message("+-------------------------------------------------------------------------------")
message("+              Fig.4e- ECM related genes                                            ")
message("+-------------------------------------------------------------------------------")

Hsa04974_genes <- strsplit(kk_results[kk_results$ID == "hsa04974",]$geneID, split = "/")  # Protein digestion and absorption
ensEMBL2id_go_ECM <- ensEMBL2id_go
ensEMBL2id_go_ECM <- subset(ensEMBL2id_go_ECM, ensEMBL2id_go_ECM$go_id == "GO:0022617" | ensEMBL2id_go_ECM$go_id == "GO:0030198" |ensEMBL2id_go_ECM$go_id == "GO:0017090"  )
ECM_genes <- unique(ensEMBL2id_go_ECM$external_gene_name)
ECM_genes_all <- c(ECM_genes, unlist(Hsa04974_genes))
selected_genes <- resSig.ann[resSig.ann$external_gene_name %in% ECM_genes_all,]


message("+-------------------------------------------------------------------------------")
message("+             Fig. 5b  WNT genes                                                ")
message("+-------------------------------------------------------------------------------")

ego_terms_go_wnt <- ensEMBL2id_go
ego_terms_go_wnt$wnt <- grepl("Wnt", ego_terms_go_wnt$name_1006)
ego_terms_go_wnt <- subset(ego_terms_go_wnt, ego_terms_go_wnt$wnt == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% ego_terms_go_wnt$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)

message("+-------------------------------------------------------------------------------")
message("+             Fig. 5b  CC genes                                                ")
message("+-------------------------------------------------------------------------------")

cellcycle_genes <- ensEMBL2id_go
cellcycle_genes$cc <- grepl("cell cycle", cellcycle_genes$name_1006)
cellcycle_genes$cc2 <- grepl("cyclin", cellcycle_genes$name_1006)
cellcycle_genes <- subset(cellcycle_genes, cellcycle_genes$cc == TRUE | cellcycle_genes$cc2 == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% cellcycle_genes$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)


message("+-------------------------------------------------------------------------------")
message("+             Fig. 5b  Glycolytic genes                                                ")
message("+-------------------------------------------------------------------------------")

glycolytic_genes <- ensEMBL2id_go
glycolytic_genes$cc <- grepl("glycolytic", glycolytic_genes$name_1006)
glycolytic_genes$cc2 <- grepl("glycolysis", glycolytic_genes$name_1006)
glycolytic_genes <- subset(glycolytic_genes, glycolytic_genes$cc == TRUE | glycolytic_genes$cc2 == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% glycolytic_genes$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)






message("+-------------------------------------------------------------------------------")
message("+           Run on server::::::::                                               ")
message("+-------------------------------------------------------------------------------")


library(GENIE3)
set.seed(123) # For reproducibility of results

Base.dir <- "/storage/CTR-Projects/CTR_gjb2/CTR_gjb2_0001/BACK_UP_2020/CTR_gjb2_0001/bulk_vs_scRNA_seq/"
Base.dir <- "/storage/CTR-Projects/CTR_gjb2/CTR_gjb2_0001/BACK_UP_2020/CTR_gjb2_0001/bulk_vs_scRNA_seq/"
setwd(Base.dir)

#resdata_ann <- read.csv("/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex/CTR_gjb2_0001_STAR_DESeq2_shr_BlockSex_deseq2_resdata_ann.csv")
resdata_ann <- read.csv("/storage/CTR-Projects/CTR_gjb2/CTR_gjb2_0001/BACK_UP_2020/CTR_gjb2_0001/bulk_vs_scRNA_seq/CTR_gjb2_0001_STAR_DESeq2_shr_BlockSex_deseq2_resdata_ann.csv")
exprMatr <- unique(resdata_ann[,c(23,9:22)])

exprMatr[1:5,1:5]
idx <- rowSums(exprMatr[,-1]) > 10
exprMatr <- exprMatr[idx,]

exprMatr <- exprMatr[order(rowSums(exprMatr[,-1]), decreasing = TRUE),]
exprMatr <- exprMatr[!duplicated(exprMatr$external_gene_name),]
rownames(exprMatr) <- exprMatr$external_gene_name
exprMatr <- exprMatr[,-1]


DEGs <- unique(resdata_ann[resdata_ann$padj < 0.05 & abs(resdata_ann$log2FoldChange) > 0,]$external_gene_name)
exprMatr_DEG <- exprMatr[rownames(exprMatr) %in%  DEGs,]



#TF_db <- read.csv("/Users/malwina/Documents/CTR-Groups/Graham_Burton/2019_11_07__tfcheckpoint_database.csv")
#TF_genes <- as.character(TF_db[TF_db$Gene_Name %in%  DEGs, ]$Gene_Name)
TF_genes <- unique(motifAnnotations_hgnc$TF)[unique(motifAnnotations_hgnc$TF) %in%  DEGs ]

weightMat <- GENIE3(as.matrix(exprMatr_DEG), regulators =  as.character(rownames(exprMatr_DEG)), nCores = 1, verbose = TRUE )

# Tree method: RF
# K: sqrt
# Number of trees: 1000


dim(weightMat)
weightMat[1:5,1:5]

# Get the list of the regulatory links
#You can obtain the list of all the regulatory links (from most likely to least likely) with this command:
linkList <- getLinkList(weightMat, threshold=0.05) # reportMax=5 , threshold=0.1
edge_listsi <- linkList[!duplicated(linkList),]

dim(linkList)
 
head(linkList)
##   regulatoryGene targetGene    weight
#1         HOXC13    RPL3P12 0.1152300
#2         ZNF423 AC103810.2 0.1088218
#3         HIVEP3      KCNN2 0.1069245
#4         MLXIPL      SGSM1 0.1051512
#5         HOXC13 AC010240.1 0.1032896

linkList_0.07 <- linkList[linkList$weight > 0.07,]
linkList_0.08 <- linkList[linkList$weight > 0.08,]
linkList_0.05 <- linkList[linkList$weight > 0.05,]

weightMat[1:5, 1:5]

dim(weightMat)
weightMat_for_ht <- weightMat[rownames(weightMat) %in% linkList_0.08$regulatoryGene,colnames(weightMat) %in% linkList_0.08$targetGene]
dim(weightMat_for_ht)


weightMat_for_ht[weightMat_for_ht == 0] <- NA
ComplexHeatmap::Heatmap(weightMat_for_ht,na_col = "grey15")













library(RcisTarget)
library("AUCell")
library(SCENIC)

#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
#             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")

#for(featherURL in dbFiles)
#{
#  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
#}


# Select motif database to use (i.e. organism and distance around TSS)
data(motifAnnotations_hgnc)
motifRankingsTSS <- importRankings("hg19-tss-centered-10kb-7species.mc9nr.feather")
motifRankings <- importRankings("hg19-500bp-upstream-7species.mc9nr.feather")


# Load gene sets to analyze. e.g.:
#geneList1 <- as.character(DEGs)
geneList1 <- as.character(RESULTS_1$external_gene_name)
geneList1 <- geneList1[geneList1 %in% colnames(motifRankings@rankings)]
geneLists <- list(geneListName=geneList1)
# geneLists <- GSEABase::GeneSet(genes, setName="geneListName") # alternative


# Motif enrichment analysis for DEGs (l2fc>1):
motifEnrichmentTable_wGenes1 <- cisTarget(geneLists, motifRankings,
                                         motifAnnot=motifAnnotations_hgnc)
motifEnrichmentTable_wGenes1TSS <- cisTarget(geneLists, motifRankingsTSS,
                                         motifAnnot=motifAnnotations_hgnc)


motifEnrich_SPI1 <- unlist(strsplit( motifEnrichmentTable_wGenes[motifEnrichmentTable_wGenes$motif == "hocomoco__SPI1_HUMAN.H11MO.0.A",]$enrichedGenes, split = ";" ))
motifEnrich_CBX3 <- unlist(strsplit( motifEnrichmentTable_wGenes[motifEnrichmentTable_wGenes$motif == "dbcorrdb__CBX3__ENCSR000BRT_1__m9",]$enrichedGenes, split = ";" ))
motifEnrich_SPIB <- unlist(strsplit( motifEnrichmentTable_wGenes[motifEnrichmentTable_wGenes$motif == "hocomoco__SPIB_HUMAN.H11MO.0.A",]$enrichedGenes, split = ";" ))

motifEnrich_lc_ERG_FLI1 <- unlist(strsplit( motifEnrichmentTable_wGenes[motifEnrichmentTable_wGenes$motif == "tfdimers__MD00145",]$enrichedGenes, split = ";" ))

 

motifEnrich_CBX3_TSS <- unlist(strsplit( motifEnrichmentTable_wGenesTSS[motifEnrichmentTable_wGenesTSS$motif == "dbcorrdb__CBX3__ENCSR000BRT_1__m9",]$enrichedGenes, split = ";" ))

motifEnrich_lc_BCL6_tss <- unlist(strsplit( motifEnrichmentTable_wGenesTSS[motifEnrichmentTable_wGenesTSS$motif == "transfac_pro__M09104",]$enrichedGenes, split = ";" ))


motifEnrich_genes <- motifEnrich_CBX3

#motifEnrich_genes <- DEG_CBX3$external_gene_name[DEG_CBX3$external_gene_name %in% DEG_BCL6_tss$external_gene_name] # 75
#motifEnrich_genes <- DEG_CBX3tss$external_gene_name[DEG_CBX3tss$external_gene_name %in% DEG_BCL6_tss$external_gene_name] # 240
motifEnrich_genes <- c("PTGFR","ASTN1","ITGA8","NELL1","CLDN10","SLC24A4","HS3ST3B1","ABCC3","CIDEA","KIRREL2","LY75","GAD1","EREG","LEF1","TLL1","FAT1","GDNF","TRIM36","PPP2R2B","ADRB2","SEMA3A","WNT2") # top CBX3 in DEGs & has DMR


enrichR_res <- enrichr(as.character(motifEnrich_genes[motifEnrich_genes %in% resdata_ann[abs(resdata_ann$log2FoldChange) > 1,]$external_gene_name]), databases = c("KEGG_2019_Human", "WikiPathways_2019_Human","Reactome_2016",  "TF_Perturbations_Followed_by_Expression", "TRRUST_Transcription_Factors_2019","Chromosome_Location", "Epigenomics_Roadmap_HM_ChIP-seq", "BioCarta_2016","InterPro_Domains_2019", "Pfam_Domains_2019","Elsevier_Pathway_Collection","GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018","Enrichr_Submissions_TF-Gene_Coocurrence"))


DEG_CBX3 <- resSig.ann[resSig.ann$external_gene_name %in% motifEnrich_CBX3,]
DEG_CBX3 <- DEG_CBX3[abs(DEG_CBX3$log2FoldChange) > 1,]
nrow(DEG_CBX3[(DEG_CBX3$log2FoldChange) > 1,]) # 171  --- 522
nrow(DEG_CBX3[(DEG_CBX3$log2FoldChange) < -1,]) # 66  --- 175
nrow(resSig.ann[(resSig.ann$log2FoldChange) > 1,]) #  2123
nrow(resSig.ann[(resSig.ann$log2FoldChange) < -1,]) #  1137
nrow(DEG_CBX3) # 133
length(unique(resSig.ann[abs(resSig.ann$log2FoldChange) >1,]$external_gene_name))

DEG_CBX3tss <- resSig.ann[resSig.ann$external_gene_name %in% motifEnrich_CBX3_TSS,]
DEG_CBX3tss <- DEG_CBX3tss[abs(DEG_CBX3tss$log2FoldChange) > 1,]
nrow(DEG_CBX3tss[(DEG_CBX3tss$log2FoldChange) > 1,]) #  522
nrow(DEG_CBX3tss[(DEG_CBX3tss$log2FoldChange) < -1,]) #  175

DEG_ERG_FLI1 <- resSig.ann[resSig.ann$external_gene_name %in% motifEnrich_lc_ERG_FLI1,]
DEG_ERG_FLI1 <- DEG_ERG_FLI1[abs(DEG_ERG_FLI1$log2FoldChange) > 1,]
nrow(DEG_ERG_FLI1[(DEG_ERG_FLI1$log2FoldChange) > 1,]) # 99
nrow(DEG_ERG_FLI1[(DEG_ERG_FLI1$log2FoldChange) < -1,]) # 34
nrow(DEG_ERG_FLI1) # 133

DEG_BCL6_tss <- resSig.ann[resSig.ann$external_gene_name %in% motifEnrich_lc_BCL6_tss,]
DEG_BCL6_tss <- DEG_BCL6_tss[abs(DEG_BCL6_tss$log2FoldChange) > 1,]
nrow(DEG_BCL6_tss[(DEG_BCL6_tss$log2FoldChange) > 1,]) # 264
nrow(DEG_BCL6_tss[(DEG_BCL6_tss$log2FoldChange) < -1,]) # 82
nrow(DEG_BCL6_tss) # 346

length(DEG_CBX3[DEG_CBX3$external_gene_name %in% DEG_BCL6_tss$external_gene_name,]$external_gene_name) # 75/237 
length(DEG_CBX3tss[DEG_CBX3tss$external_gene_name %in% DEG_BCL6_tss$external_gene_name,]$external_gene_name) # 240 / 697



enrichR_res <- enrichr(motifEnrich_genes, databases = c("KEGG_2019_Human", "WikiPathways_2019_Human","Reactome_2016",  "TF_Perturbations_Followed_by_Expression", "TRRUST_Transcription_Factors_2019","Chromosome_Location", "Epigenomics_Roadmap_HM_ChIP-seq", "BioCarta_2016","InterPro_Domains_2019", "Pfam_Domains_2019","Elsevier_Pathway_Collection","GO_Biological_Process_2018", "GO_Cellular_Component_2018", "GO_Molecular_Function_2018","Enrichr_Submissions_TF-Gene_Coocurrence"))



enrichr_Wiki <- enrichR_res$WikiPathways_2019_Human
enrichr_Wiki <- enrichr_Wiki[enrichr_Wiki$Adjusted.P.value < 0.05,]
# Wnt Signaling WP428

enrichr_React <- enrichR_res$Reactome_2016
enrichr_React <- enrichr_React[enrichr_React$Adjusted.P.value < 0.05,]

enrichr_Kegg <- enrichR_res$KEGG_2019_Human
enrichr_Kegg <- enrichr_Kegg[enrichr_Kegg$Adjusted.P.value < 0.05,]
# Wnt signaling pathway
#
# BCL6 :: Vascular smooth muscle contraction

enrichr_GOBP <- enrichR_res$GO_Biological_Process_2018
enrichr_GOBP <- enrichr_GOBP[enrichr_GOBP$Adjusted.P.value < 0.05,]

enrichr_GOMF <- enrichR_res$GO_Molecular_Function_2018
enrichr_GOMF <- enrichr_GOMF[enrichr_GOMF$Adjusted.P.value < 0.05,]

enrichr_GOCC <- enrichR_res$GO_Cellular_Component_2018
enrichr_GOCC <- enrichr_GOCC[enrichr_GOCC$Adjusted.P.value < 0.05,]

enrichr_BioCarta <- enrichR_res$BioCarta_2016
enrichr_BioCarta <- enrichr_BioCarta[enrichr_BioCarta$Adjusted.P.value < 0.05,]
# Ion Channels and Their Functional Role in Vascular Endothelium

# motifEnrich_CBX3 & BCL6 - biocarta-top::: Ion Channels and Their Functional Role in Vascular Endothelium




message("+   BARPLOT ENRICH R CBX3 BCL6 :::                          ")




enrichr_Wiki$db     <- "WikiPathways2019"
enrichr_BioCarta$db <- "BioCarta_2016"
enrichr_Kegg$db     <- "KEGG_2019_Human"

enrichR_CBX3 <- rbind(enrichr_Wiki, enrichr_)
enrichR_CBX3 <- rbind(enrichR_CBX3, enrichr_Kegg)

selected_terms <- c("Wnt Signaling WP428", "Vascular smooth muscle contraction", "Wnt signaling pathway","Ion Channels and Their Functional Role in Vascular Endothelium Homo sapiens h raccPathway", "cGMP-PKG signaling pathway", "Differentiation Pathway WP2848", "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling WP3850", "Eicosanoid Metabolism Homo sapiens h eicosanoidPathway	", "cGMP-PKG signaling pathway")


enrichR_CBX3_gg <- enrichR_CBX3[enrichR_CBX3$Term %in% selected_terms, ]
enrichR_CBX3_gg$gene_no <- gsub("/.*", "", enrichR_CBX3_gg$Overlap)
enrichR_CBX3_gg$Term <- gsub(" Homo sapiens h raccPathway", "", enrichR_CBX3_gg$Term)
enrichR_CBX3_gg$Term <- gsub("and Their Functional Role", "and Funct. Role", enrichR_CBX3_gg$Term)
enrichR_CBX3_gg <- enrichR_CBX3_gg[-3,]
enrichR_CBX3_gg$Term <- gsub("Wnt Signaling WP428", "Wnt Signaling", enrichR_CBX3_gg$Term)

#Ion Channels and Their Functional Role in Vascular Endothelium Homo sapiens h raccPathway	

p_CBX3_gg <- ggplot(enrichR_CBX3_gg, aes(x=reorder(Term, -Adjusted.P.value), y=as.numeric(gene_no))) + geom_bar(stat="identity", aes(alpha = Combined.Score) , fill = col_2nd) +
  coord_flip() +  xlab("enrichR terms for CBX3 enriched DEGs") +
  #scale_fill_manual( values = c(col_2nd)) +  
  ylab("Gene count") +
  #ylim(-max(abs(enrichKegg_molten$genes_UP)), max(abs(enrichKegg_molten$genes_UP)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1)) +
  scale_y_continuous(breaks=c(0,5,10))


p_CBX3_gg



enrichr_Wiki$db     <- "WikiPathways2019"
enrichr_BioCarta$db <- "BioCarta_2016"
enrichr_Kegg$db     <- "KEGG_2019_Human"

enrichR_BCL6 <- rbind(enrichr_Wiki, enrichr_BioCarta)
enrichR_BCL6 <- rbind(enrichR_BCL6, enrichr_Kegg)

selected_terms <- c("Wnt Signaling WP428", "Vascular smooth muscle contraction", "Wnt signaling pathway","Ion Channels and Their Functional Role in Vascular Endothelium Homo sapiens h raccPathway", "cGMP-PKG signaling pathway", "Differentiation Pathway WP2848", "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling WP3850", "Eicosanoid Metabolism Homo sapiens h eicosanoidPathway	", "Differentiation Pathway WP2848")

enrichR_BCL6_gg <- enrichR_BCL6[enrichR_BCL6$Term %in% selected_terms, ]
enrichR_BCL6_gg$gene_no <- gsub("/.*", "", enrichR_BCL6_gg$Overlap)
enrichR_BCL6_gg$Term <- gsub(" Homo sapiens h raccPathway", "", enrichR_BCL6_gg$Term)
enrichR_BCL6_gg$Term <- gsub("and Their Functional Role", "and Funct. Role", enrichR_BCL6_gg$Term)
enrichR_BCL6_gg$Term[enrichR_BCL6_gg$Term == "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling WP3850"] <- "Factors & pathways affecting IGF1-Akt sig. WP3850"
enrichR_BCL6_gg$Term <- gsub(" WP.*", "", enrichR_BCL6_gg$Term)


#Ion Channels and Their Functional Role in Vascular Endothelium Homo sapiens h raccPathway	

p_BCL6_gg <- ggplot(enrichR_BCL6_gg, aes(x=reorder(Term, -Adjusted.P.value), y=as.numeric(gene_no))) + geom_bar(stat="identity", aes(alpha = Combined.Score) , fill = col_2nd) +
  coord_flip() +  xlab("enrichR terms for BCL6 enriched DEGs") +
  #scale_fill_manual( values = c(col_2nd)) +  
  ylab("Gene count") +
  #ylim(-max(abs(enrichKegg_molten$genes_UP)), max(abs(enrichKegg_molten$genes_UP)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))

p_CBX3_gg
p_BCL6_gg




enrichr_Wiki$db     <- "WikiPathways"
enrichr_BioCarta$db <- "BioCarta"
enrichr_Kegg$db     <- "KEGG"
enrichr_GOBP$db     <- "GO_BP"
enrichr_React$db     <- "Reactome"

enrichR_CBX3_BCL6 <- rbind(enrichr_Wiki, enrichr_GOBP)
enrichR_CBX3_BCL6 <- rbind(enrichR_CBX3_BCL6, enrichr_Kegg)
enrichR_CBX3_BCL6 <- rbind(enrichR_CBX3_BCL6, enrichr_GOBP)
enrichR_CBX3_BCL6 <- rbind(enrichR_CBX3_BCL6, enrichr_React)

selected_terms <- c("Wnt Signaling WP428", "Vascular smooth muscle contraction", "Wnt signaling pathway","Ion Channels and Their Functional Role in Vascular Endothelium Homo sapiens h raccPathway", "cGMP-PKG signaling pathway", "Differentiation Pathway WP2848", "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling WP3850", "Eicosanoid Metabolism Homo sapiens h eicosanoidPathway", "cGMP-PKG signaling pathway", "Vitamin D Receptor Pathway WP2877","Cell adhesion molecules (CAMs)", "Parathyroid hormone synthesis, secretion and action", "Muscle contraction Homo sapiens R-HSA-397014", "regulation of vascular smooth muscle cell proliferation (GO:1904705)")


enrichR_CBX3_BCL6 <- enrichR_CBX3_BCL6[enrichR_CBX3_BCL6$Term %in% selected_terms, ]
enrichR_CBX3_BCL6$gene_no <- gsub("/.*", "", enrichR_CBX3_BCL6$Overlap)
enrichR_CBX3_BCL6$Term <- gsub(" Homo sapiens h.*", "", enrichR_CBX3_BCL6$Term)
enrichR_CBX3_BCL6$Term <- gsub("and Their Functional Role", "and Funct. Role", enrichR_CBX3_BCL6$Term)
enrichR_CBX3_BCL6$Term <- gsub(" WP.*", "", enrichR_CBX3_BCL6$Term)
enrichR_CBX3_BCL6$Term <- gsub("regulation", "Reg.", enrichR_CBX3_BCL6$Term)
enrichR_CBX3_BCL6$Term <- gsub(" \\(GO.*", "", enrichR_CBX3_BCL6$Term)
enrichR_CBX3_BCL6$Term <- gsub("R-HSA.*", "", enrichR_CBX3_BCL6$Term)
enrichR_CBX3_BCL6$Term[enrichR_CBX3_BCL6$Term ==  "Factors and pathways affecting insulin-like growth factor (IGF1)-Akt signaling"] <- "Factors & pathways affecting IGF1-Akt sig."


p_CBX3_BCL6 <- ggplot(enrichR_CBX3_BCL6, aes(x=reorder(Term, -Adjusted.P.value), y=as.numeric(gene_no))) + geom_bar(stat="identity", aes(alpha = Combined.Score) , fill = col_2nd) +
  coord_flip() +  xlab("enrichR terms for CBX3 & BCL6 enriched DEGs") +
  #scale_fill_manual( values = c(col_2nd)) +  
  ylab("Gene count") +
  #ylim(-max(abs(enrichKegg_molten$genes_UP)), max(abs(enrichKegg_molten$genes_UP)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1)) +
  scale_y_continuous(breaks=c(0,5,10))


p_CBX3_BCL6





plot_grid(p_CBX3_gg, p_BCL6_gg, nrow = 2, hjust = T, vjust = T)

pdf(paste(Project, "barplot_CBX3_BCL6",  ".pdf", sep="_"), onefile=FALSE, width=7, height=10) 
par(bg=NA)
plot_grid(p_CBX3_gg, p_BCL6_gg, nrow = 2, hjust = T, vjust = T)
dev.off()

pdf(paste(Project, "barplot_CBX3",  ".pdf", sep="_"), onefile=FALSE, width=7, height=2.5) 
par(bg=NA)
p_CBX3_gg
dev.off()

pdf(paste(Project, "barplot_BCL6",  ".pdf", sep="_"), onefile=FALSE, width=7, height=3) 
par(bg=NA)
p_BCL6_gg
dev.off()

# https://academic.oup.com/nar/article/44/9/4080/2462463
# Systematic identification of gene family regulators in mouse and human embryonic stem cells 
# CNOT, REST, SETDB1, GATA4, TIP60?, REX1?, DMAP1, JARID1a?, MACAF?, ZNF143, ZNF348?, ZFP322a, SMAD2 and CBX3
#      GATA4 - up,   ZNF143 SMAD2 CNOT7P1 - down like CBX3



motifEnrichmentTable_wGenes_wLogo <- addLogo(motifEnrichmentTable_wGenes)
motifEnrichmentTable_wGenes_wLogoTSS <- addLogo(motifEnrichmentTable_wGenesTSS)
resultsSubset <- motifEnrichmentTable_wGenes_wLogo[1:10,]
resultsSubsetTSS <- motifEnrichmentTable_wGenes_wLogoTSS[1:10,]

library(DT)
datatable(resultsSubsetTSS[,-c("enrichedGenes", "TF_lowConf"), with=TRUE], 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=5))

datatable(resultsSubsetTSS, 
          escape = FALSE, # To show the logo
          filter="top", options=list(pageLength=10))






anotatedTfs <- lapply(split(motifEnrichmentTable_wGenes$enrichedGenes, motifEnrichmentTable_wGenes$motif),
                      function(x) {
                        #genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                        genesSplit <- unique(unlist(strsplit(x, ";")))
                        return(genesSplit)
                      })

anotatedTfs$dbcorrdb__CBX3__ENCSR000BRT_1__m9


#Building a network
signifMotifNames <- motifEnrichmentTable_wGenesTSS$motif[1:5]
signifMotifNames <- c("dbcorrdb__CBX3__ENCSR000BRT_1__m9", "transfac_pro__M09104")

getSignificantGenes(geneLists, 
                    motifRankings,
                    signifRankingNames=signifMotifNames,
                    plotCurve=TRUE, maxRank=5000, genesFormat="none",
                    method="iCisTarget")




incidenceMatrix <- getSignificantGenes(geneLists, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="iCisTarget")$incidMatrix # iCisTarget or aprox

rownames(incidenceMatrix) <- c("CBX3", "BCL6")

library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
#Output not shown:
  
library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visConfigure(enabled = TRUE)  %>% visOptions(highlightNearest = TRUE,  nodesIdSelection = TRUE)
head(edges)
head(nodes)

server <- function(input, output) {
  output$network <- renderVisNetwork({
    # minimal example
    visNetwork(nodes, edges)  %>% visOptions(highlightNearest = TRUE,  nodesIdSelection = TRUE)
  })
}

ui <- fluidPage( visNetworkOutput("network"))

shinyApp(ui = ui, server = server)


#network <- visNetwork(nodes, edges, width = "100%")
#network %>% visSave(file = "network.html")
# same as
#visSave(network, file = "network.html")
# or
#htmlwidgets::saveWidget(network, "network.html")

pdf(paste(Project, "visNetwork",  ".pdf", sep="_"), onefile=FALSE, width=10, height=10) 
par(bg=NA)
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE,  nodesIdSelection = TRUE)

dev.off()







library(httr)
library(jsonlite)

#res = GET("http://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/CBX3/ENCODE+Transcription+Factor+Targets")

CBX3_EncodeTFtargets <- read.table("amp.pharm.mssm.edu_CBX3_EncodeTFtargets.txt", sep = ",")
#scan("file.txt", what = "", quiet=TRUE) 
CBX3_EncodeTFtargets <- strsplit(as.character(CBX3_EncodeTFtargets), split = "symbol")
CBX3_EncodeTFtargets[[6]]
CBX3_EncodeTFtargets <- CBX3_EncodeTFtargets[grep("href:/api/1.0/gene/",CBX3_EncodeTFtargets)]
CBX3_EncodeTFtargets <- gsub( "href:/api/1.0/gene/", "", CBX3_EncodeTFtargets)
CBX3_EncodeTFtargets <- gsub( "}", "", CBX3_EncodeTFtargets)
length(CBX3_EncodeTFtargets) # 11894

length(DEG_CBX3$external_gene_name[DEG_CBX3$external_gene_name %in% CBX3_EncodeTFtargets]) # 68 /237
length(DEG_CBX3tss$external_gene_name[DEG_CBX3tss$external_gene_name %in% CBX3_EncodeTFtargets]) # 208 /697

DEG_CBX3_BCL6 <- DEG_CBX3tss$external_gene_name[DEG_CBX3tss$external_gene_name %in% DEG_BCL6_tss$external_gene_name]
length(DEG_CBX3_BCL6[DEG_CBX3_BCL6 %in% CBX3_EncodeTFtargets]) # 65 / 346








message("+              Epi  db                            ")

Epi_db <- read.csv("EpiGenes_main_cur_2018_11_21.csv", header = TRUE)
Epi_db$external_gene_name <-  Epi_db$HGNC_symbol

Epi_genes <- resSig.ann[resSig.ann$external_gene_name %in% Epi_db$external_gene_name  ,]
Epi_DEGs <- Epi_db[Epi_db$external_gene_name %in% resSig.ann$external_gene_name  ,]
Epi_DEGs <- resSig.ann[resSig.ann$external_gene_name %in% Epi_db$external_gene_name  ,]
Epi_DEGs <- Epi_DEGs[abs(Epi_genes$log2FoldChange)>1,]

length(DEG_CBX3[DEG_CBX3$external_gene_name %in% Epi_DEGs$external_gene_name,]$external_gene_name) # 6/237 PRDM6   GADD45B ZBTB16  SCML4   PRDM8   PRKCB  
length(DEG_CBX3tss[DEG_CBX3tss$external_gene_name %in% Epi_DEGs$external_gene_name,]$external_gene_name) # 9 / 697  PRDM6  ZBTB16 SCML4  PRDM8  GFI1B  PRKCB  HR     SATB1  IKZF1 


motifEnrich_genes
motifEnrich_genes[motifEnrich_genes %in% Epi_genes$external_gene_name] # 0





message("+            CEMiTool                   ")


# module load R/3.6.2

#‘mnormt’, ‘conquer’

library("CEMiTool")
Project       <- "CTR_gjb2_0001_STAR_DESeq2_shrinkage_BlockSex"
significance  <- 0.05
l2fc          <- 1
col_1st       <- "firebrick2"
col_2nd       <- "steelblue3"

Base.dir      <- "/storage/CTR-Projects/CTR_gjb2/CTR_gjb2_0001/BACK_UP_2020/CTR_gjb2_0001"
setwd(Base.dir)

resSig.ann <- read.csv("CTR_gjb2_0001_STAR_DESeq2_shr_BlockSex_deseq2_DEGs_padj0.05.csv")
#expr0 <- readRDS("CTR_gjb2_0001_exprMatr_DEG.Rds")
#expr0 <- readRDS("CTR_gjb2_0001_exprMat.Rds")

expr0 <- readRDS("rld_df_ann.Rds")
expr0 <- expr0[,-1]
expr0$rowsums <- rowSums(expr0[,c(1:14)])
expr0 <- expr0[order(expr0$rowsums, decreasing = TRUE),]
expr0 <- expr0[!duplicated(expr0$external_gene_name),]
rownames(expr0) <- expr0$external_gene_name
expr0 <- expr0[,c(1:14)]





head(expr0)
sample_annot <- as.data.frame(colnames(expr0))
colnames(sample_annot)[1] <- "SampleName"
sample_annot$Class <- gsub( "First_trimester.*", "FT",sample_annot$SampleName)
sample_annot$Class <- gsub( "Second_trimester.*", "ST",sample_annot$Class)


# run cemitool with sample annotation

#expr0 <- expr0[rownames(expr0) %in% resSig.ann$external_gene_name,]
#cem <- readRDS("cem.Rds")

head(sample_annot)
cem <- cemitool(expr0, sample_annot, force_beta = TRUE)
cem

# Module inspection
nmodules(cem)

head(module_genes(cem))

options(bitmapType='cairo')
generate_report(cem, force=TRUE)
write_files(cem, force=TRUE)
save_plots(cem, "all", force=TRUE)


# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")

#Expression patterns in modules
#You can generate a plot that displays the expression of each gene within a module using the plot_profile function:
# plot gene expression within each module
cem <- plot_profile(cem)
x11(type ="cairo")

plots <- show_plot(cem, "profile")
plots[1]


show_plot(cem, "gsea")
plots_interactions <- show_plot(cem, "interaction")
#"profile", "gsea", "ora", "interaction", "beta_r2", "mean_k", "sample_tree", "mean_var", "hist", "qq".
plots_ora <- show_plot(cem, "ora")
plots_ora[1]


#Adding ORA analysis
x11(type ="cairo")

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

# plot ora results
cem <- plot_ora(cem)
plots_ora <- show_plot(cem, "ora")
plots_ora[1]
#ggsave("M1_ora.png", device = CairoPNG)

#plots <- readRDS("plots_ora.Rds")
#plots[1]

# Adding interactions

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)


# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots_inter[1]





# run cemitool ALLL ::::
# 
library(ggplot2)
cem <- cemitool(expr0, sample_annot, gmt_in, interactions=int_df, filter=TRUE, plot=TRUE, verbose=TRUE)
# create report as html document
generate_report(cem, directory="./Report",force=TRUE)

# write analysis results into files
write_files(cem, directory="./Tables")

# save all plots
save_plots(cem, "all", directory="./Plots")



modules <- read.table("CemiTools/Tables/module.tsv", header = T)
modules_list <- split(modules$genes, modules$modules)
# interesting modules: 
# 1 -  > 2000genes
# 8 - hemoglobin/immune sstem
# 7- ER, unfolded protein resposnse, preeclampsia genes
# 9/6/5- not interesting
# 4 - signaling activity, receptors, GPCR activity, TAP complex
# 3 - anatomical structure development/ wnt pahtway, sequence -speific DNA binding activity (lots of TFs)
# 2 - ECM, lipid syntheiss/transport, MAPK/ , cell surface interactions at vascular wall
# 10 - innate immune response, chemotaxis, IL17 signaling pathway, lots of secreted proteins, chemokine sig path, cytokine-cytokne receptor interactions
# 







sessionInfo()
