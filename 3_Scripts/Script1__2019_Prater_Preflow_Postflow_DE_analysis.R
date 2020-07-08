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

Base.dir      <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/"
setwd(Base.dir)
Res.dir      <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR/BlockSex"

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
resSig.annl2fc1 <- subset(resSig.ann, abs(resSig.ann$log2FoldChange) > 1)


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
#  2nd: "steelblue3"
#  1st: "firebrick2"

#resdata_ann

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

Fig_1_volcano

pdf(paste("Fig_1b_volcano", Project,  "corrLabs.pdf", sep="_"), onefile=FALSE, width=7, height=7) 
par(bg=NA)
Fig_1_volcano
dev.off()


message("+-------------------------------------------------------------------------------")
message("+                             Starting GO/KEGG                                "  )
message("+-------------------------------------------------------------------------------")

RESULTS_1 <- subset(results.df, abs(results.df$log2FoldChange) > l2fc & results.df$padj < significance)
RESULTS_0.6 <- subset(results.df, abs(results.df$log2FoldChange) > 0.6 & results.df$padj < significance)
resdata_simplified <- resdata_ann
resdata_simplified$mean_1st <- rowMeans(resdata_simplified[,c(8:15)])
resdata_simplified$mean_2nd <- rowMeans(resdata_simplified[,c(16:21)])
resdata_simplified <- resdata_simplified[,c(1,2,3,7,22:25)]


message("+-------------------------------------------------------------------------------")
message("+                                  enrichKegg                                   ")
message("+-------------------------------------------------------------------------------")

Kegg_genes <- na.omit(RESULTS_0.6)
head(Kegg_genes)
colnames(Kegg_genes)[colnames(Kegg_genes)=="entrezgene_id"] <- "Gene_ID" 

foldchanges = Kegg_genes$log2FoldChange
names(foldchanges) = Kegg_genes$Gene_ID
foldchanges <- sort(foldchanges, decreasing = T)
head(foldchanges)

kk <- enrichKEGG(Kegg_genes$Gene_ID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05, keyType = "kegg")
head(summary(kk))
kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kk_results <- as.data.frame(kk2)
keggresids = kk_results$ID

selected_kegg <- c("hsa04514" ,"hsa04060", "hsa04022",  "hsa04015", "hsa04141", "hsa04913", "hsa04921", "hsa04062", "hsa04350", "hsa04270", "hsa04010","hsa04911" ,"hsa03320" ,"hsa04918" ,"hsa00510",  "hsa04923", "hsa02010", "hsa04310", "hsa04141", "hsa04151", "hsa04072","hsa04020", "hsa04929", "hsa04512")   
  
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
enrichKegg_molten <- melt(enrichKegg_selected[,c(1,2,6,7,10:11)], id.vars=c("ID","Description","p.adjust", "qvalue") )

p_kegg_mlt <- ggplot(enrichKegg_molten, aes(x=reorder(Description, -qvalue), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue))) +
  coord_flip() +  xlab("KEGG Pathways") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +  ylab("Gene count") +
  ylim(-max(abs(enrichKegg_molten$value)), max(abs(enrichKegg_molten$value)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))

p_kegg_mlt



message("+-------------------------------------------------------------------------------")
message("+                              Cluster profiler                                 ")
message("+-------------------------------------------------------------------------------")

bkcg_genes <- as.character(unique(results.df$entrezgene_id))
geneListx <- unique(resSig.ann[,c(9,3)])
geneListx <- geneListx[!is.na(geneListx$entrezgene),]
geneList <- geneListx$log2FoldChange
names(geneList) <- geneListx$entrezgene
geneList <- sort(geneList, decreasing = T )
geneList <- geneList[unique(names(geneList))]
gene <- names(geneList)[abs(geneList) > 1]

# GSEA = gene set enrichment analysis
gsego_cc <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "CC", nPerm= 1000, minGSSize  = 20, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsego_bp <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "BP", nPerm= 1000, minGSSize  = 20, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsego_mf <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "MF", nPerm= 1000, minGSSize  = 20, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsecc_df <- as.data.frame(gsego_cc)
gsebp_df <- as.data.frame(gsego_bp)
gsemf_df <- as.data.frame(gsego_mf)



message("+-------------------------------------------------------------------------------")
message("+                              Plot gse   results                               ")
message("+-------------------------------------------------------------------------------")

gsebp_selected_terms <-  c("ER to Golgi vesicle-mediated transport", "response to endoplasmic reticulum stress", "cytokine-mediated signaling pathway",
                           "adaptive immune response", "regulation of vasculature development", "cellular response to oxygen-containing compound", 
                           "G protein-coupled receptor signaling pathway", "positive regulation of ERK1 and ERK2 cascade",  "carboxylic acid transport",
                           "positive regulation of MAPK cascade", "cellular response to lipid", "calcium ion transport", 
                           "positive regulation of peptide secretion", "regulation of transmembrane transport", "Golgi vesicle transport",  "regulation of GTPase activity", "positive regulation of cell adhesion", "regulation of cytosolic calcium ion concentration", "response to oxidative stress", "immune response-activating signal transduction",  "extracellular matrix organization", "IRE1-mediated unfolded protein response", "endoplasmic reticulum unfolded protein response", "regulation of complement activation", "negative regulation of response to endoplasmic reticulum stress", "cell redox homeostasis", 'ribosome biogenesis') 

gsemf_selected_terms <-  c("molecular transducer activity", "signaling receptor activity", "transmembrane signaling receptor activity" , "calcium ion binding" , "extracellular matrix structural constituent", "receptor ligand activity","substrate-specific channel activity" , "receptor regulator activity" , "passive transmembrane transporter activity", "cation channel activity","G protein-coupled receptor activity", "cytokine activity")

gsecc_selected_terms <- c("secretory vesicle", "cell surface" ,"vesicle membrane" , "transmembrane transporter complex", "external side of plasma membrane", "ion channel complex" ,
                          "secretory granule membrane", "transporter complex", "collagen-containing extracellular matrix" ,"plasma membrane protein complex", "ribosomal subunit", "catalytic step 2 spliceosome"  )     


#gse_BP_selected <- gsebp_df
gse_BP_selected <- gsebp_df[gsebp_df$Description %in% gsebp_selected_terms,]
gse_BP_selected$Description <- gsub("G protein-coupled receptor", "GPCR" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("regulation", "reg." , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("oxygen", "O2" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("nucleotide", "nt" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("signaling pathway", "signaling" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("growth factor stimulus", "growth factor" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("endoplasmic reticulum", "ER" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("immune response-activating signal transduction", "immune response-activ. signal transduction" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains", "adaptive immune response" , gse_BP_selected$Description)

gene_counts <- strsplit(gse_BP_selected$core_enrichment, split = "/")
gene_counts_list <- list()
for (i in seq_along(gene_counts)){ gene_counts_list[[i]] <- length(gene_counts[[i]])}
gse_BP_selected$gene_count <- gene_counts_list
gse_BP_selected$group <- ifelse(gse_BP_selected$enrichmentScore < 0 , "Preflow", "Postflow")

p_bp <- ggplot(gse_BP_selected, aes(x=reorder(Description, -qvalues), y=gene_count, fill=group)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues))) +
  coord_flip() + xlab("GSE: Biological Processes") +  scale_fill_manual( values = c( "Postflow"=col_2nd, "Preflow" =col_1st)) +
  ylab("Gene count") +  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1)) +  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))


list_up <- list()
list_down <- list()
for (i in 1:nrow(gse_BP_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$entrezgene_id %in% unlist(strsplit(gse_BP_selected[i, "core_enrichment"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <-tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <-tmp_down
  
}
gse_BP_selected$genes_UP   <- as.numeric(list_up)
gse_BP_selected$genes_DOWN <- -as.numeric(list_down)
BP_molten <- melt(gse_BP_selected[,c(1,2,7,8,14:15)], id.vars=c("ID","Description","p.adjust", "qvalues") )

p_BP_mlt <- ggplot(BP_molten, aes(x=reorder(Description, -qvalues), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues))) +
  coord_flip() +  xlab("GSE: Biological Processes") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +  ylab("Gene count") +
  ylim(-max(abs(BP_molten$value)), max(abs(BP_molten$value)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))




gse_MF_selected <- gsemf_df[gsemf_df$Description %in% gsemf_selected_terms,]
gse_MF_selected$Description <- gsub("transmembrane", "transmem." , gse_MF_selected$Description)
gse_MF_selected$Description <- gsub("G protein-coupled receptor", "GPCR" , gse_MF_selected$Description)

gene_counts <- strsplit(gse_MF_selected$core_enrichment, split = "/")
gene_counts_list <- list()
for (i in seq_along(gene_counts)){ gene_counts_list[[i]] <- length(gene_counts[[i]])}
gse_MF_selected$gene_count <- gene_counts_list
gse_MF_selected$group <- ifelse(gse_MF_selected$enrichmentScore < 0 , "Preflow", "Postflow")

p_mf <- ggplot(gse_MF_selected, aes(x=reorder(Description, -qvalues), y=gene_count, fill=group)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues))) + coord_flip()+ xlab("GSE: Molecular Function") + ylab("Gene count") +
  scale_fill_manual( values = c( "Postflow"=col_2nd, "Preflow" =col_1st)) + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_alpha_continuous( range = c(0.5, 1)) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))


list_up <- list()
list_down <- list()
for (i in 1:nrow(gse_MF_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$entrezgene_id %in% unlist(strsplit(gse_MF_selected[i, "core_enrichment"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <-tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <-tmp_down
  
}

gse_MF_selected$genes_UP <- as.numeric(list_up)
gse_MF_selected$genes_DOWN <- -as.numeric(list_down)
MF_molten <- melt(gse_MF_selected[,c(1,2,7,8,14:15)], id.vars=c("ID","Description","p.adjust", "qvalues") )

p_MF_mlt <- ggplot(MF_molten, aes(x=reorder(Description, -qvalues), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues))) +
  coord_flip() +  xlab("GSE: Molecular Function") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +  ylab("Gene count") +
  ylim(-max(abs(MF_molten$value)), max(abs(MF_molten$value)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))




gse_CC_selected <- gsecc_df[gsecc_df$Description %in% gsecc_selected_terms,]
gse_CC_selected$Description <- gsub("extracellular matrix", "ECM" , gse_CC_selected$Description)

gene_counts <- strsplit(gse_CC_selected$core_enrichment, split = "/")
gene_counts_list <- list()
for (i in seq_along(gene_counts)){ gene_counts_list[[i]] <- length(gene_counts[[i]])}
gse_CC_selected$gene_count <- gene_counts_list
gse_CC_selected$group <- ifelse(gse_CC_selected$enrichmentScore < 0 , "Preflow", "Postflow")

p_cc <- ggplot(gse_CC_selected, aes(x=reorder(Description, -qvalues), y=gene_count, fill=group)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues)))+
  coord_flip() + xlab("GSE: Cellular Component") +
  scale_fill_manual( values = c( "Postflow"=col_2nd, "Preflow" =col_1st)) +  ylab("Gene count") +
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1)) +   theme(axis.text=element_text(size=12), axis.title=element_text(size=14))

list_down <- list()
list_up <- list()
for (i in 1:nrow(gse_CC_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$entrezgene_id %in% unlist(strsplit(gse_CC_selected[i, "core_enrichment"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <-tmp_up
  tmp_down <- length(subset(df_tmp[,3], df_tmp[,3] < 0))
  list_down[[i]] <-tmp_down
  
}
gse_CC_selected$genes_UP <- as.numeric(list_up)
gse_CC_selected$genes_DOWN <- -as.numeric(list_down)
CC_molten <- melt(gse_CC_selected[,c(1,2,7,8,14:15)], id.vars=c("ID","Description","p.adjust", "qvalues") )

p_CC_mlt <- ggplot(CC_molten, aes(x=reorder(Description, -qvalues), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues))) +
  coord_flip() +  xlab("GSE: Cellular Component") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +  ylab("Gene count") +
  ylim(-max(abs(CC_molten$value)), max(abs(CC_molten$value)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))


#plot_GO_Kegg_2sided <- plot_grid(p_kegg_mlt, p_BP_mlt, p_MF_mlt, p_CC_mlt,  labels=c("", "", "", ""), align="hv", ncol = 2, nrow = 2, rel_heights = c(1, 0.6))



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


sessionInfo()
