#!/usr/local/bin/Rscript
rm(list=ls())

##############################################################################################################
#   Project: CTR_gjb2_0001_placenta_1st_vs_2nd_trimester   
#                 
#   Malwina Prater (mn367@cam.ac.uk), 2019                     
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

Project       <- "CTR_gjb2_0001_STAR_DESeq2_shrinkage"
significance  <- 0.05
l2fc          <- 1
col_1st <- "firebrick2"
col_2nd <- "steelblue3"

Base.dir      <- "/Users/malwina/Documents/CTR-Groups/Graham_Burton/CTR_gjb2_0001/STAR"
setwd(Base.dir)

message("+-------------------------------------------------------------------------------")
message("+                       Prepare sample table                                    ")
message("+-------------------------------------------------------------------------------")

sampleTable  <- read.csv("sampleTable.csv")
sample_id    <- sampleTable$sample
HTSeq.dir    <- paste(Base.dir,"/HTSEQ-COUNTS_CTR_gjb2_0001", sep="")
list.files(HTSeq.dir)

message("+-------------------------------------------------------------------------------")
message("+                       Retrieve ensEMBL annotations                            ")
message("+-------------------------------------------------------------------------------")

ensembl    =  useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description'), mart = ensembl)          
ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'description', 'go_id', 'name_1006' ), mart = ensembl)    
message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ condition)
dds <- DESeq(ddsHTSeq)

message("+-------------------------------------------------------------------------------")
message("+  Get the RESULTS    ")
message("+-------------------------------------------------------------------------------")

res <- lfcShrink(dds, coef="condition_Second_trimester_vs_First_trimester", type="normal")
res <-res[order(res$padj),]

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
write.csv(resSig.ann[order((resSig.ann$log2FoldChange), decreasing = TRUE),],   file=paste(Project, '_deseq2_DEGs_padj0.05.csv', sep=""))

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized=T)), by='row.names',sort=F)
names(resdata)[1] <- 'ensembl_gene_id'
resdata_ann <- merge(resdata, ensEMBL2id[,c(1:3)], by.x= "ensembl_gene_id")


message("+-------------------------------------------------------------------------------")
message("+ Run transformations")
message("+-------------------------------------------------------------------------------")

rld <- DESeq2::rlogTransformation(dds, blind=T)

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

Supp_Fig_1_A <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sampleTable$condition)) )) + 
  geom_point(size = 5 ) + 
  xlab(pc1lab) + ylab(pc2lab) + 
  geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) + 
  scale_fill_manual(name="Stage", values = c(col_1st, col_2nd)) +
  scale_colour_manual(name="Stage", values = c(col_1st, col_2nd)) +
  theme(text = element_text(size=elementTextSize)) + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=14))

Supp_Fig_1_A_labs <- ggplot(scores, aes(x = PC1, y = PC2, col = (factor(sampleTable$condition)) )) + 
  geom_point(size = 5 ) + geom_text_repel(aes(label=sampleTable$Sample_short), col = "black") +
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

pdf(paste(Project, "_SuppFig1_PCA_Loadings_clustering_all_genes.pdf", sep=""),width=10,height=14)
par(bg=NA)
plot_grid(Supp_Fig_1_A_labs, Supp_Fig_1_B_C, labels=c("A", ""), ncol = 1, nrow = 2, rel_widths = c(1.2, 1))
dev.off()


message("+-------------------------------------------------------------------------------")
message("+                                volcano plot                                   ")
message("+-------------------------------------------------------------------------------")
#  2nd: "steelblue3"
#  1st: "firebrick2"

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
  xlab(expression(log[2]("First_trimester" / "Second_trimester"))) + 
  ylab(expression(-log[10]("adjusted p-value"))) +   ggtitle(label = "" , subtitle = "") +  
  geom_vline(xintercept = 0, colour = "black") +  geom_hline(yintercept = 1.3, colour = "black") + 
  annotate(geom = "text", label = "First trimester", x = -4, y = 55, size = 7, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 4, y = 55, size = 7, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) +  
  theme(text = element_text(size=elementTextSize)) 

Fig_1_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 3 & padj < 0.00000001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"))

pdf(paste("Fig_1b", "_preflow_postflow_volcano_plot_padj0.05", ".pdf", sep=""), width=7, height=7, onefile=FALSE)
par(bg=NA)
volcano_plot2
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

selected_kegg <- c("hsa04514" ,"hsa04060", "hsa04022",  "hsa04015", "hsa04141", "hsa04913", "hsa04921", "hsa04062", "hsa04350", "hsa04270", "hsa04010","hsa04911" ,"hsa03320" ,"hsa04918" ,"hsa00510",  "hsa04923", "hsa02010", "hsa04310", "hsa04141", "hsa04151")

enrichKegg_selected <- kk_results[kk_results$ID %in% selected_kegg,]
enrichKegg_selected$Description <- gsub("endoplasmic reticulum", "ER"  , enrichKegg_selected$Description)

list_up <- list()
for (i in 1:nrow(enrichKegg_selected)){
  df_tmp <- resdata_simplified[resdata_simplified$external_gene_name %in% unlist(strsplit(enrichKegg_selected[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,3], df_tmp[,3] > 0))
  list_up[[i]] <-tmp_up
}
enrichKegg_selected$genes_UP <- as.numeric(list_up)
enrichKegg_selected$genes_DOWN <- enrichKegg_selected$Count - enrichKegg_selected$genes_UP
enrichKegg_selected$genes_DOWN <- -enrichKegg_selected$genes_DOWN

enrichKegg_molten <- melt(enrichKegg_selected[,c(1,2,6,7,10:11)], id.vars=c("ID","Description","p.adjust", "qvalue") )


p_kegg_mlt <- ggplot(enrichKegg_molten, aes(x=reorder(Description, -qvalue), y=value, fill=variable)) + geom_bar(stat="identity", aes(alpha = -log2(qvalue))) +
  coord_flip() +  xlab("KEGG Pathways") +
  scale_fill_manual( values = c(col_2nd, col_1st)) +  ylab("Gene count") +
  ylim(-max(abs(enrichKegg_molten$value)), max(abs(enrichKegg_molten$value)))+
  theme(legend.position = "none") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1))

pdf(paste("Fig_1c", "enrichKEGG_selected", "barplot20_padj0.05_l2fc0.6_MLT", ".pdf", sep="_"), onefile=FALSE, width=10, height=7) 
par(bg=NA)
p_kegg_mlt 
dev.off()


message("+-------------------------------------------------------------------------------")
message("+                              Cluster profiler                                 ")
message("+-------------------------------------------------------------------------------")

bkcg_genes <- as.character(unique(results.df$entrezgene_id))
geneListx <- unique(resSig.ann[,c(9,3)])
geneListx <- geneListx[!is.na(geneListx$entrezgene),]
geneList <- geneListx$log2FoldChange
names(geneList) <- geneListx$entrezgene
geneList <- sort(geneList, decreasing = T )
gene <- names(geneList)[abs(geneList) > 1]

# GSEA = gene set enrichment analysis
gsego_cc <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "CC", nPerm= 1000, minGSSize  = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsego_bp <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "BP", nPerm= 1000, minGSSize  = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsego_mf <- gseGO(geneList = geneList, OrgDb= org.Hs.eg.db, ont = "MF", nPerm= 1000, minGSSize  = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
gsecc_df <- as.data.frame(gsego_cc)
gsebp_df <- as.data.frame(gsego_bp)
gsemf_df <- as.data.frame(gsego_mf)

# GO over-representation test
ego_bp <- enrichGO(gene = gene, universe = bkcg_genes, OrgDb = org.Hs.eg.db, ont = "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
ego_mf <- enrichGO(gene = gene, universe = bkcg_genes, OrgDb = org.Hs.eg.db, ont = "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
ego_cc <- enrichGO(gene = gene, universe = bkcg_genes, OrgDb = org.Hs.eg.db, ont = "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable      = TRUE)
ego <- rbind(as.data.frame(ego_bp), as.data.frame(ego_cc), as.data.frame(ego_mf) )

message("+-------------------------------------------------------------------------------")
message("+                              Plot gse   results                               ")
message("+-------------------------------------------------------------------------------")

gsebp_selected_terms <-  c("ER to Golgi vesicle-mediated transport", "response to endoplasmic reticulum stress", "ncRNA processing",
                           "adaptive immune response", "regulation of vasculature development", "cellular response to oxygen-containing compound", 
                           "leukocyte chemotaxis", "G protein-coupled receptor signaling pathway", "positive regulation of ERK1 and ERK2 cascade", 
                           "positive regulation of MAPK cascade", "cellular response to lipid", "cytokine secretion", 
                           "positive regulation of protein secretion", "regulation of transmembrane transport", "Golgi vesicle transport",
                           "carboxylic acid transport", "regulation of Wnt signaling pathway", "cytokine-mediated signaling pathway") 
gsemf_selected_terms <-  c("molecular transducer activity","signaling receptor activity" , "transmembrane signaling receptor activity" , 
                           "calcium ion binding" , "receptor regulator activity" ,"receptor ligand activity", "channel activity" , 
                           "passive transmembrane transporter activity" , "substrate-specific channel activity" , "gated channel activity" ,
                           "cation channel activity" , "ion gated channel activity" , "G protein-coupled receptor activity" ,
                           "G protein-coupled receptor binding", "cytokine receptor binding" , "glycosaminoglycan binding" ,"sulfur compound binding",
                           "metal ion transmembrane transporter activity")
gsecc_selected_terms <- c("secretory vesicle","secretory granule" , "vesicle membrane" , "cell surface" ,"cytoplasmic vesicle membrane" ,"plasma membrane protein complex",
                          "extracellular matrix" ,  "side of membrane"  , "collagen-containing extracellular matrix", "secretory granule membrane"  , 
                          "external side of plasma membrane" ,"transporter complex", "transmembrane transporter complex" ,  "ion channel complex" ,
                          "synapse part" , "mitochondrial protein complex" , "receptor complex" , "endocytic vesicle" )     


gse_BP_selected <- gsebp_df[gsebp_df$Description %in% gsebp_selected_terms,]
gse_BP_selected$Description <- gsub("G protein-coupled receptor", "GPCR" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("regulation", "reg." , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("oxygen", "O2" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("nucleotide", "nt" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("cellular", "*" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("signaling pathway", "signaling" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("growth factor stimulus", "growth factor" , gse_BP_selected$Description)
gse_BP_selected$Description <- gsub("endoplasmic reticulum", "ER" , gse_BP_selected$Description)

gene_counts <- strsplit(gse_BP_selected$core_enrichment, split = "/")
gene_counts_list <- list()
for (i in seq_along(gene_counts)){ gene_counts_list[[i]] <- length(gene_counts[[i]])}
gse_BP_selected$gene_count <- gene_counts_list
gse_BP_selected$group <- ifelse(gse_BP_selected$enrichmentScore < 0 , "Preflow", "Postflow")

p_bp <- ggplot(gse_BP_selected, aes(x=reorder(Description, -qvalues), y=gene_count, fill=group)) + geom_bar(stat="identity", aes(alpha = -log2(qvalues))) +
  coord_flip() + xlab("GSE: Biological Processes") +  scale_fill_manual( values = c( "Postflow"=col_2nd, "Preflow" =col_1st)) +
  ylab("Gene count") +  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_alpha_continuous( range = c(0.5, 1)) +  theme(axis.text=element_text(size=12), axis.title=element_text(size=14))



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




Fig_1AB_volcano <- plot_grid(Supp_Fig_1_A, Fig_1_volcano, labels=c("A", "B"),  ncol = 2, nrow = 1, scale = c(0.98, 0.98), rel_widths = c(1.25, 1), align="h")

plot_GO_Kegg_2sided <- plot_grid(p_kegg_mlt, p_bp, p_mf, p_cc,  labels=c("C", "D", "E", "F"), align="hv", ncol = 2, nrow = 2, rel_heights = c(1, 1))


pdf(paste("Fig_1",Project, "2nd_1st_plac", "padj", significance,"l2fc1", "PCA_volcano_GSE_and_Kegg_barplots_2sidedxx", "alpha.pdf", sep="_"), onefile=FALSE, width=15, height=16) 
par(bg=NA)
plot_grid(Fig_1AB_volcano, plot_GO_Kegg_2sided, labels=c(" ", " "), align="hv", ncol = 1, nrow = 2, scale = 0.95, rel_heights = c(1, 1.4))
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
saveRDS(results_transport$external_gene_name, "All_Transport_GO_genes.rds")

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
  xlab(expression(log[2]("First_trimester" / "Second_trimester"))) +  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_vline(xintercept = 0, colour = "black") +  geom_hline(yintercept = 1.3, colour = "black") + 
  annotate(geom = "text", label = "First trimester", x = -4, y = 30, size = 6, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 3, y = 30, size = 6, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) + 
  theme(text = element_text(size=elementTextSize)) 

Fig_4_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 2 & padj < 0.00000001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.iter = 10, force = TRUE)

pdf(paste("Fig_4_", "volcano_plot_TRANSPORT", ".pdf", sep=""), width=7, height=7, onefile=FALSE)
par(bg=NA)
Fig_4_volcano
dev.off()


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
  annotate(geom = "text", label = "First trimester", x = -3, y = 35, size = 6, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 3, y = 35, size = 6, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) +  
  theme(text = element_text(size=elementTextSize)) 

Fig_5_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 2.2 & padj < 0.00000001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.iter = 10, force = T)

pdf(paste("Fig_5_", "volcano_plot_PROLIF_DIFF", ".pdf", sep=""), width=6, height=6, onefile=FALSE)
par(bg=NA)
Fig_5_volcano
dev.off()




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
  annotate(geom = "text", label = "First trimester", x = -3, y = 22, size = 6, colour = "black") + 
  annotate(geom = "text", label = "Second trimester", x = 3, y = 22, size = 6, colour = "black") + 
  scale_color_manual(values = c("Second_trimester" = col_2nd, "First_trimester" = col_1st, "none" = "#636363")) +  
  theme(text = element_text(size=elementTextSize)) 

Fig_3_volcano <- colored + geom_text_repel(data=subset(data, abs(lfc) > 1 & padj < 0.001), mapping = aes(label = symbol), size = 4, color = 'black', box.padding = unit(0.3, "lines"), point.padding = unit(0.5, "lines"), max.iter = 10)

pdf(paste(Project, "_preflow_postflow___HORMONE___volcano_plot_padj0.05", ".pdf", sep=""), width=7, height=7, onefile=FALSE)
par(bg=NA)
Fig_3_volcano
dev.off()



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
#saveRDS(selected_genes, "ER_genes_stress_kegg_go.rds")

message("+-------------------------------------------------------------------------------")
message("+          Fig. 2b     secretory Genes                                          ")
message("+-------------------------------------------------------------------------------")

ego_terms_go_secretion <- ensEMBL2id_go
ego_terms_go_secretion$secretion <- grepl("secret", ego_terms_go_secretion$name_1006)
ego_terms_go_secretion <- subset(ego_terms_go_secretion, ego_terms_go_secretion$secretion == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% ego_terms_go_secretion$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)
#saveRDS(selected_genes2$external_gene_name, "secretory_genes_ego.rds")


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
selected_genes <- subset(selected_genes, abs(selected_genes$log2FoldChange) > 1)
#saveRDS(selected_genes$external_gene_name, "oxy_genes_ego.rds")

message("+-------------------------------------------------------------------------------")
message("+             Fig. TFs :                                            ")
message("+-------------------------------------------------------------------------------")

TF_db <- read.csv("/Users/malwina/Documents/CTR-Groups/Graham_Burton/2019_11_07__tfcheckpoint_database.csv")
TF_genes <- resSig.ann[resSig.ann$external_gene_name %in% TF_db$Gene_Name,]
TF_genes <- subset(TF_genes, abs(TF_genes$log2FoldChange) > 0.6)
as.character(TF_genes$external_gene_name)
#saveRDS(TF_genes$external_gene_name, "TF_genes_tf_checkpoint_db.rds")


message("+-------------------------------------------------------------------------------")
message("+              Fig.4e- ECM related genes                                            ")
message("+-------------------------------------------------------------------------------")

Hsa04974_genes <- strsplit(kk_results[kk_results$ID == "hsa04974",]$geneID, split = "/")
ensEMBL2id_go_ECM <- ensEMBL2id_go
ensEMBL2id_go_ECM <- subset(ensEMBL2id_go_ECM, ensEMBL2id_go_ECM$go_id == "GO:0022617" | ensEMBL2id_go_ECM$go_id == "GO:0030198" |ensEMBL2id_go_ECM$go_id == "GO:0017090"  )
ECM_genes <- unique(ensEMBL2id_go_ECM$external_gene_name)
ECM_genes_all <- c(ECM_genes, unlist(Hsa04974_genes))
selected_genes <- resSig.ann[resSig.ann$external_gene_name %in% ECM_genes_all,]
#saveRDS(selected_genes$external_gene_name, "ECM_genes.rds")


message("+-------------------------------------------------------------------------------")
message("+             Fig. 5b  WNT genes                                                ")
message("+-------------------------------------------------------------------------------")

ego_terms_go_wnt <- ensEMBL2id_go
ego_terms_go_wnt$wnt <- grepl("Wnt", ego_terms_go_wnt$name_1006)
ego_terms_go_wnt <- subset(ego_terms_go_wnt, ego_terms_go_wnt$wnt == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% ego_terms_go_wnt$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)
#saveRDS(selected_genes2$external_gene_name, "WNT_signaling_resSig.rds")

message("+-------------------------------------------------------------------------------")
message("+             Fig. 5b  WNT genes                                                ")
message("+-------------------------------------------------------------------------------")

cellcycle_genes <- ensEMBL2id_go
cellcycle_genes$cc <- grepl("cell cycle", cellcycle_genes$name_1006)
cellcycle_genes$cc2 <- grepl("cyclin", cellcycle_genes$name_1006)
cellcycle_genes <- subset(cellcycle_genes, cellcycle_genes$cc == TRUE | cellcycle_genes$cc2 == TRUE)
selected_genes2 <- resSig.ann[resSig.ann$external_gene_name %in% cellcycle_genes$external_gene_name,]
selected_genes2 <- subset(selected_genes2, abs(selected_genes2$log2FoldChange) > 0.6)
#saveRDS(selected_genes2$external_gene_name, "cell_cycle_genes_resSig.rds")



sessionInfo()
