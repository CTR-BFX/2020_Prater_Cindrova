### RNA-Seq Analysis

Paired-end sequencing was performed on Illumina NextSeq Direct High Output with read lengths of 100 bp. 
QC of sequencing was assessed using FastQC, fastq_screen and Picard, summarised with MultiQC (v0.9dev). 
Reads were trimmed with TrimGalore! and aligned to the human genome (GRCh38) with STAR aligner, with 91.2% reads uniquely mapped and mean of 53.4M paired reads/sample. 
Gene quantification was determined with HTSeq-Counts (v0.6.1p1). 
Additional quality control was performed with rRNA and mtRNA counts script. 
Counts extracted with htseq-counts were used to perform differential gene analysis in R (version 3.5.2) using package DESeq2 (v.1.22.2). 
Read counts were normalised on the estimated size factors. 
Principal component analysis (PCA) was performed on rlog-transformed count data for all genes. 
Gene Ontology and Kegg pathway analysis was performed using clusterProfiler package (v.3.10.1). 
The data matrix for scRNA-seq data were obtained from the Wang lab (GEO accession number GSE89497).  
Heatmaps were generated with ‘ComplexHeatmap' R package (v 1.20.0). Karyoplot was generated with karyoploteR (v1.8.8).
The RNA-sequencing data is accessible through the ArrayExpress series accession number:  E-MTAB-6683 (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6683). 


### DNA Methylation Analysis

DNA Methylation Analysis was performed using Infinium MethylationEPIC array.
Genomic DNA was isolated by QIAamp DNA mini kit (Qiagen, cat. no. 51304) following manufacturer’s instructions. 
Buffer AL (200 μl) was added to the sample, mixed by pulse-vortexing for 15 sec, before incubating at 70oC for 10 min. 
Absolute Ethanol (200 l) was then added to the sample, and mixed by pulse-vortexing for 15 sec before transferring to the QIAamp Mini spin column and centrifuged at 6000 g for 1 min. 
The Mini spin column was washed once with Buffer AW1 (500 μl) following by Buffer AW2 (500 μl) before centrifuging at full speed for 1 min. 
For elution of genomic DNA, DNase-free water (100 μl) was added and incubated for 1 min before centrifuging at 6000 g for 1 min. 
The step repeated one more time with another 100 μl DNase-free water. 
DNA concentration of the samples were quantified by NanoDrop and the DNA quality was checked by resolving in 0.8% agarose gel, in which there was a major band visualized at around 10 kbp without obvious smear below, indicating good quality DNA. 
Genomic DNA oxidative bisulfite (oxBS) conversion was performed using the CEGX TrueMethyl kit (Cambridge Epigenetix / NuGEN,  cat. no. CEGXTMS) and used for microarray-based DNA methylation analysis, performed at GenomeScan (GenomeScan B.V., Leiden, The Netherlands) on the HumanMethylation850 BeadChip (Illumina, Inc., San Diego, CA, U.S.A). 
The EPIC arrau interrogates approximately 865,000 CpG sites representing about 99% of the RefSeq genes.  
The resulting iDAT files were imported and analysed using ChAMP (v2.9.10)1,2. 
Samples were processed filtering for a probe detection p-value <= 0.01, probes with a beadcount <3 in at least 5% of samples, no CpG and known SNPs3 at probe starts, probes aligning to multiple locations,  and QC using the on array control probes. 
In total, 750150 probes on the array passed the filtering and QC steps. 
The BMIQ4 method was used to normalise the two probe types (I and II) present on the array. 
Beta methylation values from the EPIC array range from 0 (unmethylated) to 1 (methylated) and are equivalent of percentage methylation. 

EPIC methylation array data have been deposited in the ArrayExpress database at EMBL-EBI under accession number E-MTAB-XXXX (https://www.ebi.ac.uk/arrayexpress/experiments/E-MTABXXXX). 

