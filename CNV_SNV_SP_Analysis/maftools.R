############# KIRC_maftools ##################
# MAF files from Kidney Renal Clear Cell Carcinoma

# maftools: Summarize, Analyze and Visualize Mutation Anotated Files (MAF) Files
# URL: https://www.bioconductor.org/packages/release/bioc/html/maftools.html
# version 2.4.05

#### Installing packages ####
packages_bioconductor <- c("TCGAbiolinks","maftools","BSgenome.Hsapiens.UCSC.hg38","SummarizedExperiment")
packages_cran <- c("DT", "tidyverse", "stringr", "data.table", "pheatmap","NMF")

#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed from Bioconductor and loaded
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})
package.check <- lapply(packages_cran, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

rm(packages_cran, packages_bioconductor, package.check)

setwd("~/MESTRADO/CNV_ANALYSIS_MAFTOOLS")

## KIRC Maf files ---------------------------
load("kirc_maf_mafclin.RData")

## Reading gistic or CNV -------------------------

all.lesions_95_01 <- read_file("./data/all_lesions.conf_95_01.txt")
amp.genes_95_01 <- read_file("./data/amp_genes.conf_95_01.txt")
del.genes_95_01 <- read_file("./data/del_genes.conf_95_01.txt")
scores.gis_95_01 <- read_file("./data/scores_95_01.gistic")

KIRC_95_01.gistic = readGistic(gisticAllLesionsFile = all.lesions_95_01, gisticAmpGenesFile = amp.genes_95_01,
                            gisticDelGenesFile = del.genes_95_01, gisticScoresFile = scores.gis_95_01, isTCGA = TRUE)

# GISTIC object
KIRC_95_01.gistic

# checking 

getSampleSummary(KIRC_95_01.gistic) # 534 samples (Tumor_Sample_Barcode)
getGeneSummary(KIRC_95_01.gistic) # 1062 genes (hugo)
getCytobandSummary(KIRC_95_01.gistic) # 365 cytobands
write.GisticSummary(gistic = KIRC_95_01.gistic, basename = 'KIRC_95_01.gistic')

# Save gistic data
save(KIRC_95_01.gistic,file = "KIRC_gistic.RData")
## Vizualizing Gistic ---------------------------------

# genome plot
gisticChromPlot(gistic = KIRC_95_01.gistic, markBands = "all")

# bubble plot
gisticBubblePlot(gistic = KIRC_95_01.gistic)

# oncoplot 
dev.new()
gisticOncoPlot(gistic = KIRC_95_01.gistic, clinicalData = kirc.mafclin, 
               clinicalFeatures = 'ajcc_pathologic_stage', sortByAnnotation = TRUE, top = 10)  

gisticOncoPlot(gistic = KIRC_95_01.gistic, top = 10)


## Enrichment Analysis ####
gene_95 <- as.data.frame(KIRC_95_01.gistic@gene.summary)

library("gprofiler2")
library("clusterProfiler")
library("GDCRNATools")
library("org.Hs.eg.db")
library("annotate")

gene2 <- gene_95$Hugo_Symbol
gene2 <- gene2 %>% str_replace_all("-","hifen")
gene2 <- gene2 %>% str_replace_all("_","under")
gene2 <- gene2 %>% str_replace_all("\\.","point")
gene2 <- gsub("[^[:alnum:]]","",as.character(gene2))
gene2 <- gene2 %>% str_replace_all("hifen","-")
gene2 <- gene2 %>% str_replace_all("under","_")
gene2 <- gene2 %>% str_replace_all("point","\\.")

gene2 <- as.data.frame(gene2)
colnames(gene2) <- "Hugo_Symbol"

gff <- import.gff("gencode.v41.annotation.gff3")
gff<- as.data.frame(gff@elementMetadata) 
#Filtrar tres colunas do gff
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
#Tirar linhas duplicadas
gff <- distinct(gff)

gff$gene_id <- sub("\\..*", "", gff$gene_id)

mer <- merge(gene2,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)

for(i in 1:length(mer$Hugo_Symbol[])){
  if (is.na(mer[i,"gene_id"])){
    mer[i,"gene_id"] <- mer[i,"Hugo_Symbol"]
  }
}

##gprofiler2
enrich <- gost(mer$gene_id,"hsapiens")
pathways_gp <- enrich[["result"]][["term_name"]]
pathways_gp <- pathways_gp %>% as.data.frame()

gostplot(enrich, capped = TRUE, interactive = TRUE)

## clusterProfiler
enrich2 <- gdcEnrichAnalysis(mer$gene_id)
gdcEnrichPlot(enrich2, type = 'bubble', category = "DO", num.terms = 10) 

DO <- enrich2[enrich2$Category == "DO",]
KEEG <- enrich2[enrich2$Category == "KEGG",]
GO_BP <- enrich2[enrich2$Category == "GO_BP",]
GO_MF <- enrich2[enrich2$Category == "GO_MF",]
GO_CC <- enrich2[enrich2$Category == "GO_CC",]


## Oncogenic Signaling Pathways --------------------
genes <- KIRC_95_01.gistic@data

genes_amp <- genes[genes$Variant_Classification == "Amp",]
genes_un_amp <- unique(genes_amp$Hugo_Symbol)
genes_un_amp <- genes_un_amp %>% str_replace_all("-","hifen")
genes_un_amp <- genes_un_amp %>% str_replace_all("_","under")
genes_un_amp <- genes_un_amp %>% str_replace_all("\\.","point")
genes_un_amp <- gsub("[^[:alnum:]]","",as.character(genes_un_amp))
genes_un_amp <- genes_un_amp %>% str_replace_all("hifen","-")
genes_un_amp <- genes_un_amp %>% str_replace_all("under","_")
genes_un_amp <- genes_un_amp %>% str_replace_all("point","\\.")

genes_del <- genes[genes$Variant_Classification == "Del",]
genes_un_del <- unique(genes_del$Hugo_Symbol)
genes_un_del <- genes_un_del %>% str_replace_all("-","hifen")
genes_un_del <- genes_un_del %>% str_replace_all("_","under")
genes_un_del <- genes_un_del %>% str_replace_all("\\.","point")
genes_un_del <- gsub("[^[:alnum:]]","",as.character(genes_un_del))
genes_un_del <- genes_un_del %>% str_replace_all("hifen","-")
genes_un_del <- genes_un_del %>% str_replace_all("under","_")
genes_un_del <- genes_un_del %>% str_replace_all("point","\\.")

OncogenicPathways(maf = kirc.maf)
PlotOncogenicPathways(maf = kirc.maf, pathways = "RTK-RAS", fullPathway = TRUE)

onco_signa_path <- read_tsv("oncogenic_sig_patwhays.tsv")
table(onco_signa_path$Pathway)

# Genes amplificados dentro da via Cell_Cycle
# Del: CDKN2A
cell_cicle_path <- onco_signa_path[onco_signa_path$Pathway == "Cell_Cycle", ]
intersect(genes_un_amp, cell_cicle_path$Gene)
intersect(genes_un_del, cell_cicle_path$Gene)

# Genes amplificados dentro da via WNT
wnt_path <- onco_signa_path[onco_signa_path$Pathway == "WNT", ]
intersect(genes_un_amp, wnt_path$Gene)
intersect(genes_un_del, wnt_path$Gene)

# Genes amplificados dentro da via Hippo
hippo_path <- onco_signa_path[onco_signa_path$Pathway == "Hippo", ] 
intersect(genes_un_amp, hippo_path$Gene)
intersect(genes_un_del, hippo_path$Gene)

# Genes amplificados dentro da via MYC
myc_path <- onco_signa_path[onco_signa_path$Pathway == "MYC", ]
intersect(genes_un_amp, myc_path$Gene)
intersect(genes_un_del, myc_path$Gene)

# Genes amplificados dentro da via NOTCH
# Del: APH1A, NOTCH2
notch_path <- onco_signa_path[onco_signa_path$Pathway == "NOTCH", ]
intersect(genes_un_amp, notch_path$Gene)
intersect(genes_un_del, notch_path$Gene)

# Genes amplificados dentro da via NRF2
nrf2_path <- onco_signa_path[onco_signa_path$Pathway == "NRF2", ]
intersect(genes_un_amp, nrf2_path$Gene)
intersect(genes_un_del, nrf2_path$Gene)

# Genes amplificados dentro da via PI3K
pi3k_path <- onco_signa_path[onco_signa_path$Pathway == "PI3K", ]
intersect(genes_un_amp, pi3k_path$Gene)
intersect(genes_un_del, pi3k_path$Gene)

# Genes amplificados dentro da via RTK-RAS
# Amp: JAK2
# Del: NRAS
rtk_ras_path <- onco_signa_path[onco_signa_path$Pathway == "RTK-RAS", ]
intersect(genes_un_amp, rtk_ras_path$Gene)
intersect(genes_un_del, rtk_ras_path$Gene)

# Genes amplificados dentro da via TGF-Beta
tgf_beta_path <- onco_signa_path[onco_signa_path$Pathway == "TGF-Beta", ]
intersect(genes_un_amp, tgf_beta_path$Gene)
intersect(genes_un_del, tgf_beta_path$Gene)

# Genes amplificados dentro da via TP53
tp53_path <- onco_signa_path[onco_signa_path$Pathway == "TP53", ]
intersect(genes_un_amp, tp53_path$Gene)
intersect(genes_un_del, tp53_path$Gene)

## mRNAs from ceRNAs inside the oncogenic signaling pathways ####
load("GDCRNATools_KIRC.RData")
mRNA_ceRNA <- nodes[nodes$type == "pc",]
rm(ceOutput2,ceOutput,deALL,deALL_MIR,DEGAll_DESeq2,DEGAll_DESeq2_MIR,deLNC,dePC,
   edges,nodes)

# Genes dentro da via Cell_Cycle
# CCND2
intersect(mRNA_ceRNA$symbol, cell_cicle_path$Gene)

# Genes dentro da via WNT
# WNT5A
intersect(mRNA_ceRNA$symbol, wnt_path$Gene)

# Genes dentro da via Hippo
# CSNK1E
intersect(mRNA_ceRNA$symbol, hippo_path$Gene)

# Genes dentro da via MYC
# MXD3
intersect(mRNA_ceRNA$symbol, myc_path$Gene)

# Genes dentro da via NOTCH
intersect(mRNA_ceRNA$symbol, notch_path$Gene)

# Genes dentro da via NRF2
intersect(mRNA_ceRNA$symbol, nrf2_path$Gene)

# Genes dentro da via PI3K
intersect(mRNA_ceRNA$symbol, pi3k_path$Gene)

# Genes dentro da via RTK-RAS
# FGFR2 INSR
intersect(mRNA_ceRNA$symbol, rtk_ras_path$Gene)

# Genes dentro da via TGF-Beta
intersect(mRNA_ceRNA$symbol, tgf_beta_path$Gene)

# Genes dentro da via TP53
intersect(mRNA_ceRNA$symbol, tp53_path$Gene)

## Amplified and Del genes in the ceRNA ? -------------

gene_amp_inters <- intersect(genes_un_amp,mRNA_ceRNA$symbol) # RSRP1 CBFA2T3 RFLNB
un_amp_gene_ceRNA <- genes_amp[genes_amp$Hugo_Symbol %in% gene_amp_inters,]
unique(un_amp_gene_ceRNA$Cytoband)

gene_del_inters <- intersect(genes_un_del,mRNA_ceRNA$symbol) # RSRP1 ATP1A1
un_del_gene_ceRNA <- genes_del[genes_del$Hugo_Symbol %in% gene_del_inters,]
unique(un_del_gene_ceRNA$Cytoband)


## Druggable Interactions ------------
dgi.gene_up <- drugInteractions(genes = genes_un_amp, drugs = TRUE)
a <- unique(dgi.gene_up$drug_name)
dgi.mRNA_ceRNA <- drugInteractions(genes = mRNA_ceRNA$symbol, drugs = TRUE)
b <- unique(dgi.mRNA_ceRNA$drug_name)

c <- intersect(a,b)

#c_amp <- dgi.gene_up[dgi.gene_up$drug_name == c,]

#for (j in 1:length(c)){
#  for (i in 1:length(dgi.gene_up$drug_name)){
#    if(dgi.gene_up[dgi.gene_up$drug_name[i],] == c[j]){
#      c_amp[i,] <- dgi.gene_up[i,i]
#    }
#  }
#}

## References -----------------

# citation("maftools")
# Mayakonda A, Lin DC, Assenov Y, Plass C, Koeffler HP. 2018. Maftools: efficient and comprehensive analysis of somatic variants in cancer. Genome Resarch PMID: 30341162

# about MAF files: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/
# about download data: https://www.bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/query.html
# about costumize oncoplots: http://bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/oncoplots.html
# about validated signatures: https://cancer.sanger.ac.uk/cosmic/signatures

# sessionInfo()
