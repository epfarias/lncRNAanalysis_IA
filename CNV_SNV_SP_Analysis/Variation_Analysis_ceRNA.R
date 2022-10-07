############# Variations Analysis ##################

## Installing packages ####
packages_bioconductor <- c("TCGAbiolinks","maftools","BSgenome.Hsapiens.UCSC.hg38","SummarizedExperiment",
                           "ComplexHeatmap","rtracklayer","clusterProfiler", "org.Hs.eg.db","enrichplot")
packages_cran <- c("DT", "tidyverse", "stringr", "data.table", "pheatmap","NMF","gprofiler2")

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

## Load data  ---------------------------
load("kirc_maf_mafclin.RData")
load("GDCRNATools_KIRC.RData")
load("KIRC_gistic.RData")

## Anottation -------------------------
gff <- import.gff("gencode.v41.annotation.gff3")
gff<- as.data.frame(gff@elementMetadata) 
#Filtrar tres colunas do gff
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
#Tirar linhas duplicadas
gff <- dplyr::distinct(gff)

gff$gene_id <- sub("\\..*", "", gff$gene_id)

## Obtaning the unique genes in each dataset -------------------------
# Oncogenic Signaling Pathway
genes_osp <- readr::read_tsv("oncogenic_sig_patwhays.tsv")
genes_osp <- dplyr::rename(genes_osp, "Hugo_Symbol" = "Gene")
genes_osp <- merge(genes_osp,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)
genes_osp <- dplyr::distinct(genes_osp)

genes_vhl_path <- c("EPAS1", "EGLN2", "VHL", "RBX1", "ELOC", "ELOB","CUL2",
                    "ARNT","SLC2A1", "VEGFA", "TGFB1", "PDGFB", "TGFA","CREBBP")

genes_vhl_path <- gff[gff$gene_name %in% genes_vhl_path,]

genes_vhl_path$Pathway <- "VHL_path"

genes_vhl_path$OG_TSG <- "NA"

genes_vhl_path <- dplyr::rename(genes_vhl_path, Hugo_Symbol = gene_name)

for(i in 1:length(genes_osp$Hugo_Symbol[])){
  if (is.na(genes_osp[i,"gene_id"])){
    genes_osp[i,"gene_id"] <- genes_osp[i,"Hugo_Symbol"]
  }
}

genes_osp <- rbind(genes_osp,genes_vhl_path)

# MAF
genes_maf <- kirc.maf@data[,c(1,124)]
genes_maf <- merge(genes_maf,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)
genes_maf <- genes_maf[!duplicated(genes_maf$Hugo_Symbol),]
genes_maf <- merge(genes_maf,kirc.maf@gene.summary, by.y = "Hugo_Symbol", by.x = "Hugo_Symbol", all.x = TRUE) 
genes_maf <- genes_maf[,c(1:4,14)]
qnt_pct <- length(kirc.maf@clinical.data[["Tumor_Sample_Barcode"]])
genes_maf$alteration_percentual <- (genes_maf$total/qnt_pct)*100

for(i in 1:length(genes_maf$Hugo_Symbol[])){
  if (is.na(genes_maf[i,"gene_id"])){
    genes_maf[i,"gene_id"] <- genes_maf[i,"Hugo_Symbol"]
  }
}

genes_maf <- genes_maf[genes_maf$alteration_percentual >= 1,]


# DEA
genes_DEA <- deALL[,c(1,7)]
genes_DEA <- genes_DEA[!duplicated(genes_DEA$symbol),]
genes_DEA <- dplyr::rename(genes_DEA, "Hugo_Symbol" = "symbol")
genes_DEA$VARIANT_CLASS <- "DEA"
genes_DEA <- rownames_to_column(genes_DEA)
genes_DEA <- dplyr::rename(genes_DEA, "gene_id" = "rowname")
genes_DEA <- merge(genes_DEA,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)
genes_DEA <- genes_DEA[,-c(3,5)]
genes_DEA <- dplyr::rename(genes_DEA, "gene_id" = "gene_id.x")

for(i in 1:length(genes_DEA$Hugo_Symbol[])){
  if (is.na(genes_DEA[i,"gene_id"])){
    genes_DEA[i,"gene_id"] <- genes_DEA[i,"Hugo_Symbol"]
  }
}

# Gistic
genes_GISTIC <- KIRC_95_01.gistic@data[,c(1,3,4)]
genes_GISTIC$Hugo_Symbol <- gsub("\\[|\\]", "", genes_GISTIC$Hugo_Symbol)
genes_GISTIC <- dplyr::rename(genes_GISTIC, "VARIANT_CLASS" = "Variant_Type")
genes_GISTIC <- genes_GISTIC[!duplicated(genes_GISTIC$Hugo_Symbol),]
genes_GISTIC <- merge(genes_GISTIC,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)
genes_GISTIC <- genes_GISTIC[,-3]

for(i in 1:length(genes_GISTIC$Hugo_Symbol[])){
  if (is.na(genes_GISTIC[i,"gene_id"])){
    genes_GISTIC[i,"gene_id"] <- genes_GISTIC[i,"Hugo_Symbol"]
  }
}

# ceRNA Network
genes_ceRNA <- nodes[,c(2,3)]
genes_ceRNA <- dplyr::rename(genes_ceRNA, "Hugo_Symbol" = "symbol")
genes_ceRNA <- merge(genes_ceRNA,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)

for(i in 1:length(genes_ceRNA$Hugo_Symbol[])){
  if (is.na(genes_ceRNA[i,"gene_id"])){
    genes_ceRNA[i,"gene_id"] <- genes_ceRNA[i,"Hugo_Symbol"]
  }
}

## Constructing enrichment analysis data
pc_gene <- genes_ceRNA[genes_ceRNA$type == "pc",]
gene <- bitr(pc_gene$Hugo_Symbol, fromType = "SYMBOL", 
             toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

pc_genes <- merge(gene,dePC,by.x = "SYMBOL", by.y = "symbol",all.x = TRUE)

pc_geneList = pc_genes[,5]
names(pc_geneList) = as.character(pc_genes[,2])
pc_geneList = sort(pc_geneList, decreasing = TRUE)

genes_kegg <- bitr_kegg(pc_genes$ENTREZID, fromType = "ncbi-geneid", 
                        toType = "kegg", organism = "hsa")
pc_genes_kegg <- merge(gene,dePC,by.x = "SYMBOL", by.y = "symbol",all.x = TRUE)
pc_geneslist_kegg = pc_genes_kegg[,5]
names(pc_geneslist_kegg) = as.character(pc_genes_kegg[,2])
pc_geneslist_kegg = sort(pc_geneslist_kegg,decreasing = TRUE)

# Removing extra data
rm(ceOutput,ceOutput2, deALL, deALL_MIR, DEGAll_DESeq2, DEGAll_DESeq2_MIR,deLNC, dePC,
   edges, KIRC_95_01.gistic,kirc.maf,kirc.mafclin, nodes,gff)

## Using UpSetplot -------------------
lt <-  list(GISTIC = genes_GISTIC$gene_id,
            MAF = genes_maf$gene_id,
            DEA = genes_DEA$gene_id,
            OSP = genes_osp$gene_id,
            CERNA = genes_ceRNA$gene_id)

lt <- list_to_matrix(lt)

m1 <- make_comb_mat(lt, mode = "intersect")

UpSet(m1, top_annotation = upset_top_annotation(m1, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m1, add_numbers = TRUE))

m1

MAF_DEA_GISTIC_CERNA <- extract_comb(m1, "11101")
MAF_DEA_GISTIC_CERNA2 <- genes_maf[genes_maf$gene_id %in% MAF_DEA_GISTIC_CERNA,]

MAF_DEA_OSP_CERNA <- extract_comb(m1, "01111")
MAF_DEA_OSP_CERNA2 <- genes_maf[genes_maf$gene_id %in% MAF_DEA_OSP_CERNA,]

MAF_DEA_CERNA <- extract_comb(m1, "01101")
MAF_DEA_CERNA2 <- genes_maf[genes_maf$gene_id %in% MAF_DEA_CERNA,]

MAF_GISTIC_CERNA <- extract_comb(m1,"11001")
MAF_GISTIC_CERNA2 <- genes_maf[genes_maf$gene_id %in% MAF_GISTIC_CERNA,]

MAF_OSP_CERNA <- extract_comb(m1, "01011")
MAF_OSP_CERNA2 <- genes_maf[genes_maf$gene_id %in% MAF_OSP_CERNA,]

DEA_GISTIC_CERNA <- extract_comb(m1,"10101")
DEA_GISTIC_CERNA2 <- genes_DEA[genes_DEA$gene_id %in% DEA_GISTIC_CERNA,]

DEA_OSP_CERNA <- extract_comb(m1, "00111")
DEA_OSP_CERNA2 <- genes_DEA[genes_DEA$gene_id %in% DEA_OSP_CERNA,]

MAF_CERNA <- extract_comb(m1, "01001")
MAF_CERNA2 <- genes_maf[genes_maf$gene_id %in% MAF_CERNA,]

DEA_CERNA <- extract_comb(m1, "00101")
DEA_CERNA2 <- genes_DEA[genes_DEA$gene_id %in% DEA_CERNA,]

GISTIC_CERNA <- extract_comb(m1, "10001")
GISTIC_CERNA2 <- genes_GISTIC[genes_GISTIC$gene_id %in% GISTIC_CERNA,]

OSP_CERNA <- extract_comb(m1, "00011")
OSP_CERNA2 <- genes_osp[genes_osp$gene_id %in% OSP_CERNA,]

## genes_ceRNA in other datasets ------------------
genes_ceRNA2 <- genes_maf[genes_maf$Hugo_Symbol %in% genes_ceRNA$Hugo_Symbol, ]
genes_ceRNA3 <- genes_GISTIC[genes_GISTIC$Hugo_Symbol %in% genes_ceRNA$Hugo_Symbol, ]
genes_ceRNA4 <- genes_DEA[genes_DEA$Hugo_Symbol %in% genes_ceRNA$Hugo_Symbol, ]
genes_ceRNA5 <- genes_osp[genes_osp$Hugo_Symbol %in% genes_ceRNA$Hugo_Symbol, ]
uniques <- merge(genes_ceRNA, genes_ceRNA2, by.y = "Hugo_Symbol", by.x = "Hugo_Symbol", all.x = TRUE)
uniques <- uniques[,-c(6,7)]
uniques <- merge(uniques, genes_ceRNA3, by.y = "Hugo_Symbol", by.x = "Hugo_Symbol", all.x = TRUE)
uniques <- uniques[,-c(9,10)]
uniques <- merge(uniques, genes_ceRNA4, by.y = "Hugo_Symbol", by.x = "Hugo_Symbol", all.x = TRUE)
uniques <- uniques[,-c(3,4)]
uniques <- merge(uniques, genes_ceRNA5, by.y = "Hugo_Symbol", by.x = "Hugo_Symbol", all.x = TRUE)
uniques <- uniques[,-c(12,13)]

uniques <- dplyr::rename(uniques, "MAF_DATA" = "VARIANT_CLASS.x")
uniques <- dplyr::rename(uniques, "GISTIC_DATA" = "Variant_Classification")
uniques <- dplyr::rename(uniques, "DEA_DATA" = "VARIANT_CLASS.y")
uniques <- dplyr::rename(uniques, "gene_id" = "gene_id.x")

## Enrichment Analysis Using gProfiler2 -------------------
# genes ceRNA 
gostres_ceRNA <- gost(query = genes_ceRNA$Hugo_Symbol, organism = "hsapiens", sources = "GO:BP")
gostplot(gostres_ceRNA, capped = TRUE, interactive = TRUE)

# genes DEA 
gostres_DEA <- gost(query = genes_DEA$Hugo_Symbol, organism = "hsapiens", sources = "GO:BP")
gostplot(gostres_DEA, capped = TRUE, interactive = TRUE)

# genes maf_dea_gistic
maf_dea_gistic <- extract_comb(m1, "11100")

gostres_maf_dea_gistic <- gost(query = maf_dea_gistic, organism = "hsapiens", sources = "GO:BP")
gostplot(gostres_maf_dea_gistic, capped = TRUE, interactive = TRUE)

gostres_maf<- gost(query = genes_maf$Hugo_Symbol, organism = "hsapiens", sources = "GO:BP")
gostplot(gostres_maf, capped = TRUE, interactive = TRUE)

gostres_gistic <- gost(query = genes_GISTIC$Hugo_Symbol, organism = "hsapiens", sources = "GO:BP")
gostplot(gostres_gistic, capped = TRUE, interactive = TRUE)

## Enrichment Analysis Using clusterProfiler -------------------
library("clusterProfiler")
library("org.Hs.eg.db")

# coding genes from ceRNAs
# Gene Ontology BP
de <- names(pc_geneList)
edo <- enrichGO(de,'org.Hs.eg.db', 'ENTREZID', ont = "BP")
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
tplot_GO <- treeplot(edox2)
plot(tplot_GO)

# KEGG
de_kegg <- names(pc_geneslist_kegg)
edo_kegg <- enrichKEGG(de_kegg, organism = "hsa")
edox_kegg <- setReadable(edo_kegg, 'org.Hs.eg.db', 'ENTREZID')
edox2_kegg <- pairwise_termsim(edox_kegg)
tplot_GO <- treeplot(edox2_kegg)
plot(tplot_GO)

## Searching miRNAs from ceRNA in the LncSEA data base
miRNAs_LncSEA <- read_table2("MicroRNA.csv") # Informação obtida no LncSEA inserindo o 18 lncRNAs da minha ceRNA

mir <- genes_ceRNA[genes_ceRNA$type == "mir",] ## 75 miRNAs

miRNAs_ceRNAs <- miRNAs_LncSEA[ miRNAs_LncSEA$Set %in% mir$Hugo_Symbol,]
table(miRNAs_ceRNAs$Sub_Class)

miRNAs_ceRNAs2 <-miRNAs_ceRNAs[!duplicated(miRNAs_ceRNAs$Set),]
