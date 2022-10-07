library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("GOSemSim")
library("rtracklayer")

## Anottation -------------------------
gff <- import.gff("gencode.v41.annotation.gff3")
gff<- as.data.frame(gff@elementMetadata) 
#Filtrar tres colunas do gff
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
#Tirar linhas duplicadas
gff <- dplyr::distinct(gff)

gff$gene_id <- sub("\\..*", "", gff$gene_id)

## Genes -------------------------
load("GDCRNATools_KIRC.RData")

##ceRNA
genes_ceRNA <- nodes[,c(2,3)]
genes_ceRNA <- dplyr::rename(genes_ceRNA, "Hugo_Symbol" = "symbol")
genes_ceRNA <- merge(genes_ceRNA,gff, by.y = "gene_name", by.x = "Hugo_Symbol", all.x = TRUE)
genes_ceRNA <- genes_ceRNA[,-c(3)]

pc_gene <- genes_ceRNA[genes_ceRNA$type == "pc",]
gene <- bitr(pc_gene$Hugo_Symbol, fromType = "SYMBOL", 
             toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

pc_genes <- merge(gene,dePC,by.x = "SYMBOL", by.y = "symbol",all.x = TRUE)

pc_geneList = pc_genes[,5]
names(pc_geneList) = as.character(pc_genes[,2])
pc_geneList = sort(pc_geneList, decreasing = TRUE)

## Removing extra information ----------------
rm(ceOutput,ceOutput2,deALL,deALL_MIR,DEGAll_DESeq2,DEGAll_DESeq2_MIR,deLNC,dePC,
   edges,nodes)

## Enrichment Analysis ---------------
# Gene Ontology BP
de <- names(pc_geneList)
edo <- enrichGO(de,'org.Hs.eg.db', 'ENTREZID', ont = "BP")
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
tplot_GO <- treeplot(edox2)
??treeplot()
plot(tplot_KEGG)

# KEGG
edo_kegg <- enrichKEGG(de, organism = "hsa")
edox_kegg <- setReadable(edo_kegg, 'org.Hs.eg.db', 'ENTREZID')
edox2_kegg <- pairwise_termsim(edox_kegg)
tplot_GO <- treeplot(edox2_kegg)
plot(tplot_GO)

# Reactome
edo_reac <- enrichPathway(de, "human")
edox_reac <- setReadable(edo_reac, 'org.Hs.eg.db', 'ENTREZID')
edox2_reac <- pairwise_termsim(edox_reac)
tplot_reac <- heatplot(edox2_reac)
plot(tplot_reac)

## enricher wikiPathways, MSigDB

cell_marker_data <- read_excel("Cell_marker_Human.xlsx")

## instead of `cellName`, users can use other features (e.g. `cancerType`)
cells <- cell_marker_data %>%
  dplyr::select(cancer_type, GeneID) %>%
  dplyr::mutate(geneID = strsplit(GeneID, ', ')) %>%
  tidyr::unnest()

x <- enricher(de, TERM2GENE = cells)
heatplot(x)
