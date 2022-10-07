BiocManager::install("CEMiTool")
library("CEMiTool")

## Carregar Dados -------------------------
load("kirc_maf_mafclin.RData")
load("kirc_GDCRNATools.RData")
load("GDCRNATools_KIRC.RData")

## Anottation -------------------------
gff <- import.gff("gencode.v41.annotation.gff3")
gff<- as.data.frame(gff@elementMetadata) 
#Filtrar tres colunas do gff
gff <- gff[,c("gene_id" ,"gene_type", "gene_name")]
#Tirar linhas duplicadas
gff <- dplyr::distinct(gff)

gff$gene_id <- sub("\\..*", "", gff$gene_id)

## Expression Data ----------------------
exp_rna <- as.data.frame(rnaExpr)
exp_rna <- tibble::rownames_to_column(exp_rna)

mer_rna <- merge(exp_rna,gff, by.y = "gene_id", by.x = "rowname", all.x = TRUE)

for(i in 1:length(mer_rna$rowname)){
  if (is.na(mer_rna[i,"gene_name"])){
    mer_rna[i,"gene_name"] <- mer_rna[i,"rowname"]
  }
}

dup <- unique(mer_rna[duplicated(mer_rna$gene_name),"gene_name"])
local_dup <- row.names(mer_rna[mer_rna$gene_name %in% dup,])

for (x in local_dup){
  mer_rna[x, "gene_name"] <- mer_rna[x, "rowname"]
} 

rownames(mer_rna) <- mer_rna$gene_name
mer_rna <- mer_rna[-c(1:5),-c(1,604,605)]

colnames(mer_rna) <- substr(colnames(mer_rna), 1,12)

## Clinical Anottation -------------------------
clinical_anotation <- kirc.mafclin@clinical.data
clinical_anotation <- clinical_anotation[,c(70,3)] 

## Igualar tamanho dos dados de expressão e anotação clínica -------------------
rows <- clinical_anotation$bcr_patient_barcode
rows <- rows[!rows %in% c("TCGA-CJ-4913", "TCGA-AK-3444","TCGA-CJ-4923","TCGA-BP-4988")]

mer_rna2 <- dplyr::select(mer_rna,rows)

clinical_anotation <- clinical_anotation[clinical_anotation$bcr_patient_barcode %in% rows,]
clinical_anotation <-  dplyr::rename(clinical_anotation, SampleName = bcr_patient_barcode )
clinical_anotation <-  dplyr::rename(clinical_anotation, Class = ajcc_pathologic_stage)

clinical_anotation <- clinical_anotation %>%
  mutate(Class = forcats::fct_collapse(Class, S1 = "Stage I", S2= "Stage II", 
                                          S3 = "Stage III", S4 = "Stage IV" ))
clinical_anotation <- as.data.frame(clinical_anotation) 

# Remove Data 
rm(metaMatrix.MIR,metaMatrix.RNA,mirCounts,rnaCounts,rnaExpr)
rm(kirc.maf,kirc.mafclin)
rm(ceOutput,ceOutput2,deALL,deALL_MIR,DEGAll_DESeq2,DEGAll_DESeq2_MIR,deLNC,dePC,edges,nodes)

## Using CEMiTool ----------------------
cem <- cemitool(mer_rna2,clinical_anotation, sample_name_column = "SampleName", class_column = "Class")

cem

save_plots(cem, "all")
