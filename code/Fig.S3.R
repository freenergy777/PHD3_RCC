suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(ggplot2)
  library(survival)
  library(tibble)
  library(openxlsx)
  library(DESeq2)
  library(R.utils)
  library(org.Hs.eg.db)
  library(pheatmap)
  library(tidyverse)
  library(pheatmap)
  library(ggpubr)
  library(edgeR)
  library(biomaRt)
})

KIRC = GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling", # for expression data
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts", # for raw count
  data.type = "Gene Expression Quantification",
  access = "open")
GDCdownload(KIRC)

KIRC_assay <- GDCprepare(KIRC) # transform the downloaded data into summerizedExperimental object or data frame
KIRC_raw_count <- assay(KIRC_assay,1) # select raw integer count data
keep_KIRC <- rowSums(cpm(KIRC_raw_count)>1)>=2 # keep genes with CPM values greater than one in at least two samples 
cnt_KIRC <- as.matrix(KIRC_raw_count[keep_KIRC,])
cnt_KIRC <- as.data.frame(cnt_KIRC)
cnt_KIRC <- cnt_KIRC[,-c(47,114)] # rule-out inappropriate case
cnt_KIRC <- as.matrix(cnt_KIRC)

group <- substring(colnames(cnt_KIRC),14,15)
group[group%in%'01'] <- 'Tumor'
group[group%in%'11'] <- 'Normal'
group <- factor(group, levels = c('Normal', 'Tumor'))
design <- model.matrix(~group)
dim(cnt_KIRC) # dimesion check for cnt_KIRC
dim(design) # dimesion check for design matrix
rownames(design) <- colnames(cnt_KIRC)

#################################################################
KIRC_TMM <- read_excel("KIRC_TMM.xlsx")
KIRC_TMM_BC_domain <- as.data.frame(KIRC_TMM)
KIRC_TMM_BC_domain <- column_to_rownames(KIRC_TMM_BC_domain, var = 'Label')
KIRC_TMM_BC_domain_t <- t(KIRC_TMM_BC_domain)

annotation_col <- data.frame(
  CellType = c(rep('Solid Tissue Normal', 73), rep('Primary solid Tumor', 540)),
  row.names = colnames(KIRC_TMM_BC_domain_t))

plot <- pheatmap(KIRC_TMM_BC_domain_t,
                 scale = 'column',
                 color = colorRampPalette(c('steelblue3','white','tomato3'))(100),
                 main = 'TMM data retrieved from patient RNA-seq data',
                 show_colnames = F,
                 cluster_cols = F,
                 annotation_col = annotation_col)
tiff(filename = "heatmap_S3f_1.tiff", units= 'cm', width = 15, height = 6, res = 300)
print(plot)
dev.off()

#################################################################
# TMM normalization
dge <- DGEList(counts = cnt_KIRC, group = group)
dge <- calcNormFactors(dge, method = 'TMM')

# perform limma-trend method based differential expression analysis
logdge <- cpm(dge, log = T, prior.count = 3) # calculate the CPM value
fit <- lmFit(logdge, design) # fit linear model to predict the data or infer the relationship between variables
fit <- eBayes(fit, trend = T) # calculate T value, F value and log-odds based on Bayesian
res_limma <-as.data.frame(topTable(fit, n=Inf)) # extract result table
head(res_limma)

res_limma$sig <- ifelse(res_limma$adj.P.Val < 0.05 & res_limma$logFC >= 1, 'up',
                        ifelse(res_limma$adj.P.Val < 0.05 & res_limma$logFC <= -1, 'down', 'NA'))
res_limma$sig <- as.factor(res_limma$sig)
summary(res_limma$sig)
res_limma <- rownames_to_column(res_limma, var = 'ensembl_gene_id')
res_limma$ensembl_gene_id <- gsub('\\..*','',res_limma$ensembl_gene_id) 

listEnsembl()  
ensembl <- useEnsembl(biomart = 'genes')  
datasets <- listDatasets(ensembl)
ensembl.con <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
attr <- listAttributes(ensembl.con)
filters <- listFilters(ensembl.con)
gene_ID <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = res_limma$ensembl_gene_id,
                 mart = ensembl.con)
res_limma_with_gene_ID <- merge(gene_ID, res_limma, by = 'ensembl_gene_id')

rows_to_move <- c(5942,6882,27253,1379,6267,1321,28833,4140,13087,12830)
res_limma_with_gene_ID_egln <-res_limma_with_gene_ID[rows_to_move,]
res_limma_with_gene_ID_remain <- res_limma_with_gene_ID[!res_limma_with_gene_ID$external_gene_name %in% 
                                                          c('EGLN3', 'EGLN2','EGLN1','PCCA', 'PCCB', 'MCCC1','PC','ACACA','ACACB'),]
res_limma_with_gene_ID_final <- rbind(res_limma_with_gene_ID_remain, res_limma_with_gene_ID_egln)
rownames(res_limma_with_gene_ID_final) <- NULL

keyvals.colour <- ifelse(is.na(res_limma_with_gene_ID_final$external_gene_name), 'gray100',
                         ifelse(res_limma_with_gene_ID_final$external_gene_name %in% c('EGLN3', 'EGLN2'), 'tomato3', 
                                ifelse(res_limma_with_gene_ID_final$external_gene_name %in% c('PCCA', 'PCCB', 'MCCC1','PC', 'ACACB'), 'steelblue3','gray100')))

keyvals.colour[is.na(keyvals.colour)] <- 'gray100' <- 'NA'
names(keyvals.colour)[keyvals.colour == 'tomato3'] <- 'up-regulated gene expression'
names(keyvals.colour)[keyvals.colour == 'steelblue3'] <- 'down-regulated gene expression'

p <- EnhancedVolcano(res_limma_with_gene_ID_final,
                     lab = NA,
                     x = 'logFC',
                     y = 'adj.P.Val', 
                     title = 'Differenrial mRNA expression',
                     subtitle = 'bulk RNA-seq',
                     pCutoff = 1e-2, 
                     FCcutoff = 1.0,
                     colAlpha = 0.8,
                     labSize = 1,
                     colCustom = keyvals.colour,
                     cutoffLineType = 'twodash',
                     legendPosition = 'right')

p <- p + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank()) 
p <- p + scale_x_continuous(
  limits = c(-15, 15),
  breaks = c(-10, -5, -1, 0, 1, 5, 10))
print(p)
tiff(filename = "volcanoplot_S3E.tiff", units= 'cm', width = 30, height = 20, res = 300)
print(p)
dev.off()
