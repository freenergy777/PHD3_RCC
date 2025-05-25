suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(ggplot2)
  library(survival)
  library(tibble)
  library(openxlsx)
  library(DESeq2)
  library(R.utils)
  library(org.Hs.eg.db)
  library(tidyverse)
  library(pheatmap)
  library(ggpubr)
  library(edgeR)
  library(biomaRt)
})

KIRC_clinic <- GDCquery_clinic('TCGA-KIRC')
KIRC_clinic$deceased <- ifelse(KIRC_clinic$vital_status == "Alive", FALSE, TRUE)
KIRC_clinic$overall_survival <- ifelse(KIRC_clinic$vital_status == "Alive", KIRC_clinic$days_to_last_follow_up, KIRC_clinic$days_to_death)

free1_all = GDCquery(
  project = 'TCGA-KIRC',
  data.category = "Transcriptome Profiling", # Exp data
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  access = "open")
GDCdownload(free1_all)

# get counts using GDCprepare as extracting gene expression value from GDC portal
# GDCprepare is function that transform the downloaded data into a summarizedExperiment object or a data frame
free1_data <- GDCprepare(free1_all, summarizedExperiment = TRUE)
free1_matrix <- assay(free1_data, "unstranded") # count data #60,660 entries (gene), 541 total columns (sample)

# extract gene and sample metadata from summarizedExperiment object
# The SummarizedExperiment class is a matrix-like container where rows represent features of interest (e.g. genes, transcripts, exons, etc...) and columns represent samples (with sample data summarized as a DataFrame)
gene_metadata <- as.data.frame(rowData(free1_data)) # gene data # 60660 obs. of 10 variables
coldata <- as.data.frame(colData(free1_data)) # sample data # 541 entries, 73 total columns

# setting up countData object for converting into a format suitable for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = free1_matrix, # gene x sample
                              colData = coldata, # sample x various dummy column
                              design = ~ 1) 

# removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# variance-stabilizing transformation (vst) to set a cut-off standard for the amount of expression
# If the variance is not fixed, the vst is performed.
vsd <- DESeq2::vst(dds, blind = FALSE)
free1_matrix_vst <- assay(vsd)

# get data for SIRT3 gene and add gene metadata information to it using pipe operator
# since the gene_matadata contains the gene name corresponding to the gene_id, move it to the vst matrix
# gather column into key-value pairs via gather function
free1_EGLN3 <- free1_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "EGLN3")

free1_PC <- free1_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "PC")

free1_PDH <- free1_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(., gene_metadata, by = "gene_id") %>%
  filter(gene_name == "PDHA1")

# time, status, strata information acquired → put all information into one data frame
free1_EGLN3$case_id <- gsub('-01.*', '', free1_EGLN3$case_id)
free1_EGLN3 <- merge(free1_EGLN3, KIRC_clinic, by.x = 'case_id', by.y = 'submitter_id')
free1_EGLN3_filterd <- subset(free1_EGLN3, ajcc_pathologic_m !='MX')
free1_EGLN3$ajcc_pathologic_stage <- ifelse(free1_EGLN3$ajcc_pathologic_stage %in% c('Stage I', 'Stage II'), 'early', 'late')

free1_PC$case_id <- gsub('-01.*', '', free1_PC$case_id)
free1_PC <- merge(free1_PC, KIRC_clinic, by.x = 'case_id', by.y = 'submitter_id')
free1_PC_filterd <- subset(free1_PC, ajcc_pathologic_m !='MX')
free1_PC$ajcc_pathologic_stage <- ifelse(free1_PC$ajcc_pathologic_stage %in% c('Stage I', 'Stage II'), 'early', 'late')

free1_PDH$case_id <- gsub('-01.*', '', free1_PDH$case_id)
free1_PDH <- merge(free1_PDH, KIRC_clinic, by.x = 'case_id', by.y = 'submitter_id')
free1_PDH_filterd <- subset(free1_PDH, ajcc_pathologic_m !='MX')
free1_PDH$ajcc_pathologic_stage <- ifelse(free1_PDH$ajcc_pathologic_stage %in% c('Stage I', 'Stage II'), 'early', 'late')

# boxplot
free1_EGLN3_filterd$ajcc_pathologic_m <- factor(free1_EGLN3_filterd$ajcc_pathologic_m, levels = c('M0', 'M1'))
free1_PC_filterd$ajcc_pathologic_m <- factor(free1_PC_filterd$ajcc_pathologic_m, levels = c('M0', 'M1'))
free1_PDH_filterd$ajcc_pathologic_m <- factor(free1_PDH_filterd$ajcc_pathologic_m, levels = c('M0', 'M1'))

q <- ggboxplot(free1_EGLN3_filterd, x='ajcc_pathologic_m',y='counts', outlier.shape = 1, outlier.size = 0.5, 
               ylab = 'VST-normalized expression values', width = 0.4,
               fill = 'ajcc_pathologic_m')
q1 <- q + scale_fill_manual(values = c('M0'='#E7B800', 'M1'='#2E9FDF'))
q2 <- q1 + stat_compare_means(method = 'wilcox.test', label.x = 1.4)
ggarrange(q2)

q <- ggboxplot(free1_PC_filterd, x='ajcc_pathologic_m',y='counts', outlier.shape = 1, outlier.size = 0.5, 
               ylab = 'VST-normalized expression values', width = 0.4,
               fill = 'ajcc_pathologic_m')
q1 <- q + scale_fill_manual(values = c('M0'='#E7B800', 'M1'='#2E9FDF'))
q2 <- q1 + stat_compare_means(method = 'wilcox.test', label.x = 1.4)
ggarrange(q2)

q <- ggboxplot(free1_PDH_filterd, x='ajcc_pathologic_m',y='counts', outlier.shape = 1, outlier.size = 0.5, 
               ylab = 'VST-normalized expression values', width = 0.4,
               fill = 'ajcc_pathologic_m')
q1 <- q + scale_fill_manual(values = c('M0'='#E7B800', 'M1'='#2E9FDF'))
q2 <- q1 + stat_compare_means(method = 'wilcox.test', label.x = 1.4)
ggarrange(q2)

# fitting survival curve
median_value <- median(free1_EGLN3_filterd$counts)
quantile_90 <- quantile(free1_EGLN3_filterd$counts, 0.9)
quantile_80 <- quantile(free1_EGLN3_filterd$counts, 0.8)
quantile_75 <- quantile(free1_EGLN3_filterd$counts, 0.75)
quantile_70 <- quantile(free1_EGLN3_filterd$counts, 0.70)
mean_value <- mean(free1_EGLN3_filterd$counts)

free1_EGLN3_filterd$strata <- ifelse(free1_EGLN3_filterd$counts >= median_value, "HIGH", "LOW")
free1_EGLN3_filterd$strata <- ifelse(free1_EGLN3_filterd$counts >= quantile_90, "HIGH", "LOW")
free1_EGLN3_filterd$strata <- ifelse(free1_EGLN3_filterd$counts >= quantile_80, "HIGH", "LOW")
free1_EGLN3_filterd$strata <- ifelse(free1_EGLN3_filterd$counts >= quantile_75, "HIGH", "LOW")
free1_EGLN3_filterd$strata <- ifelse(free1_EGLN3_filterd$counts >= quantile_70, "HIGH", "LOW")
free1_EGLN3_filterd$case_id <- gsub('-01.*', '', free1_EGLN3_filterd$case_id)
fit1 <- survfit(Surv(overall_survival, deceased) ~ strata, data = free1_EGLN3_filterd)
fit1  

p1 <- ggsurvplot(fit1, 
                 size = 0.2,
                 data = free1_EGLN3_filterd,
                 pval = T, pval.size = 4,
                 linetype = 'strata', surv.median.line = 'hv',
                 ggtheme = theme_bw(),
                 xscale = 'd_y', xlab='Time (years)',
                 break.x.by = 365.25,
                 palette = c('red3', 'black'))
p1$plot <- p1$plot + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) 
print(p1)
