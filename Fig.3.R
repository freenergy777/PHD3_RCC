suppressPackageStartupMessages({
  library(rhdf5)
  library(zellkonverter)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(monocle3)
  library(tibble)
  library(RColorBrewer)
  library(openxlsx)
  library(SeuratWrappers)
  library(pheatmap)
  library(harmony)
  library(sctransform)
  library(scater)
  library(cepo)
  library(glmGamPoi)
  library(future)
  library(reticulate)
  library(future)
  library(future.apply)
})

f = "RCC_upload_final_raw_counts.h5ad"
h5ls(f)
data <- readH5AD(file = f,
                 verbose=T, layers=F, varm=F,
                 obsm=F, varp=F, obsp=F, uns=F)
counts(data) = data@assays@data$X

# Major level annotation
table(data$broad_type)
# Minor level annotation
table(data$annotation)
# Patients
table(data$patient)
# Tissue
table(data$summaryDescription)

data <- data[, !data$summaryDescription %in% c("Blood", "Fat", "Thrombus")]
data <- data[,data@colData@listData[["percent.mt"]] < 10]
# this brings down the number of cells from 270855 to 187200

# edit gene expression data
data@assays@data$X = NULL 
data = data[!grepl("^RP", rownames(data)),] # delete ribosomal gene data
data = scater::logNormCounts(data)
data$tissue_major = paste0(data$summaryDescription, "__", 
                           data$broad_type)
keep = names(which(table(paste0(data$summaryDescription, "__", data$broad_type)) > 100)) # cell type sorting 
data = data[, data$tissue_major %in% keep]
# this brings down the number of cells from 187200 to 186610
table(data$broad_type)

#################################################################
# pre-processing & normalization
data_seurat = do.call(cbind, lapply(unique(data$broad_type), function(x) {
  subset = data[, data$broad_type == x]
  if (ncol(subset) > 15000) {
    set.seed(1)
    subset = subset[, sample(1:ncol(subset), 15000)]
  }
  return(subset)
}))
# batch effect normalization, non-linear dimensionality reduction
plan("multisession", workers=4)
set.seed(123)
options(future.seed = TRUE)
data_seurat <- Seurat::as.Seurat(data_seurat, counts='counts', data = 'logcounts') %>% # count as raw data; data as log transformed data
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  SCTransform(assay = "originalexp", vars.to.regress = "percent.mt", verbose = T, method 
              ='glmGamPoi', n_genes=3000,
              conserve.memory = TRUE, return.only.var.genes = TRUE) %>%
  RunPCA(npcs=10, approx = TRUE) %>%
  RunHarmony(group.by.vars = "patient", assay.use = "SCT") %>%
  FindNeighbors(dims = 1:10, reduction = 'harmony') %>% # for SNN
  RunUMAP(dims = 1:10, umap.method = "uwot", reduction = 'harmony') %>%
  FindClusters(random.seed = 123, algorithm = 4) # algorithm 4 for leiden cluster
saveRDS(data_seurat, file = 'D:/2025.04/monocle3_practice/20250420_3k/20250420_3k_SCTransform_suera.rds')

DimPlot(data_seurat, reduction = "umap", group.by = "broad_type") 
DimPlot(data_seurat, reduction = "umap", group.by = "summaryDescription")

#################################################################
rd = Embeddings(data_seurat, "umap")
dfumap = as.data.frame(rd)
dfumap$tissue = data_seurat$summaryDescription
dfumap$tissue_major = data_seurat$tissue_major
dfumap$broad = data_seurat$broad_type
dfumap$egln3 = data_seurat@assays$originalexp@data["EGLN3",]
dfumap$EGLN3_expression = ifelse(dfumap$egln3>0, 'Positive', 'Negative') 
dfumap$pc_sig = colMeans(data_seurat@assays$originalexp@data[PC_genes,])
dfumap$pc = data_seurat@assays$originalexp@data["PC",]
dfumap$coexp = dfumap$egln3 > 0 & dfumap$pc_sig > 0

dfumap_RCCs <- dfumap %>%
  filter(tissue=='Tumour'& tissue_major=='Tumour__RCC'& broad=='RCC')
dfumap_RCCs_filtered <- dfumap %>%
  filter(tissue=='Tumour'& tissue_major=='Tumour__RCC'& broad=='RCC'& pc_sig>0)
dfumap_elgn3_positive_RCCs <- dfumap %>%
  filter(tissue=='Tumour'& tissue_major=='Tumour__RCC'& broad=='RCC'& EGLN3_expression=='Positive')
dfumap_elgn3_negative_RCCs <- dfumap %>%
  filter(tissue=='Tumour'& tissue_major=='Tumour__RCC'& broad=='RCC'& EGLN3_expression=='Negative')

dfumap_RCCs$EGLN3_expression <- factor(dfumap_RCCs$EGLN3_expression, levels = c('Negative', 'Positive'))
dfumap_RCCs_filtered$EGLN3_expression <- factor(dfumap_RCCs_filtered$EGLN3_expression, levels = c('Negative', 'Positive'))

p1 <- ggboxplot(dfumap_RCCs_filtered, x='EGLN3_expression', y='pc_sig', width = 0.4, outlier.shape = 1, size = 0.4, outlier.size = 0.5,
               color = 'black', notch = T, xlab = 'EGLN3 expression', ylab = 'PC related protein \n transcriptomic signatuare',
               fill = 'EGLN3_expression') +
  stat_compare_means(method = "t.test", label = "p.signif", comparisons = list(c('Negative', 'Positive')), bracket.size = 0.3) 
p1 <- p1 + scale_fill_manual(values = c('Negative'='skyblue3', 'Positive'='red1'))
ggarrange(p1)

# Gene expression
dfqc = do.call(rbind, lapply(keep, function(x) {
  idx = data$tissue_major == x
  data_subset = logcounts(data[c(goi, PC_genes), idx])
  
  df = data.frame(
    label = x,
    gene = rownames(data_subset),
    mean = rowMeans(data_subset),
    prop = rowSums(data_subset > 0)/ncol(data_subset)
  )
  return(df)
}))

dfqc$variable = as.character(dfqc$label)
dfqc$tissue = sapply(strsplit(dfqc$label, "__"), "[[", 1)
dfqc$major = sapply(strsplit(dfqc$label, "__"), "[[", 2)
dfqc$type = dfqc$tissue
dfqc$type[grepl("Normal", dfqc$type)] = "Normal"
dfqc$type[grepl("Tumour.normal", dfqc$type)] = "Adjacent"
dfqc$type[grepl("Tumour", dfqc$type)] = "Tumour"
dfqc$tissue = gsub("Normal.adrenal", "Adrenal", dfqc$tissue)
dfqc$tissue = gsub("Normal.kidney", "Kidney", dfqc$tissue)

dfqc$type = factor(dfqc$type, levels = c("Normal", "Adjacent","Tumour", "Metastasis"))
dfqc$type = factor(dfqc$type, levels = c("Normal", "Adjacent","Tumour", "Metastasis"))

dfqc$major = factor(dfqc$major, levels = c("RCC", "Epi_PT", "Epi_non-PT", "EC", "Fibro", "T-cell", "B-cell",
                                           "Myeloid", "NK", "Plasma", "Mast","pDC"))
p1 = dfqc %>%
  filter(gene != "EGLN3") %>%
  ggplot(aes(x = type, y = gene, col = mean, size = prop)) + 
  geom_point() + scale_color_distiller(palette = "Spectral") + 
  theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~major, ncol = 12)

p2 = dfqc %>%
  filter(gene == "EGLN3") %>%
  ggplot(aes(x = type, y = gene, col = mean, size = prop)) + 
  geom_point() + scale_color_distiller(palette = "Spectral") + 
  theme(axis.text.x = element_text(angle = 90)) + facet_wrap(~major, ncol = 12)

tiff(filename = "dotplot_w_egln3.tiff", units= 'cm', width = 30, height = 10, res = 300)
print(p2)
dev.off()
