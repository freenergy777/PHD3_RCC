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

PC_genes = c("PC","PCCA","PCCB","MCCC1","MCCC2","ACACA","ACACB")
data$egln3 = logcounts(data)["EGLN3",] # genes are stored in row & cells are stored in column
data$pc = logcounts(data)["PC",]
data$pc_sig = colMeans(logcounts(data)[PC_genes,])
data$coexp = data$egln3 > 0 & data$pc_sig > 0

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
p1 = ggplot(dfumap, aes(umap_1, umap_2, col = broad)) + 
  ggrastr::geom_point_rast(alpha = 1, size = 0.005) + 
  scale_color_manual(values = broad_col) + 
  theme_classic() + ggtitle("Cell type annotations") + 
  theme(aspect.ratio = 1)

p2 = ggplot(dfumap, aes(umap_1, umap_2, col = egln3)) + 
  ggrastr::geom_point_rast(alpha = 1, size = 0.005) + 
  scale_color_viridis_c(option = "C") + 
  theme_classic() + ggtitle("EGLN3 expression") + 
  theme(aspect.ratio = 1)

tiff(filename = "umap_celltype.tiff", units= 'cm', width = 30, height = 20, res = 600)
tiff(filename = "umap_egln3.tiff", units= 'cm', width = 30, height = 20, res = 600)
print(p1)
dev.off()

#################################################################
dfqc[dfqc$gene == "EGLN3" & dfqc$major == "RCC",]
mean(dfqc[dfqc$gene == "EGLN3" & dfqc$major %in% c("Epi_PT", "Epi_non-PT"),"proportion"])

dfqc2 = dfqc %>%
  group_by(major, gene) %>%
  summarise(proportion = mean(proportion))

propviz = dfqc2 %>%
  filter(gene == "EGLN3") %>%
  ggplot(aes(reorder(major, proportion, max), proportion, fill = major)) + 
  geom_col() + ylab("Proportion of cells expressing EGLN3") + 
  scale_fill_manual(values = broad_col) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90))
print(propviz)
