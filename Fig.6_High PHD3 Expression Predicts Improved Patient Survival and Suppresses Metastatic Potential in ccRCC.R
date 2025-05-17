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
