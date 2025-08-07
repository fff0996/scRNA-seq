library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(future)

# -----------------------------------------
# 1. Read and Create Seurat Objects
# -----------------------------------------
sample_names <- c("YON21_W0", "YON21_C3D1", "YON21_PD", "YON95_W0", "YON95_C3D1", "YON95_PD")
seurat_list <- list()

for (sample in sample_names) {
  path <- paste0("/path/to/filtered_feature_bc_matrix/", sample)
  mtx <- Read10X(path)
  obj <- CreateSeuratObject(mtx, project = sample, min.cells = 3, min.features = 200)
  obj$sample <- sample
  seurat_list[[sample]] <- obj
}

# -----------------------------------------
# 2. Add QC Metrics
# -----------------------------------------
for (i in 1:length(seurat_list)) {
  seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^MT-")
  seurat_list[[i]][["percent.ribo"]] <- PercentageFeatureSet(seurat_list[[i]], pattern = "^RPS|^RPL")
}

# -----------------------------------------
# 3. Filter by QC thresholds
# -----------------------------------------
for (i in 1:length(seurat_list)) {
  seurat_list[[i]] <- subset(seurat_list[[i]],
                             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 &
                                      percent.mt < 15 & percent.ribo < 60)
}

# -----------------------------------------
# 4. Normalize with SCTransform
# -----------------------------------------
seurat_list <- lapply(seurat_list, function(x) SCTransform(x, verbose = FALSE))

features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
library(future)
plan("sequential")  # 병렬 처리 대신 순차 실행
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
#병렬처리 (이렇게 하면 터짐...)
#plan("multisession", workers = 4)  # 예: 4-core 사용
#options(future.globals.maxSize = 100 * 1024^3)
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
                      

combined <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# 8. Set default assay to integrated
DefaultAssay(combined) <- "integrated"

# 9. Run PCA
combined <- RunPCA(combined)

# 10. Run UMAP
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

samples <- c("YON21_W0", "YON21_C3D1", "YON21_PD", "YON95_W0", "YON95_C3D1", "YON95_PD")

# 1. subset each sample from integrated object
sample_list <- lapply(samples, function(s) subset(combined, subset = orig.ident == s))
names(sample_list) <- samples

# 결과 저장할 리스트
clean_list <- list()

# 2. 각 샘플에 대해 반복 실행
for (sample_name in samples) {
  cat("Processing:", sample_name, "\n")
  
  sample <- sample_list[[sample_name]]
  DefaultAssay(sample) <- "SCT"
  
  # [1] Clustering 먼저
  sample <- FindNeighbors(sample, dims = 1:20, graph.name = "SCT_snn", verbose = FALSE)
sample <- FindClusters(sample, graph.name = "SCT_snn", resolution = 0.5, verbose = FALSE)
  
  # [2] Sweep to find optimal pK
  sweep.res.list <- paramSweep(sample, PCs = 1:20, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal.pK <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
  cat("Optimal pK for", sample_name, ":", optimal.pK, "\n")
  
  # [3] Estimate expected doublets
  annotations <- sample$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * nrow(sample@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # [4] Run DoubletFinder
  sample <- doubletFinder(
    sample,
    PCs = 1:20,
    pN = 0.25,
    pK = optimal.pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = TRUE
  )

 df_col <- grep("DF.classifications", colnames(sample@meta.data), value = TRUE)
  df_col <- tail(df_col, 1)  # 가장 마지막 것만 사용

  Idents(sample) <- sample[[df_col]][, 1]
  doublet_cells <- WhichCells(sample, idents = "Doublet")
  sample_clean <- subset(sample, cells = setdiff(Cells(sample), doublet_cells))
  
  clean_list[[sample_name]] <- sample_clean
}

DimPlot(combined, group.by = "seurat_clusters", label = TRUE)
