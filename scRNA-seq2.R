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

clean_list[["YON21_W0"]]
sample_full <- sample_list[["YON21_W0"]]   
merged_obj <- merge(clean_list[[1]], y = clean_list[2:6], 
                    add.cell.ids = names(clean_list), 
                    project = "merged_clean")

merged_clean <- merge(clean_list[[1]], y = clean_list[2:6], add.cell.ids = names(clean_list), project = "merged")

DefaultAssay(merged_clean) <- "SCT"
merged_clean <- RunPCA(merged_clean, verbose = FALSE)
merged_clean <- FindNeighbors(merged_clean, dims = 1:20, verbose = FALSE)
merged_clean <- FindClusters(merged_clean, resolution = 0.5, verbose = FALSE)
merged_clean <- RunUMAP(merged_clean, dims = 1:20, verbose = FALSE)
merged_clean <- RunPCA(merged_clean, verbose = FALSE)
merged_clean <- FindVariableFeatures(merged_clean)
merged_clean <- RunPCA(merged_clean, verbose = FALSE)
merged_clean <- FindNeighbors(merged_clean, dims = 1:20, verbose = FALSE)
merged_clean <- FindClusters(merged_clean, resolution = 0.5, verbose = FALSE)
merged_clean <- RunUMAP(merged_clean, dims = 1:20, verbose = FALSE)

Idents(merged_clean) <- "seurat_clusters"
DimPlot(merged_clean, reduction = "umap", label = TRUE)

                      


                      
# 1. Reference 불러오기
reference <- readRDS("pbmc_multimodal_2023.rds")    

map_one <- function(qobj, ref){
DefaultAssay(qobj) <- "SCT"
qobj <- tryCatch(JoinLayers(qobj), error = function(e) qobj)  # 안전용

 # 1) 각 객체의 SCT scale.data에 실제로 존재하는 feature 집합
rf <- rownames(ref[["SCT"]]@scale.data)
qf <- rownames(qobj[["SCT"]]@scale.data)
common <- intersect(rf, qf)

# 2) 후보 선정: SelectIntegrationFeatures → 실제 공통 feature로 필터
feats0 <- SelectIntegrationFeatures(object.list = list(ref, qobj),
                                   nfeatures = min(3000, length(common)))
feats  <- intersect(feats0, common)

#(안전장치) 만약 개수가 너무 적으면 VariableFeatures 교차로 보충
if (length(feats) < 2000) {
vf <- intersect(VariableFeatures(ref), VariableFeatures(qobj))
feats <- unique(c(feats, intersect(vf, common)))
}
# (더 안전장치) 그래도 부족하면 그냥 common에서 상위 n개 사용
if (length(feats) < 2000) {
feats <- common[seq_len(min(2000, length(common)))]
}
# 3) PrepSCTIntegration (둘 다 같은 anchor.features로)
ref_prep <- PrepSCTIntegration(list(ref), anchor.features = feats)[[1]]
q_prep   <- PrepSCTIntegration(list(qobj), anchor.features = feats)[[1]]
 # 4) Anchors & Map
anc <- FindTransferAnchors(
 reference = ref_prep,
query     = q_prep,
normalization.method = "SCT",
reference.reduction  = "pca",   # 필요 시 "spca"
dims = 1:50
)
q_mapped <- MapQuery(
anchorset = anc,
query = q_prep,
reference = ref_prep,
refdata = list(
celltype.l1 = "celltype.l1",
celltype.l2 = "celltype.l2"
),
reference.reduction = "pca",
reduction.model = "wnn.umap"
)
return(q_mapped)
}

                 
objs <- SplitObject(merged_clean, split.by = "orig.ident")
mapped_list <- lapply(objs, map_one, ref = reference)



                 library(Seurat)
library(dplyr)
library(ggplot2)

## 0) 안전 세팅 (v5 레이어 정리 + RNA 기본)
DefaultAssay(merged) <- "RNA"
merged <- tryCatch(JoinLayers(merged), error = function(e) merged)

## 1) CD8 TEM만 추출
Idents(merged) <- "predicted.celltype.l2"
cd8tem <- subset(merged, idents = "CD8 TEM")
DefaultAssay(cd8tem) <- "RNA"
cd8tem <- tryCatch(JoinLayers(cd8tem), error = function(e) cd8tem)

## 2) 재클러스터링
cd8tem <- NormalizeData(cd8tem, verbose = FALSE)
cd8tem <- FindVariableFeatures(cd8tem, nfeatures = 3000, verbose = FALSE)
cd8tem <- ScaleData(cd8tem, verbose = FALSE)
cd8tem <- RunPCA(cd8tem, npcs = 50, verbose = FALSE)
cd8tem <- FindNeighbors(cd8tem, dims = 1:30)
cd8tem <- FindClusters(cd8tem, resolution = 0.3)  # 0.2~0.6 조절 가능
cd8tem <- RunUMAP(cd8tem, dims = 1:30)

## 3) 시그니처 패널 정의
genes_exh <- c("PDCD1","CTLA4","LAG3","TIGIT","HAVCR2","TOX","ENTPD1","CXCL13",
               "EOMES","BATF","IRF4","NFATC1","NFATC2","PRDM1","CD160","CD101","CD244","CD38","CD69")
genes_cyto <- c("GZMA","GZMB","GZMH","GZMK","GZMM","GNLY","PRF1","FASLG","IFNG","TNF","IL2","IL2RA")
genes_naive<- c("SELL","CCR7","TCF7","LEF1","IL7R","CD27","CD28","KLF2","BACH2","S1PR1","FOXP1")
keep <- function(gs, obj) intersect(gs, rownames(obj))

## 4) 모듈 점수 계산
cd8tem <- AddModuleScore(cd8tem, list(keep(genes_exh, cd8tem)),  name = "exh_")   # exh_1
cd8tem <- AddModuleScore(cd8tem, list(keep(genes_cyto, cd8tem)), name = "cyto_")  # cyto_1
cd8tem <- AddModuleScore(cd8tem, list(keep(genes_naive,cd8tem)), name = "naive_") # naive_1
#####만약 1개의 cell에서 특정유전자발현 보고 싶다면 
####GetAssayData(obj, assay=assay_use, slot="data")["GZMB", cell_id]


                   
## 5) (권장) 샘플/시점 내 표준화 후 상태 라벨링
# orig.ident이 있으면 그 기준으로 z-score; 없으면 전체 z-score
z <- function(x) as.numeric(scale(x))
if (!is.null(cd8tem$orig.ident)) {
  cd8tem$exh_z   <- ave(cd8tem$exh_1,   cd8tem$orig.ident, FUN = z)
  cd8tem$cyto_z  <- ave(cd8tem$cyto_1,  cd8tem$orig.ident, FUN = z)
  cd8tem$naive_z <- ave(cd8tem$naive_1, cd8tem$orig.ident, FUN = z)
} else {
  cd8tem$exh_z   <- z(cd8tem$exh_1)
  cd8tem$cyto_z  <- z(cd8tem$cyto_1)
  cd8tem$naive_z <- z(cd8tem$naive_1)
}

thr <- 0.10  # 경계 여유 (데이터 보고 0.05~0.2 튜닝)
lab <- rep("CD8_undetermined", ncol(cd8tem))
lab[(cd8tem$cyto_z > cd8tem$exh_z) & (cd8tem$cyto_z > cd8tem$naive_z + thr)] <- "CD8_cytotoxic"
lab[(cd8tem$exh_z  > cd8tem$cyto_z) & (cd8tem$exh_z  > cd8tem$naive_z + thr)] <- "CD8_exhausted"
lab[(cd8tem$naive_z > pmax(cd8tem$cyto_z, cd8tem$exh_z) + thr)]               <- "CD8_naive_like"
cd8tem$CD8_state <- lab

## 6) 플롯
p1 <- DimPlot(cd8tem, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle("CD8 TEM – reclustering")
p2 <- DimPlot(cd8tem, group.by = "CD8_state", label = TRUE, repel = TRUE) +
  ggtitle("CD8 TEM – Naive / Cytotoxic / Exhausted / Undetermined")

print(p1); print(p2)

# 마커 검증 DotPlot
markers <- c("CCR7","TCF7","IL7R","GZMB","PRF1","GNLY","PDCD1","LAG3","TOX")
print(DotPlot(cd8tem, features = markers, group.by = "CD8_state") + RotatedAxis())







#KNN 이웃 뽑기
cell_id <- "YON21_C3D1_AAACCATTCCCGACTC-1_1"
k <- 30

## 0) 모양 확인
stopifnot(inherits(obj,"Seurat"), cell_id %in% colnames(obj))
redn <- intersect(c("pca","umap"), names(obj@reductions)); stopifnot(length(redn)>=1)

## 1) 좌표 가져와서 이웃 k개 계산(자체 kNN)
emb <- Embeddings(obj, redn[1])
stopifnot(cell_id %in% rownames(emb))
X <- emb[, 1:min(20, ncol(emb)), drop=FALSE]
x <- X[cell_id, , drop=FALSE]
d2 <- rowSums((X - matrix(x, nrow(X), ncol(X), byrow=TRUE))^2)
nbrs <- names(sort(d2))[2:(k+1)]  # 
