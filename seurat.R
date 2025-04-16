library(Seurat)
library(future)
library(SeuratDisk)


w0_all <- Read10X(w0_dir)
c3d1_all <- Read10X(c3d1_dir)
pd_all <- Read10X(pd_dir)
w0 <- CreateSeuratObject(w0_all, project = "W0")
c3d1 <- CreateSeuratObject(c3d1_all, project = "C3D1")
pd <- CreateSeuratObject(pd_all, project = "PD")


sc_all <- merge(
  w0,
  y = list(c3d1, pd),
  add.cell.ids = c("W0", "C3D1", "PD"),  # prefix 자동 붙음
  project = "ALL"
)

# ---------------------------------------------------------------
# 0. 라이브러리
library(Seurat)
library(DoubletFinder)
library(ggplot2)

# ---------------------------------------------------------------
# 1. 기본 Seurat 객체 생성
seurat_obj <- CreateSeuratObject(counts = raw_counts, project = "scQC")

# ---------------------------------------------------------------
# 2. 기본 QC 지표 계산
sc_all[["percent.mt"]]   <- PercentageFeatureSet(sc_all, pattern = "^MT-")
sc_all[["percent.ribo"]] <- PercentageFeatureSet(sc_all, pattern = "^RPS|^RPL")

sc_all_f1 <- subset(sc_all,
                 subset = nFeature_RNA > 200 &
                          nFeature_RNA < 5000 &
                          percent.mt < 15 &
                          percent.ribo < 60)


sc_all_f1_sc <- SCTransform(sc_all_f1, verbose = TRUE)

# 3. PCA → UMAP
sc_all_f1_sc <- RunPCA(sc_all_f1_sc)
sc_all_f1_sc <- RunUMAP(sc_all_f1_sc, dims = 1:20)


VlnPlot(sc_all_f1_sc, features = c("nFeature_RNA", "percent.mt"))
thresh <- quantile(sc_all_f1_sc$nCount_RNA, 0.99)

# 필터링
sc_all_f1_sc <- subset(sc_all_f1_sc, subset = nCount_RNA < thresh)
sc_all_f1_sc <- PrepSCTFindMarkers(sc_all_f1_sc)

# 4. 이제부터는 경고 없이 VlnPlot, FeaturePlot 가능
VlnPlot(sc_all, features = c("nFeature_RNA", "percent.mt"))
# ---------------------------------------------------------------
# 3. QC 시각화 (violin plot)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        pt.size = 0.1, ncol = 4)


library(DoubletFinder)
 options(future.globals.maxSize = 2 * 1024^3)




sc_w0 <- subset(sc_all_f1_sc, subset = orig.ident == "W0")
DefaultAssay(sc_w0) <- "SCT"

# sweep + find.pK
sweep.res.list <- paramSweep(sc_w0, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal.pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
optimal.pK


sc_w0 <- FindNeighbors(sc_w0, dims = 1:20)
sc_w0 <- FindClusters(sc_w0, resolution = 0.5) 
# 적당한 pK 고른 후
annotations <- sc_w0$seurat_clusters  # cluster가 있어야 함. 없으면 clustering 돌려야 함
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc_w0@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# doubletFinder 실행
sc_w0 <- doubletFinder(
  sc_w0,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.06,  # 너가 고른 값
  nExp = nExp_poi.adj,
  reuse.pANN = NULL,
  sct = TRUE
)


# W0
Idents(sc_w0) <- sc_w0[["DF.classifications_0.25_0.08_686"]][,1]
doublet_w0 <- WhichCells(sc_w0, idents = "Doublet")
sc_w0_clean <- subset(sc_w0, cells = setdiff(Cells(sc_w0), doublet_w0))

# C3D1
Idents(sc_c3d1) <- sc_c3d1[["DF.classifications_0.25_0.06_813"]][,1]
doublet_c3d1 <- WhichCells(sc_c3d1, idents = "Doublet")
sc_c3d1_clean <- subset(sc_c3d1, cells = setdiff(Cells(sc_c3d1), doublet_c3d1))

# PD (예: pK=0.05, nExp=761인 경우라면)
Idents(sc_pd) <- sc_pd[["DF.classifications_0.25_0.05_761"]][,1]
doublet_pd <- WhichCells(sc_pd, idents = "Doublet")
sc_pd_clean <- subset(sc_pd, cells = setdiff(Cells(sc_pd), doublet_pd))


sc_all_clean <- merge(sc_w0_clean, y = list(sc_c3d1_clean, sc_pd_clean), add.cell.ids = c("W0", "C3D1", "PD"))



sc_all_sc <- SCTransform(sc_all_clean, verbose = FALSE)

# (필요시) 병렬 처리 코어 수 조정
plan("multisession", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

# 1. Reference 불러오기
reference <- readRDS("pbmc_multimodal_2023.rds")
sc_all_sc <- SCTransform(sc_all, verbose = FALSE)

# 4. Anchor 찾기
anchors <- FindTransferAnchors(
  reference = reference,
  query = sc_all_sc,
  normalization.method = "SCT",            # SCTransform 기반 anchor
  reference.reduction = "spca",            # reference에서 제공하는 PCA 기반
  dims = 1:50
)

# 5. MapQuery()로 label transfer 수행
sc_all_sc <- MapQuery(
  anchorset = anchors,
  query = sc_all_sc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",           # coarse cell type
    celltype.l2 = "celltype.l2"            # fine-grained cell type (우리가 주로 쓰는 거!)
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"             # reference UMAP 좌표에 맞춰 시각화
)
###이상치 제거 (doublet)
sc_filtered <- subset(sc_all_sc, subset = predicted.celltype.l2 != "Doublet")

# nUMI 상위 1% 제거 (예시)
thresh <- quantile(sc_filtered$nCount_RNA, 0.99)
sc_filtered <- subset(sc_filtered, subset = nCount_RNA < thresh)


# 6. Plot UMAP
DimPlot(sc_all_sc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE)

# 7. Gene expression dot plot
DotPlot(cd8_tem,
        features = c("GIMAP7", "CXCR3", "TNFRSF1A"),
        group.by = "orig.ident") +
  RotatedAxis()
