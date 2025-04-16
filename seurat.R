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

# 6. Plot UMAP
DimPlot(sc_all_sc, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE)

# 7. Gene expression dot plot
DotPlot(cd8_tem,
        features = c("GIMAP7", "CXCR3", "TNFRSF1A"),
        group.by = "orig.ident") +
  RotatedAxis()
