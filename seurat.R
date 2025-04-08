library(Seurat)
library(future)
library(SeuratDisk)

# 💥 병렬 처리 글로벌 사이즈 제한 풀기 (10GB)
options(future.globals.maxSize = 10000 * 1024^2)

# (필요시) 병렬 처리 코어 수 조정
plan("multisession", workers = 4)

# 📌 1. Reference 불러오기
reference <- readRDS("pbmc_multimodal_2023.rds")

# 📌 2. Query 데이터 준비 (예: sc_all or all_tcr 등)
# --> 여기선 all_tcr이 merge된 full object라고 가정

# 🔁 3. SCTransform (normalized counts 기반으로 anchor 찾기용)
all_tcr <- SCTransform(all_tcr, verbose = FALSE)

# 📌 4. Anchor 찾기
anchors <- FindTransferAnchors(
  reference = reference,
  query = all_tcr,
  normalization.method = "SCT",            # SCTransform 기반 anchor
  reference.reduction = "spca",            # reference에서 제공하는 PCA 기반
  dims = 1:50
)

# 📌 5. MapQuery()로 label transfer 수행
all_tcr <- MapQuery(
  anchorset = anchors,
  query = all_tcr,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",           # coarse cell type
    celltype.l2 = "celltype.l2"            # fine-grained cell type (우리가 주로 쓰는 거!)
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"             # reference UMAP 좌표에 맞춰 시각화
)

# ✅ 6. Plot UMAP
DimPlot(all_tcr, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, repel = TRUE)
