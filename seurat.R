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


#VlnPlot(sc_all_f1_sc, features = c("nFeature_RNA", "percent.mt"))

# 이상치 필터링
thresh <- quantile(sc_all_f1_sc$nCount_RNA, 0.99)
sc_all_f1_sc <- subset(sc_all_f1_sc, subset = nCount_RNA < thresh)
#sc_all_f1_sc <- PrepSCTFindMarkers(sc_all_f1_sc)

# 4. 이제부터는 경고 없이 VlnPlot, FeaturePlot 가능
#VlnPlot(sc_all, features = c("nFeature_RNA", "percent.mt"))
# ---------------------------------------------------------------
# 3. QC 시각화 (violin plot)
#VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
        pt.size = 0.1, ncol = 4)


library(DoubletFinder)
options(future.globals.maxSize = 2 * 1024^3)



# W0
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
  pK = 0.05,  # 너가 고른 값
  nExp = nExp_poi.adj,
  reuse.pANN = NULL,
  sct = TRUE
)

# C3D1
sc_c3d1 <- subset(sc_all_f1_sc, subset = orig.ident == "C3D1")
DefaultAssay(sc_c3d1) <- "SCT"

# sweep + find.pK
sweep.res.list <- paramSweep(sc_c3d1, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal.pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
optimal.pK


sc_c3d1 <- FindNeighbors(sc_c3d1, dims = 1:20)
sc_c3d1 <- FindClusters(sc_c3d1, resolution = 0.5) 
# 적당한 pK 고른 후
annotations <- sc_c3d1$seurat_clusters  # cluster가 있어야 함. 없으면 clustering 돌려야 함
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc_c3d1@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# doubletFinder 실행
sc_c3d1 <- doubletFinder(
  sc_c3d1,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.05,  # 너가 고른 값
  nExp = nExp_poi.adj,
  reuse.pANN = NULL,
  sct = TRUE
)


# PD
sc_pd <- subset(sc_all_f1_sc, subset = orig.ident == "PD")
DefaultAssay(sc_pd) <- "SCT"

# sweep + find.pK
sweep.res.list <- paramSweep(sc_pd, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
optimal.pK <- bcmvn[which.max(bcmvn$BCmetric), "pK"]
optimal.pK


sc_pd <- FindNeighbors(sc_pd, dims = 1:20)
sc_pd <- FindClusters(sc_pd, resolution = 0.5) 
# 적당한 pK 고른 후
annotations <- sc_pd$seurat_clusters  # cluster가 있어야 함. 없으면 clustering 돌려야 함
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075 * nrow(sc_pd@meta.data))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# doubletFinder 실행
sc_pd <- doubletFinder(
  sc_pd,
  PCs = 1:20,
  pN = 0.25,
  pK = 0.09,  # 너가 고른 값
  nExp = nExp_poi.adj,
  reuse.pANN = NULL,
  sct = TRUE
)




# W0
Idents(sc_w0) <- sc_w0[["DF.classifications_0.25_0.05_686"]][,1]
doublet_w0 <- WhichCells(sc_w0, idents = "Doublet")
sc_w0_clean <- subset(sc_w0, cells = setdiff(Cells(sc_w0), doublet_w0))

# C3D1
Idents(sc_c3d1) <- sc_c3d1[["DF.classifications_0.25_0.05_813"]][,1]
doublet_c3d1 <- WhichCells(sc_c3d1, idents = "Doublet")
sc_c3d1_clean <- subset(sc_c3d1, cells = setdiff(Cells(sc_c3d1), doublet_c3d1))

# PD (예: pK=0.05, nExp=761인 경우라면)
Idents(sc_pd) <- sc_pd[["DF.classifications_0.25_0.09_595"]][,1]
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


##TCR clonotype분석
# 1. TCR annotation 불러오기
w0_tcr <- read.csv("W0/vdj_t/filtered_contig_annotations.csv")
c3d1_tcr <- read.csv("C3D1/vdj_t/filtered_contig_annotations.csv")
pd_tcr <- read.csv("PD/vdj_t/filtered_contig_annotations.csv")

# 2. TCR-positive barcode 추출
w0_barcodes <- unique(w0_tcr$barcode)
c3d1_barcodes <- unique(c3d1_tcr$barcode)
pd_barcodes <- unique(pd_tcr$barcode)

tmp <- c3d1_tcr[,c(1,5:10,23,24,27)]
tmp$barcode <- paste("PD",tmp$barcode,sep="_")
pd_tcr_TRA <- tmp[tmp$chain== "TRA",]
pd_tcr_TRB <- tmp[tmp$chain== "TRB",]
tra_unique <- pd_tcr_TRA %>%
  group_by(barcode) %>%
  filter(n() == 1) %>%
  ungroup()
trb_unique <- pd_tcr_TRB %>%
  group_by(barcode) %>%
  filter(n() == 1) %>%
  ungroup()

joined_tcr_strict <- inner_join(tra_unique, trb_unique, by = "barcode")

pd_tcr_m <- as.data.frame(joined_tcr_strict)
pd_tcr_m$cdr3 <- paste(pd_tcr_m$cdr3.x,pd_tcr_m$cdr3.y,sep="_")
pd_tcr_m$cdr3_nt <- paste(pd_tcr_m$cdr3_nt.x,pd_tcr_m$cdr3_nt.y,sep="_")
pd_tcr_m$v_gene <- paste(pd_tcr_m$v_gene.x,pd_tcr_m$v_gene.y,sep="_")
pd_tcr_m$j_gene <- paste(pd_tcr_m$j_gene.x,pd_tcr_m$j_gene.y,sep="_")
pd_tcr_m$clone <- paste(pd_tcr_m$cdr3,pd_tcr_m$cdr3_nt,pd_tcr_m$v_gene,pd_tcr_m$j_gene,sep="__")
head(pd_tcr_m)



set.seed(42)
cd8_tem_cells_sub <- RunPCA(cd8_tem_cells_sub)
cd8_tem_cells_sub <- FindNeighbors(cd8_tem_cells_sub, dims = 1:20)
cd8_tem_cells_sub <- FindClusters(cd8_tem_cells_sub, resolution = 0.2) 
cd8_tem_cells_sub <- RunUMAP(cd8_tem_cells_sub, dims = 1:20)  # 🔥 중요: 새로 돌림!
 DimPlot(cd8_tem_cells_sub,reduction="umap",group.by="seurat_clusters",label=T,repel=T)
# 2. SCE 객체 생성 및 UMAP/cluster 옮기기
sce <- SingleCellExperiment(
  assays = list(logcounts = as.matrix(GetAssayData(cd8_tem_cells_sub, layer = "counts")))
)
reducedDim(sce, "UMAP") <- Embeddings(cd8_tem_cells_sub, "umap")
sce$cluster <- cd8_tem_cells_sub$seurat_clusters # 또는 seurat_clusters

# 3. Slingshot
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP',
 start.clus = "4",
approx_points = 300,
  stretch = 2,
  omega = TRUE)
# 4. Pseudotime 시각화용 dataframe
pt <- slingPseudotime(sce)

df <- data.frame(
  UMAP_1 = reducedDim(sce, "UMAP")[, 1],
  UMAP_2 = reducedDim(sce, "UMAP")[, 2],
  pseudotime = pt[, 1],
  cluster = sce$cluster
)

# 5. 시각화
centroids <- df %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1) +
  scale_color_viridis_c(na.value = "lightgray") +
  geom_text(data = centroids, aes(label = cluster), color = "black", size = 4, fontface = "bold") +
  theme_minimal(base_size = 14) +
  labs(title = "Slingshot Pseudotime with Cluster Labels")



Idents(cd8_tem_cells) <- "seurat_clusters"
markers <- FindAllMarkers(cd8_tem_cells,
                          only.pos = TRUE,       # upregulated만
                          min.pct = 0.25,        # 최소 비율
                          logfc.threshold = 0.25 # 최소 로그 fold change
)

top_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
DoHeatmap(cd8_tem_cells, features = top_markers$gene) +
  theme_minimal()



# 1. TCR마커 유전자 제거 후 marker정리
markers <- FindAllMarkers(cd8_tem_cells_sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_filtered <- markers %>% filter(!grepl("^TRA|^TRB|^TRD|^TRG", gene))

top_markers <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(cd8_tem_cells_sub, features = unique(top_markers$gene)) + NoLegend()




# 2. 기능 스코어 계산 (Effector vs Dysfunction)
cyto_genes <- list(c("GZMB", "PRF1", "IFNG", "GNLY", "GZMA", "NKG7", "KLRG1"))
dys_genes <- list(c("PDCD1", "TOX", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "BATF"))

cd8_tem_clean <- AddModuleScore(cd8_tem_clean, features = cyto_genes, name = "CytotoxicScore")
cd8_tem_clean <- AddModuleScore(cd8_tem_clean, features = dys_genes, name = "DysfunctionScore")


# 3. score 시각화 
VlnPlot(cd8_tem_clean, features = c("CytotoxicScore1", "DysfunctionScore1"), group.by = "seurat_clusters")
FeaturePlot(cd8_tem_clean, features = c("CytotoxicScore1", "DysfunctionScore1"))



# 4. pseudotime 축에서 score 변화 
df <- FetchData(cd8_tem_clean, vars = c("CytotoxicScore1", "DysfunctionScore1", "pseudotime")) %>%
  pivot_longer(cols = c("CytotoxicScore1", "DysfunctionScore1"), names_to = "score", values_to = "value")

ggplot(df, aes(x = pseudotime, y = value, color = score)) +
  geom_point(alpha = 0.2, size = 0.3) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  theme_minimal(base_size = 14)


##분화구 여러개 보기 
library(Slingshot)
library(ggplot2)

# UMAP 좌표 추출
umap_coords <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# 클러스터 정보 붙이기 (색상용)
umap_coords$cluster <- sce$cluster

# Base: UMAP point + 슬링샷 선 겹치기
plot(reducedDim(sce, "UMAP"),
     col = scales::hue_pal()(length(unique(sce$cluster)))[sce$cluster],
     pch = 16, asp = 1, main = "Slingshot Lineages on UMAP")

lines(SlingshotDataSet(sce), lwd = 2, col = "black")



##lineage 별 GZMB, TOX발현도 비교 
# rownames를 cell column으로 추가
df$cell <- rownames(df)
common_cells <- intersect(colnames(cd8_tem_cells_sub), df$cell)
df_sub <- df[df$cell %in% common_cells, ]

df_sub$GZMB <- GetAssayData(cd8_tem_cells_sub, slot = "data")["GZMB", common_cells]
df_sub$TOX <- GetAssayData(cd8_tem_cells_sub, slot = "data")["TOX", common_cells]

ggplot(df_sub, aes(x = pseudotime, y = GZMB)) +
  geom_point(size = 0.7, alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  theme_minimal(base_size = 14) +
  labs(title = "GZMB Expression Along Pseudotime")




df_sub$pseudotime <- pt[, 4]

df_long <- df_sub %>%
  pivot_longer(cols = c("GZMB", "TOX"), names_to = "gene", values_to = "expression")

ggplot(df_long, aes(x = pseudotime, y = expression, color = gene)) +
  geom_point(size = 0.7, alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, size = 1.2) +
  theme_minimal(base_size = 14) +
  labs(title = "GZMB vs TOX Expression Along Pseudotime",
       y = "Expression", x = "Pseudotime") +
  scale_color_manual(values = c("GZMB" = "#1f77b4", "TOX" = "#d62728"))  # 파랑/빨강



##clone dynamics
library(dplyr)
library(ggplot2)
library(tidyr)

# 1. meta 테이블 만들기
meta <- data.frame(
  cell = colnames(cd8_tem_cells_sub),
  clone = cd8_tem_cells_sub$clone,
  sample = cd8_tem_cells_sub$orig.ident
)

# 2. clone별 sample별 cell 수 카운트
clone_counts <- meta %>%
  filter(!is.na(clone)) %>%
  group_by(clone, sample) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = n, values_fill = 0)

# 3. 상태 분류
clone_counts <- clone_counts %>%
  mutate(category = case_when(
    W0 == 0 & C3D1 > 0 ~ "Novel",
    W0 > 0 & C3D1 == 0 ~ "Contracted",
    W0 > 0 & C3D1 > W0 * 2 ~ "Expanded",
    W0 > 0 & C3D1 > 0 ~ "Persistent",
    TRUE ~ NA_character_
  ))

# 4. 비율 계산 (총 비율 기준)
clone_summary <- clone_counts %>%
  count(category) %>%
  mutate(perc = 100 * n / sum(n))

# 5. 플롯
ggplot(clone_summary, aes(x = "w0→c3d1", y = perc, fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Clone Dynamics in CD8 TEM (W0→C3D1)", y = "Percentage of Clones (%)", x = "") +
  scale_fill_manual(values = c(
    "Persistent" = "gray70",
    "Contracted" = "red",
    "Expanded" = "yellow",
    "Novel" = "blue"
  )) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")



#novel clone 뽑기 
# 1. meta table
meta <- data.frame(
  cell = colnames(cd8_tem_cells_sub),
  clone = cd8_tem_cells_sub$clone,
  sample = cd8_tem_cells_sub$orig.ident
)

# 2. 각 샘플별 clone set 추출
clones_W0 <- unique(meta$clone[meta$sample == "W0"])
clones_C3 <- unique(meta$clone[meta$sample == "C3D1"])

# 3. novel 정의: C3D1에만 있고, W0에는 없는 clonotype
novel_clones <- setdiff(clones_C3, clones_W0)
novel_clones <- novel_clones[!is.na(novel_clones)]  # NA 제거

# 4. novel 중 5개 이상 cell에 존재하는 clonotype만 추리기
meta_novel <- meta[meta$clone %in% novel_clones & meta$sample == "C3D1", ]
clone_counts <- table(meta_novel$clone)
valid_clones <- names(clone_counts[clone_counts >= 5])
meta_novel <- meta_novel[meta_novel$clone %in% valid_clones, ]



umap_df <- reducedDim(sce, "UMAP") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")
colnames(umap_df) <- c("cell", "umap_1", "umap_2")

plot_df <- pt_df %>%
  left_join(umap_df, by = "cell") %>%
  left_join(meta_novel, by = "cell") %>%
  filter(!is.na(clone))


ggplot(plot_df, aes(x = umap_1, y = umap_2, color = clone)) +
  geom_point(size = 1.8, alpha = 0.9) +
  facet_wrap(~ lineage) +
  theme_minimal(base_size = 14) +
  labs(title = "Novel Clonotypes (≥5 cells) on Slingshot Lineages", color = "Clone")


ggplot(plot_df, aes(x = umap_1, y = umap_2, color = clone)) +
  geom_point(size = 1.8, alpha = 0.9) +
  facet_wrap(~ lineage) +
  theme_minimal(base_size = 14) +
  labs(title = "Novel Clonotypes (≥5 cells) on Slingshot Lineages", color = "Clone")


all_df <- pt_df %>%
  left_join(umap_df, by = "cell")

novel_cells <- meta_novel$cell
highlight_df <- all_df %>% filter(cell %in% novel_cells)

ggplot(all_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = pseudotime), size = 1, alpha = 0.7) +
  scale_color_viridis_c(option = "plasma", na.value = "lightgray") +
  geom_point(data = highlight_df, shape = 1, size = 2.5, color = "red", stroke = 1) +  # 동그라미 강조
  facet_wrap(~ lineage) +
  theme_minimal(base_size = 14) +
  labs(title = "Slingshot Pseudotime with Novel Clones Highlighted") +
  theme(legend.position = "none")


# 1. novel clone 리스트
novel_clones_c3 <- unique(meta_novel$clone)

# 2. 전체 metadata
meta_all <- data.frame(
  cell = colnames(cd8_tem_cells_sub),
  clone = cd8_tem_cells_sub$clone,
  sample = cd8_tem_cells_sub$orig.ident
)

# 3. novel clone만 필터링
novel_df <- meta_all %>%
  filter(clone %in% novel_clones_c3)

# 4. sample별 cell 수 집계
clone_counts <- novel_df %>%
  group_by(clone, sample) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = n, values_fill = 0)

# 5. fate 분류 (PD에서 5개 이상 셀로 존재해야 Persistent로 인정)
clone_counts <- clone_counts %>%
  mutate(fate = case_when(
    C3D1 > 0 & PD >= 5 ~ "Persistent",
    C3D1 > 0 & PD < 5 ~ "Contracted",
    TRUE ~ "Other"
  ))

# 6. 비율 요약
fate_summary <- clone_counts %>%
  count(fate) %>%
  mutate(perc = 100 * n / sum(n))

# 7. 시각화
ggplot(fate_summary, aes(x = "Novel Clones (C3D1)", y = perc, fill = fate)) +
  geom_bar(stat = "identity", width = 0.7) +
  labs(title = "Fate of Novel Clones from C3D1 (≥5 cells in PD)", y = "Percentage (%)", x = "") +
  scale_fill_manual(values = c("Persistent" = "gray50", "Contracted" = "red")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")



library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# 1. pseudotime, UMAP, cluster 정보 준비
pt_df <- slingPseudotime(sce) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  pivot_longer(-cell, names_to = "lineage", values_to = "pseudotime")

umap_df <- reducedDim(sce, "UMAP") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

meta_df <- data.frame(
  cell = rownames(colData(sce)),
  cluster = sce$cluster,
  clone = cd8_tem_cells_sub$clone[rownames(colData(sce))],
  sample = cd8_tem_cells_sub$orig.ident[rownames(colData(sce))]
)

plot_df <- umap_df %>%
  left_join(pt_df, by = "cell") %>%
  left_join(meta_df, by = "cell")

# 2. novel clone 목록 (기존에 만든 valid_clones 사용)
novel_cells <- plot_df %>%
  filter(clone %in% valid_clones, sample %in% c("C3D1", "PD")) %>%
  mutate(clone = factor(clone))

# 3. pseudotime plot + novel clone cell 위치 강조 (색으로 clone 구분)
ggplot(plot_df, aes(x = UMAP.1, y = UMAP.2, color = pseudotime)) +
  geom_point(size = 1, alpha = 0.3) +
  scale_color_viridis_c(option = "plasma", na.value = "lightgray") +
  geom_point(data = novel_cells, aes(x = UMAP.1, y = UMAP.2, color = clone),
             size = 2.2, alpha = 1, show.legend = TRUE) +
  facet_wrap(~ lineage) +
  theme_minimal(base_size = 14) +
  labs(title = "Novel Clone (5종류) Trajectories on Slingshot Pseudotime") +
  guides(color = guide_legend(override.aes = list(size = 3)))


ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = pseudotime), size = 1, alpha = 0.3) +
  scale_color_viridis_c(option = "plasma", na.value = "lightgray") +
  geom_point(data = novel_cells, aes(x = umap_1, y = umap_2, color = NULL, fill = clone),
             size = 2.3, shape = 21, color = "black", stroke = 0.2, show.legend = TRUE) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(novel_cells$clone)))) +
  facet_wrap(~ lineage) +
  theme_minimal(base_size = 14) +
  labs(title = "Tracking Novel Clones (C3D1→PD) on Slingshot Pseudotime") +
theme(legend.position = "none")

