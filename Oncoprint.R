library(dplyr)

library(ComplexHeatmap)

library(tools)

cnv_df <- cnv %>%
  filter(gene != "." & !is.na(gene)) %>%
  filter(depth >= 10, probes >= 10, weight >= 0.01) %>%
  mutate(
    Alteration = case_when(
      cn == 0 & log2 < -1.1 ~ "Deep_deletion",
      cn >= 5 & log2 > 1 ~ "Amplification",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Alteration)) %>%
  separate_rows(gene, sep = ",") %>%
  mutate(Sample = "S01-004TR220722") %>%
  select(Gene = gene, Sample, Alteration)

# 우선순위 정의 X
#priority_order <- c(
  "Nonsense" = 1,
  "Frameshift" = 2,
  "Splice" = 3,
  "Missense" = 4,
  "Inframe" = 5,
  "Start Lost" = 6,
  "Other" = 7
  #"Amplification" = 8,
  #"Deep_deletion" = 9
)

#priority_order_ann <- c("stop_gained", "frameshift", "splice", "missense", "inframe", "start_lost")

# VCF 파일 목록 지정
vcf_files <- list.files("your_vcf_dir", pattern = "\\.vcf$", full.names = TRUE)

# 전체 결과 저장
all_result <- list()

for (vcf_file in vcf_files) {
  sample_id <- file_path_sans_ext(basename(vcf_file))
  lines <- readLines(vcf_file)
  body <- lines[!grepl("^#", lines)]
  if (length(body) == 0) next

  df <- read.table(text = body, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(df)[8] <- "INFO"

  # ANN 필드 추출
  ann_raw <- ifelse(grepl("ANN=", df$INFO),
                    gsub(".*ANN=([^;]+).*", "\\1", df$INFO),
                    NA)
  ann_split <- strsplit(ann_raw, ",")

  result <- list()

  for (i in seq_along(ann_split)) {
    anns <- ann_split[[i]]

    if (!is.na(anns[1])) {
      for (ann in anns) {
        x <- strsplit(ann, "\\|")[[1]]
        if (length(x) >= 4) {
          effect <- tolower(x[2])
          gene <- x[4]

          # 우리가 정의한 effect 분류 (impact 사용 안함!)
          class <- case_when(
            grepl("stop_gained", effect) ~ "Nonsense",
            grepl("frameshift", effect) ~ "Frameshift",
            grepl("splice", effect) ~ "Splice",
            grepl("missense", effect) ~ "Missense",
            grepl("inframe", effect) ~ "Inframe",
            grepl("start_lost", effect) ~ "Start Lost",
            TRUE ~ NA_character_
          )

          if (!is.na(class) && gene != "") {
            result[[length(result) + 1]] <- data.frame(
              Gene = gene,
              Sample = sample_id,
              Alteration = class,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (length(result) > 0) {
    all_result[[length(all_result) + 1]] <- do.call(rbind, result)
  }
}



maf_df <- do.call(rbind, all_result)

# 우선순위 정렬
maf_clean <- maf_df %>%
  group_by(Gene, Sample) %>%
  summarise(Alteration = paste(sort(unique(Alteration)), collapse = ";"), .groups = "drop")




col <- c(
  "Missense" = "#8DD3C7",     # mint green
  "Nonsense" = "#FFFFB3",     # yellow
  "Frameshift" = "#BEBADA",   # light purple
  "Splice" = "#FB8072",       # salmon
  "Inframe" = "#80B1D3",      # blue
  "Start Lost" = "#FDB462",   # orange
  "Amplification" = "#E41A1C",# red
  "Deep_deletion" = "#377EB8" # blue
)
> alter_fun <- list(
  background = alter_graphic("rect", fill = "#CCCCCC"),
  Missense = alter_graphic("rect", fill = col["Missense"]),
  Nonsense = alter_graphic("rect", fill = col["Nonsense"]),
  Frameshift = alter_graphic("rect", fill = col["Frameshift"]),
  Splice = alter_graphic("rect", fill = col["Splice"]),
  Inframe = alter_graphic("rect", fill = col["Inframe"]),
  `Start Lost` = alter_graphic("rect", fill = col["Start Lost"]),
  Amplification = alter_graphic("rect", fill = col["Amplification"]),
  Deep_deletion = alter_graphic("rect", fill = col["Deep_deletion"])
)



alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = "#FFFFFF", col = NA))
  },
  Missense = function(x, y, w, h) {
    grid.rect(x, y + h*0.25, w-unit(2, "pt"), h*0.5,
              gp = gpar(fill = col["Missense"], col = NA))
  },
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y - h*0.25, w-unit(2, "pt"), h*0.5,
              gp = gpar(fill = col["Nonsense"], col = NA))
  },
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h*0.5,
              gp = gpar(fill = col["Frameshift"], col = NA))
  },
  Splice = function(x, y, w, h) {
    grid.circle(x, y, r = unit(2, "mm"),
                gp = gpar(fill = col["Splice"], col = NA))
  },
  Inframe = function(x, y, w, h) {
    grid.polygon(
      x = unit(c(x - w*0.4, x + w*0.4, x), "native"),
      y = unit(c(y - h*0.3, y - h*0.3, y + h*0.4), "native"),
      gp = gpar(fill = col["Inframe"], col = NA)
    )
  },
  `Start Lost` = function(x, y, w, h) {
    grid.rect(x, y, w*0.3, h*0.9,
              gp = gpar(fill = col["Start Lost"], col = NA))
  }
 # Amplification = function(x, y, w, h) {
 #   grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
   #           gp = gpar(fill = col["Amplification"], col = NA))
  #},
  #Deep_deletion = function(x, y, w, h) {
   # grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"),
              gp = gpar(fill = col["Deep_deletion"], col = NA))
  #}
)
heatmap_legend_param <- list(
  title = "Alterations",
  at = names(col),
  labels = names(col)
)
library(reshape2)


genes_of_interest <- c("NF1","TP53", "KRAS", "BRAF", "SETBP1", "CDKN2A", "MYC", "PIK3CA"
,"PIK3R2","FGFR1","FGF19","CCND1","FGF4","FGF3","FGF10","RICTOR","SUZ12","NOTCH1","KMT2A","ANKRD26","PIK3C2G","CARD11","KDM5C","ARID1B","MDC1","WT1","ERCC5","CHEK1","MYCN","ALK","EGEFR","MET","CDK6","BRAF","ATM","NRG1")
genes_of_interest <- c("TP53", "KRAS", "BRAF", "SETBP1", "CDKN2A", "MYC", "PIK3CA")

maf_sub <- maf_clean %>%
  filter(Gene %in% genes_of_interest)

onco_df <- maf_sub %>%
  group_by(Gene, Sample) %>%
  summarise(Alteration = paste(sort(unique(Alteration)), collapse = ";"), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Sample, values_from = Alteration, values_fill = "")
mat_top <- as.matrix(onco_df[, -1])
rownames(mat_top) <- onco_df$Gene

oncoPrint(mat_top,
  alter_fun = alter_fun,
  col = col,
  column_title = "Top 20 Mutated Genes",
  heatmap_legend_param = heatmap_legend_param)






alter_fun <- list(
  background = function(x, y, w, h) {
    h_each <- h / 6
    for (i in 1:6) {
      y_pos <- y + h/2 - h_each/2 - (i - 1) * h_each
      grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
                just = "center", gp = gpar(fill = "#CCCCCC", col = NA))
    }
  },

  Missense = function(x, y, w, h) {
    h_each <- h / 6
    y_pos <- y + h/2 - h_each/2 - 0 * h_each
    grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
              just = "center", gp = gpar(fill = col["Missense"], col = NA))
  },

  Nonsense = function(x, y, w, h) {
    h_each <- h / 6
    y_pos <- y + h/2 - h_each/2 - 1 * h_each
    grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
              just = "center", gp = gpar(fill = col["Nonsense"], col = NA))
  },

  Frameshift = function(x, y, w, h) {
    h_each <- h / 6
    y_pos <- y + h/2 - h_each/2 - 2 * h_each
    grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
              just = "center", gp = gpar(fill = col["Frameshift"], col = NA))
  },

  Splice = function(x, y, w, h) {
    h_each <- h / 6
    y_pos <- y + h/2 - h_each/2 - 3 * h_each
    grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
              just = "center", gp = gpar(fill = col["Splice"], col = NA))
  },

  Inframe = function(x, y, w, h) {
    h_each <- h / 6
    y_pos <- y + h/2 - h_each/2 - 4 * h_each
    grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
              just = "center", gp = gpar(fill = col["Inframe"], col = NA))
  },

  `Start Lost` = function(x, y, w, h) {
    h_each <- h / 6
    y_pos <- y + h/2 - h_each/2 - 5 * h_each
    grid.rect(x, y_pos, w - unit(2, "pt"), h_each,
              just = "center", gp = gpar(fill = col["Start Lost"], col = NA))
  }
)



ht <- oncoPrint(
  mat_top_filtered,
  alter_fun = alter_fun,
  alter_fun_is_vectorized = FALSE,
  col = col,
show_column_names = TRUE,
  column_title = "Top Mutated Genes",
  heatmap_legend_param = list(
    title = "Alterations",
    at = names(col),
    labels = names(col),
    legend_gp = gpar(fill = col)
  ),
  name = "oncoprint"  # 내부 heatmap 이름 지정!
)

draw(ht)

decorate_heatmap_body("oncoprint", {
  nr <- nrow(mat_top)
  for (i in 1:(nr - 1)) {
    y_line <- unit(i / nr, "npc")
    grid.lines(x = unit(c(0, 1), "npc"),
               y = y_line,
               gp = gpar(col = "black", lwd = 0.8, lty = 3))  # 점선 여기!
  }
})
