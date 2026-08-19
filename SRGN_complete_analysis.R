# ============================================================
# SRGN in Canine and Human Osteosarcoma
# Complete Analysis Pipeline
# Author: Ruwini
# Date: 2026
# R version: 4.5.1
# ============================================================

# ── Libraries ─────────────────────────────────────────────────────────
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(survival)
library(survminer)
library(Matrix)
library(readxl)
library(readr)
library(ggrepel)
library(dorothea)
library(viper)
library(patchwork)

# ── Working directory ─────────────────────────────────────────────────
setwd("C:/Users/user/OneDrive/Desktop/GSE252470")


# ============================================================
# LAYER 1: CANINE scRNA-seq ANALYSIS (GSE252470)
# ============================================================

# ── 1.1 Load raw 10x data ─────────────────────────────────────────────
sample_ids <- c("N1_1", "N1_2", "N2_1", "N2_2", "N3", "N4", "N5", "N6")

seurat_list <- lapply(sample_ids, function(sid) {
  counts <- ReadMtx(
    mtx      = file.path(sid, "matrix.mtx"),
    cells    = file.path(sid, "barcodes.tsv"),
    features = file.path(sid, "features.tsv")
  )
  CreateSeuratObject(counts,
                     project      = sid,
                     min.cells    = 3,
                     min.features = 200)
})
names(seurat_list) <- sample_ids

# ── 1.2 Merge all samples ─────────────────────────────────────────────
merged <- merge(seurat_list[[1]],
                y            = seurat_list[-1],
                add.cell.ids = sample_ids)

dim(merged)
# Result: 16,499 genes x 57,685 cells

# Confirm SRGN detected
rownames(merged)[grep("srgn", rownames(merged), ignore.case = TRUE)]
# Result: "SRGN"

# ── 1.3 Quality control ───────────────────────────────────────────────
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

# Visualise QC metrics
VlnPlot(merged,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol     = 3,
        pt.size  = 0)

# Filter: 200–6,000 genes, <50,000 UMI, <25% mitochondrial
merged_filtered <- subset(merged,
                           subset = nFeature_RNA > 200  &
                                    nFeature_RNA < 6000 &
                                    nCount_RNA   < 50000 &
                                    percent.mt   < 25)

dim(merged_filtered)
# Result: 16,499 genes x 47,078 cells

# ── 1.4 Seurat preprocessing pipeline ────────────────────────────────
merged_filtered <- merged_filtered %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)

# ── 1.5 UMAP visualisation ────────────────────────────────────────────
DimPlot(merged_filtered,
        reduction = "umap",
        label     = TRUE,
        pt.size   = 0.3) +
  ggtitle("Canine OSA - all cells")

# SRGN on UMAP
FeaturePlot(merged_filtered,
            features = "SRGN",
            cols     = c("lightgrey", "red"),
            pt.size  = 0.3) +
  ggtitle("SRGN expression across all cells")

# ── 1.6 Cell type marker DotPlot ──────────────────────────────────────
DotPlot(merged_filtered,
        features  = c("PTPRC", "CD3E", "CD79A",
                       "CD68", "LYZ", "CTSK", "ACP5",
                       "FAP", "DCN", "PECAM1",
                       "RUNX2", "SPP1", "SRGN"),
        group.by  = "seurat_clusters",
        cols      = c("lightgrey", "darkred"),
        dot.scale = 6) +
  RotatedAxis() +
  scale_y_discrete(expand = expansion(add = 1.5)) +
  ggtitle("Cell type markers by cluster") +
  xlab("") + ylab("Cluster")

ggsave("dotplot_markers_final.png", width = 12, height = 20, dpi = 300)

# ── 1.7 Cell type annotation ──────────────────────────────────────────
merged_filtered[["cell_type"]] <- dplyr::recode(
  as.character(merged_filtered$seurat_clusters),
  "0"  = "Tumor",      "1"  = "CAF",
  "2"  = "Tumor",      "3"  = "T cell",
  "4"  = "T cell",     "5"  = "Osteoclast",
  "6"  = "Tumor",      "7"  = "Tumor",
  "8"  = "Osteoblast", "9"  = "CAF",
  "10" = "Osteoclast", "11" = "Tumor",
  "12" = "CAF",        "13" = "Osteoblast",
  "14" = "Endothelial","15" = "Tumor",
  "16" = "Tumor",      "17" = "Tumor",
  "18" = "Tumor",      "19" = "Tumor",
  "20" = "Tumor",      "21" = "Tumor",
  "22" = "Macrophage", "23" = "Macrophage",
  "24" = "Tumor",      "25" = "B cell"
)

# Annotated UMAP
DimPlot(merged_filtered,
        group.by = "cell_type",
        label    = TRUE,
        pt.size  = 0.3,
        repel    = TRUE) +
  ggtitle("Canine OSA - cell types")

ggsave("umap_celltypes.png", width = 12, height = 8, dpi = 300)

# ── 1.8 SRGN by cell type ─────────────────────────────────────────────
VlnPlot(merged_filtered,
        features = "SRGN",
        group.by = "cell_type",
        pt.size  = 0,
        cols     = c("#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#FF6B00",
                     "#CC79A7", "#999999")) +
  ggtitle("SRGN expression by cell type in canine OSA") +
  xlab("") + ylab("Normalized expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SRGN_by_celltype.png", width = 10, height = 6, dpi = 300)

# ── 1.9 SRGN enrichment statistics ───────────────────────────────────
merged_filtered <- JoinLayers(merged_filtered)

# Osteoclasts
osteo_cells <- subset(merged_filtered, cell_type == "Osteoclast")
srgn_osteo  <- FetchData(osteo_cells, vars = "SRGN", layer = "counts")

# Tumor cells
tumor_cells <- subset(merged_filtered, cell_type == "Tumor")
srgn_tumor  <- FetchData(tumor_cells, vars = "SRGN", layer = "counts")

cat("Osteoclast SRGN+:", round(mean(srgn_osteo$SRGN > 0) * 100, 2), "%\n")
cat("Tumor SRGN+:",      round(mean(srgn_tumor$SRGN > 0) * 100, 2), "%\n")
# Result: 72.89% osteoclasts vs 17.88% tumor cells

# Chi-square test
contingency <- matrix(
  c(sum(srgn_osteo$SRGN > 0), sum(srgn_osteo$SRGN == 0),
    sum(srgn_tumor$SRGN > 0), sum(srgn_tumor$SRGN == 0)),
  nrow     = 2,
  dimnames = list(c("SRGN+", "SRGN-"), c("Osteoclast", "Tumor"))
)
print(contingency)
chisq.test(contingency)
# Result: chi-sq = 5399.9, p < 2.2e-16

# Save
saveRDS(merged_filtered, "GSE252470_annotated.rds")


# ============================================================
# LAYER 2: BULK RNA-seq CORRELATION (DOG2 / GSE238110)
# ============================================================

# ── 2.1 Load DOG2 bulk RNA-seq ────────────────────────────────────────
dog2 <- read.csv("GSE238110_RawCountFile_combined.csv",
                 row.names   = 1,
                 check.names = FALSE)

dim(dog2)
# Result: 37,952 genes x 198 samples

# Confirm SRGN detected
grep("SRGN", rownames(dog2), ignore.case = TRUE, value = TRUE)
# Result: "SRGN_3"

# ── 2.2 Verify no zero-expression samples ────────────────────────────
cat("Samples with SRGN = 0:", sum(as.numeric(dog2["SRGN_3",]) == 0), "\n")
cat("Samples with CTSK = 0:", sum(as.numeric(dog2["CTSK_3",]) == 0), "\n")
# All samples have detectable expression

# ── 2.3 Extract log2 expression ───────────────────────────────────────
df_cor <- data.frame(
  SRGN = log2(as.numeric(dog2["SRGN_3", ]) + 1),
  CTSK = log2(as.numeric(dog2["CTSK_3", ]) + 1),
  ACP5 = log2(as.numeric(dog2["ACP5_3", ]) + 1),
  MMP9 = log2(as.numeric(dog2["MMP9_3", ]) + 1),
  CD68 = log2(as.numeric(dog2["CD68_3", ]) + 1)
)

# ── 2.4 Spearman correlation tests ───────────────────────────────────
cat("SRGN vs CTSK:\n")
cor.test(df_cor$SRGN, df_cor$CTSK, method = "spearman")
# Result: R = 0.22, p = 0.0019

cat("\nSRGN vs ACP5:\n")
cor.test(df_cor$SRGN, df_cor$ACP5, method = "spearman")
# Result: R = 0.046, p = 0.52

cat("\nSRGN vs MMP9:\n")
cor.test(df_cor$SRGN, df_cor$MMP9, method = "spearman")
# Result: R = 0.052, p = 0.47

cat("\nSRGN vs CD68:\n")
cor.test(df_cor$SRGN, df_cor$CD68, method = "spearman")
# Result: R = 0.23, p = 0.0013

# ── 2.5 Four-panel correlation figure ────────────────────────────────
p1 <- ggplot(df_cor, aes(x = SRGN, y = CTSK)) +
  geom_point(alpha = 0.5, colour = "steelblue") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(title = "SRGN vs CTSK", x = "SRGN (log2)", y = "CTSK (log2)")

p2 <- ggplot(df_cor, aes(x = SRGN, y = ACP5)) +
  geom_point(alpha = 0.5, colour = "darkorange") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(title = "SRGN vs ACP5", x = "SRGN (log2)", y = "ACP5 (log2)")

p3 <- ggplot(df_cor, aes(x = SRGN, y = MMP9)) +
  geom_point(alpha = 0.5, colour = "darkgreen") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(title = "SRGN vs MMP9", x = "SRGN (log2)", y = "MMP9 (log2)")

p4 <- ggplot(df_cor, aes(x = SRGN, y = CD68)) +
  geom_point(alpha = 0.5, colour = "purple") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top") +
  theme_bw() +
  labs(title = "SRGN vs CD68", x = "SRGN (log2)", y = "CD68 (log2)")

combined <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
annotate_figure(combined,
  top = text_grob("SRGN correlates with osteoclast/myeloid markers in DOG2 (n=198)",
                  size = 13, face = "bold"))

ggsave("SRGN_markers_4panel.png", width = 12, height = 10, dpi = 300)


# ============================================================
# LAYER 3: CANINE SURVIVAL ANALYSIS (DOG2 CLINICAL)
# ============================================================

# ── 3.1 Load clinical data ────────────────────────────────────────────
clin <- read_excel("COTC022 SOC Data Set.xlsx")

# Calculate OS in days
clin$OS_days   <- as.numeric(difftime(clin$`Date off study`,
                                       clin$`Date of Surgery`,
                                       units = "days"))
clin$OS_status <- as.numeric(
  clin$`Cause of death (OSA, Unknown, Other)` == "OSA")

cat("Median OS:", median(clin$OS_days, na.rm = TRUE), "days\n")
cat("Events:", sum(clin$OS_status, na.rm = TRUE), "\n")

# ── 3.2 Match RNA-seq and clinical IDs ───────────────────────────────
rna_ids    <- sapply(colnames(dog2), function(x) strsplit(x, "_")[[1]][2])
common_ids <- intersect(rna_ids, clin$`Patient ID`)
cat("Matched samples:", length(common_ids), "\n")
# Result: 97 matched samples

# Extract SRGN expression for matched samples
srgn_matched <- log2(as.numeric(
  dog2["SRGN_3", names(rna_ids[rna_ids %in% common_ids])]) + 1)
names(srgn_matched) <- rna_ids[rna_ids %in% common_ids]

clin_matched <- clin[match(names(srgn_matched), clin$`Patient ID`), ]

# ── 3.3 Build survival dataframe ──────────────────────────────────────
surv_df <- data.frame(
  OS_days    = clin_matched$OS_days,
  OS_status  = clin_matched$OS_status,
  SRGN_expr  = srgn_matched,
  SRGN_group = ifelse(srgn_matched >= median(srgn_matched),
                      "SRGN High", "SRGN Low")
)

# ── 3.4 Kaplan-Meier ──────────────────────────────────────────────────
fit_canine <- survfit(Surv(OS_days, OS_status) ~ SRGN_group, data = surv_df)

ggsurvplot(fit_canine,
           data          = surv_df,
           pval          = TRUE,
           risk.table    = TRUE,
           tables.height = 0.25,
           palette       = c("#E41A1C", "#377EB8"),
           legend.title  = "",
           xlab          = "Time (days)",
           ylab          = "Overall Survival",
           title         = "SRGN and OS in canine OSA (n=97)",
           ggtheme       = theme_bw())

ggsave("SRGN_survival_KM.png", width = 8, height = 7, dpi = 300)
# Result: p = 0.94 — no significant association


# ============================================================
# LAYER 3B: HUMAN SURVIVAL ANALYSIS (GSE39055)
# ============================================================
# Moved to a standalone, self-contained script to avoid duplication:
#   srgn_survival_GSE39055.R
# That script pulls GSE39055 directly via getGEO() (no local file needed) and
# reproduces the human recurrence-free survival analysis (Fig. 3): median-split
# KM, optimal-cutpoint KM, and Cox PH. Run it as a separate step.


# ============================================================
# LAYER 4: TF REGULON ANALYSIS (DoRothEA / VIPER)
# ============================================================

# ── 4.1 Load regulons ─────────────────────────────────────────────────
data(dorothea_hs, package = "dorothea")
regulons     <- dorothea_hs %>% filter(confidence %in% c("A", "B"))
regulon_list <- df2regulon(regulons)

cat("TFs in regulon database:", length(regulon_list), "\n")

# ── 4.2 Extract osteoclast expression matrix ──────────────────────────
merged_filtered <- readRDS("GSE252470_annotated.rds")
merged_filtered <- JoinLayers(merged_filtered)

osteo_cells  <- colnames(merged_filtered)[merged_filtered$cell_type == "Osteoclast"]
expr_matrix  <- GetAssayData(merged_filtered, layer = "data")
expr_dense   <- as.matrix(expr_matrix[, osteo_cells])
cat("Osteoclast matrix:", nrow(expr_dense), "genes x",
    ncol(expr_dense), "cells\n")

# Free memory before VIPER
rm(expr_matrix); gc()

# ── 4.3 Run VIPER ─────────────────────────────────────────────────────
tf_activity <- viper(eset    = expr_dense,
                     regulon = regulon_list,
                     minsize = 4,
                     method  = "scale",
                     verbose = FALSE)

cat("VIPER complete\n")
cat("TF activity matrix:", nrow(tf_activity), "TFs x",
    ncol(tf_activity), "cells\n")
# Result: 117 TFs x 4,047 cells

# ── 4.4 Identify TFs enriched in SRGN-high osteoclasts ───────────────
srgn_osteo_vip <- as.numeric(expr_dense["SRGN", ])
srgn_group_vip <- ifelse(srgn_osteo_vip >= median(srgn_osteo_vip),
                          "SRGN_high", "SRGN_low")

tf_pvals  <- apply(tf_activity, 1, function(tf) {
  wilcox.test(tf[srgn_group_vip == "SRGN_high"],
              tf[srgn_group_vip == "SRGN_low"])$p.value
})

tf_effect <- apply(tf_activity, 1, function(tf) {
  mean(tf[srgn_group_vip == "SRGN_high"]) -
  mean(tf[srgn_group_vip == "SRGN_low"])
})

tf_results <- data.frame(
  TF   = rownames(tf_activity),
  effect = tf_effect,
  pval   = tf_pvals,
  padj   = p.adjust(tf_pvals, method = "fdr")
) %>% arrange(padj)

# Top TFs enriched in SRGN-high osteoclasts
head(tf_results %>% filter(effect > 0), 10)
# Result: FOSL2, CEBPA, STAT5B, ETS1, SRF (all FDR < 10^-130)

# ── 4.5 Volcano plot ──────────────────────────────────────────────────
tf_results$log10padj   <- -log10(tf_results$padj)
tf_results$significant <- tf_results$padj < 0.05 & abs(tf_results$effect) > 0.5

top_tfs <- c("FOSL2", "CEBPA", "STAT5B", "ETS1", "SRF", "FOXO4")
tf_results$label <- ifelse(tf_results$TF %in% top_tfs, tf_results$TF, "")

ggplot(tf_results, aes(x = effect, y = log10padj,
                        colour = significant, label = label)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_text_repel(size = 4, fontface = "bold", max.overlaps = 20) +
  scale_colour_manual(values = c("grey70", "#E41A1C")) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", colour = "grey50") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title      = element_text(face = "bold")) +
  labs(title = "TF activity in SRGN-high vs SRGN-low osteoclasts",
       x     = "Mean difference in TF activity",
       y     = "-log10(FDR)")

ggsave("TF_activity_volcano.png", width = 8, height = 6, dpi = 300)

saveRDS(tf_results, "TF_activity_results.rds")


# ============================================================
# LAYER 5: HUMAN scRNA-seq ANALYSIS (GSE162454)
# ============================================================

# ── 5.1 Extract TAR archive ───────────────────────────────────────────
untar("GSE162454_RAW.tar", exdir = "GSE162454/raw")

# ── 5.2 Organise into per-sample subfolders ───────────────────────────
sample_ids_human <- c("OS_1", "OS_2", "OS_3", "OS_4", "OS_5", "OS_6")
gsm_ids <- c("GSM4952363", "GSM4952364", "GSM4952365",
             "GSM5155198", "GSM5155199", "GSM5155200")

for (i in seq_along(sample_ids_human)) {
  dir.create(file.path("GSE162454", sample_ids_human[i]),
             showWarnings = FALSE)
  file.copy(
    file.path("GSE162454/raw",
              paste0(gsm_ids[i], "_", sample_ids_human[i],
                     "_barcodes.tsv.gz")),
    file.path("GSE162454", sample_ids_human[i], "barcodes.tsv.gz"))
  file.copy(
    file.path("GSE162454/raw",
              paste0(gsm_ids[i], "_", sample_ids_human[i],
                     "_features.tsv.gz")),
    file.path("GSE162454", sample_ids_human[i], "features.tsv.gz"))
  file.copy(
    file.path("GSE162454/raw",
              paste0(gsm_ids[i], "_", sample_ids_human[i],
                     "_matrix.mtx.gz")),
    file.path("GSE162454", sample_ids_human[i], "matrix.mtx.gz"))
}

# ── 5.3 Load all 6 human samples ─────────────────────────────────────
human_list <- lapply(sample_ids_human, function(sid) {
  counts <- Read10X(data.dir = file.path("GSE162454", sid))
  CreateSeuratObject(counts,
                     project      = sid,
                     min.cells    = 3,
                     min.features = 200)
})
names(human_list) <- sample_ids_human

# Merge
human_merged <- merge(human_list[[1]],
                      y            = human_list[-1],
                      add.cell.ids = sample_ids_human)
dim(human_merged)
# Result: 24,611 genes x 50,780 cells

# Confirm SRGN present
rownames(human_merged)[grep("SRGN", rownames(human_merged),
                             ignore.case = TRUE)]
# Result: "SRGN"

# ── 5.4 QC filtering ──────────────────────────────────────────────────
human_merged[["percent.mt"]] <- PercentageFeatureSet(human_merged,
                                                      pattern = "^MT-")

VlnPlot(human_merged,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

# Filter: 200-7,000 genes, <60,000 UMI, <20% mitochondrial
human_filtered <- subset(human_merged,
                          subset = nFeature_RNA > 200  &
                                   nFeature_RNA < 7000 &
                                   nCount_RNA   < 60000 &
                                   percent.mt   < 20)
dim(human_filtered)
# Result: 24,611 genes x 44,086 cells

# ── 5.5 Seurat preprocessing ─────────────────────────────────────────
human_filtered <- human_filtered %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(dims = 1:20)

# ── 5.6 SRGN on human UMAP ───────────────────────────────────────────
FeaturePlot(human_filtered,
            features = "SRGN",
            cols     = c("lightgrey", "red"),
            pt.size  = 0.3) +
  ggtitle("SRGN expression in human OSA (GSE162454)")

ggsave("SRGN_human_scRNAseq_UMAP.png", width = 8, height = 6, dpi = 300)

# ── 5.7 Cell type marker DotPlot ─────────────────────────────────────
DotPlot(human_filtered,
        features  = c("PTPRC", "CD3E", "CD79A",
                       "CD68", "LYZ", "CTSK", "ACP5",
                       "FAP", "DCN", "PECAM1",
                       "RUNX2", "SPP1", "SRGN"),
        group.by  = "seurat_clusters",
        cols      = c("lightgrey", "darkred"),
        dot.scale = 6) +
  RotatedAxis() +
  scale_y_discrete(expand = expansion(add = 0.8)) +
  ggtitle("Cell type markers in human OSA")

ggsave("human_dotplot_markers.png", width = 12, height = 14, dpi = 300)

# ── 5.8 Cell type annotation ──────────────────────────────────────────
human_labels <- c(
  "0"  = "Tumor",      "1"  = "Tumor",      "2"  = "Tumor",
  "3"  = "Tumor",      "4"  = "Tumor",      "5"  = "Osteoclast",
  "6"  = "Tumor",      "7"  = "Tumor",      "8"  = "Tumor",
  "9"  = "Tumor",      "10" = "Tumor",      "11" = "CAF",
  "12" = "CAF",        "13" = "CAF",        "14" = "Endothelial",
  "15" = "Macrophage", "16" = "Tumor",      "17" = "T cell",
  "18" = "Tumor",      "19" = "B cell",     "20" = "T cell",
  "21" = "Macrophage", "22" = "Tumor",      "23" = "Osteoclast",
  "24" = "Tumor"
)

human_filtered[["cell_type"]] <- dplyr::recode(
  as.character(human_filtered$seurat_clusters),
  !!!human_labels)

# ── 5.9 SRGN by cell type ─────────────────────────────────────────────
VlnPlot(human_filtered,
        features = "SRGN",
        group.by = "cell_type",
        pt.size  = 0) +
  ggtitle("SRGN expression by cell type in human OSA (GSE162454)") +
  xlab("") + ylab("Normalized expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("SRGN_human_celltype.png", width = 10, height = 6, dpi = 300)

# ── 5.10 Quantify SRGN-positive cells ────────────────────────────────
human_filtered <- JoinLayers(human_filtered)

osteo_h      <- subset(human_filtered, cell_type == "Osteoclast")
tumor_h      <- subset(human_filtered, cell_type == "Tumor")
srgn_osteo_h <- FetchData(osteo_h, vars = "SRGN", layer = "counts")
srgn_tumor_h <- FetchData(tumor_h, vars = "SRGN", layer = "counts")

cat("Human osteoclast SRGN+:",
    round(mean(srgn_osteo_h$SRGN > 0) * 100, 2), "%\n")
cat("Human tumor SRGN+:",
    round(mean(srgn_tumor_h$SRGN > 0) * 100, 2), "%\n")
# Result: Osteoclast 62.34%, Tumor 70.25%

# Chi-square test
contingency_h <- matrix(
  c(sum(srgn_osteo_h$SRGN > 0), sum(srgn_osteo_h$SRGN == 0),
    sum(srgn_tumor_h$SRGN > 0), sum(srgn_tumor_h$SRGN == 0)),
  nrow     = 2,
  dimnames = list(c("SRGN+", "SRGN-"), c("Osteoclast", "Tumor"))
)
print(contingency_h)
chisq.test(contingency_h)
# Result: chi-sq = 73.014, p < 2.2e-16

# Contamination check: is SRGN in tumor cells correlated with osteoclast markers?
tumor_data <- FetchData(tumor_h,
                         vars  = c("SRGN", "CTSK", "ACP5"),
                         layer = "counts")
cor.test(tumor_data$SRGN, tumor_data$CTSK, method = "spearman")
cor.test(tumor_data$SRGN, tumor_data$ACP5, method = "spearman")

# Save
saveRDS(human_filtered, "GSE162454_annotated.rds")


# ============================================================
# CROSS-SPECIES COMPARISON FIGURE
# ============================================================

# ── Summary dataframe ─────────────────────────────────────────────────
comparison_df <- data.frame(
  Species   = c("Canine OSA\n(GSE252470)", "Canine OSA\n(GSE252470)",
                "Human OSA\n(GSE162454)",  "Human OSA\n(GSE162454)"),
  Cell_type = c("Osteoclast", "Tumor", "Osteoclast", "Tumor"),
  SRGN_pct  = c(72.89, 17.88, 62.34, 70.25)
)

ggplot(comparison_df, aes(x = Cell_type, y = SRGN_pct, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = paste0(SRGN_pct, "%")),
            position = position_dodge(width = 0.6),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  scale_y_continuous(limits = c(0, 85), expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "top",
        axis.text       = element_text(size = 12),
        axis.title      = element_text(size = 13),
        legend.text     = element_text(size = 11),
        plot.title      = element_text(face = "bold", size = 14)) +
  labs(title = "SRGN expression across cell types: canine vs human OSA",
       x     = "Cell type",
       y     = "% SRGN-positive cells",
       fill  = "")

ggsave("SRGN_cross_species_comparison.png", width = 8, height = 6, dpi = 300)

# ── Cross-species summary table ───────────────────────────────────────
summary_table <- data.frame(
  Species      = c("Canine", "Canine", "Human", "Human"),
  Cell_type    = c("Osteoclast", "Tumor", "Osteoclast", "Tumor"),
  SRGN_pos_pct = c(72.89, 17.88, 62.34, 70.25),
  Dataset      = c("GSE252470", "GSE252470", "GSE162454", "GSE162454"),
  n_cells      = c(4047, 23067, 2663, 34556)
)
print(summary_table)
write.csv(summary_table, "cross_species_SRGN_summary.csv", row.names = FALSE)


# ============================================================
# SESSION INFO
# ============================================================
sessionInfo()
