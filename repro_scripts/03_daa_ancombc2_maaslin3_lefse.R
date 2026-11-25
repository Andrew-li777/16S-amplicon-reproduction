# ==============================================================================
# 03_daa_ancombc2_maaslin3_lefse.R — Track A-2 组成性差异丰度多方法交叉验证
#
# 核心原则:
#   - 严禁输入抽平后的计数表，必须使用未抽平的原始计数表 (01_ps_miseq.rds)。
#   - 在属 (Genus) 水平并行对比三大主流/前沿方法与传统基准:
#       1. ANCOMBC::ancombc2 (抽样分数偏差校正 + 结构零 + 伪计数敏感性)
#       2. maaslin3::maaslin3 (Nature Methods 2025/2026: 丰度+流行率联合模型)
#       3. microeco::trans_diff(method = "lefse") (经典 LDA 效应量基准)
#       4. 传统 TSS 相对丰度 Wilcoxon (BH FDR, 用作基线参照)
#   - 输出一致性对比表与效应量双联图。
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(maaslin3)
  library(microeco)
  library(file2meco)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

get_project_root <- function() {
  c_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", c_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
    if (basename(script_dir) == "repro_scripts") {
      return(normalizePath(file.path(script_dir, ".."), mustWork = TRUE))
    }
    return(script_dir)
  }
  normalizePath(getwd(), mustWork = TRUE)
}

proj_root <- get_project_root()
setwd(proj_root)

SEED <- 20261125L
dir.create("data/processed/03_maaslin3_out", recursive = TRUE, showWarnings = FALSE)
dir.create("repro_results", recursive = TRUE, showWarnings = FALSE)

# 1. 加载未抽平 phyloseq 对象并聚合并规范属水平名称 ----------------------------
ps <- readRDS("data/processed/01_ps_miseq.rds")
ps_gen <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

raw_genera <- as.character(tax_table(ps_gen)[, "Genus"])
clean_genera <- make.unique(gsub("[^A-Za-z0-9_.-]", "_", raw_genera), sep = "_")
taxa_names(ps_gen) <- clean_genera
tax_table(ps_gen)[, "Genus"] <- clean_genera

# 流行率过滤 (>= 10% 样本检出)
otu_gen <- as(otu_table(ps_gen), "matrix")
if (taxa_are_rows(ps_gen)) otu_gen <- t(otu_gen)
prev_vec <- colSums(otu_gen > 0) / nrow(otu_gen)
keep_gen <- names(prev_vec)[prev_vec >= 0.10]
ps_gen <- prune_taxa(keep_gen, ps_gen)
otu_gen <- otu_gen[, keep_gen, drop = FALSE]
meta_df <- data.frame(sample_data(ps_gen))
meta_df$When <- factor(meta_df$When, levels = c("Early", "Late"))

cat(sprintf("[DAA] 未抽平属水平筛选后保留 %d 个属 (流行率 >= 10%%, %d 样本)\n",
            length(keep_gen), nrow(otu_gen)))

# 2. Method 1: ANCOM-BC2 (Windows 必须保持 n_cl = 1) ---------------------------
cat("[METHOD 1] 运行 ANCOM-BC2 (pseudo_sens = TRUE, struc_zero = TRUE, n_cl = 1)...\n")
set.seed(SEED)
anc_out <- ANCOMBC::ancombc2(
  data = ps_gen,
  taxa_are_rows = FALSE,
  fix_formula = "When",
  p_adj_method = "BH",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 0,
  group = "When",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE
)
anc_dt <- as.data.table(anc_out$res)
anc_clean <- anc_dt[, .(
  Genus = taxon,
  ANCOMBC2_LFC = round(lfc_WhenLate, 4),
  ANCOMBC2_SE = round(se_WhenLate, 4),
  ANCOMBC2_Q = q_WhenLate,
  ANCOMBC2_Passed_SS = passed_ss_WhenLate,
  ANCOMBC2_Sig = (diff_WhenLate == TRUE)
)]

# 3. Method 2: MaAsLin 3 (丰度 + 流行率联合广义线性模型) -----------------------
cat("[METHOD 2] 运行 MaAsLin 3 (TSS + LOG + median_comparison_abundance = TRUE)...\n")
set.seed(SEED)
m3_out <- maaslin3::maaslin3(
  input_data = as.data.frame(otu_gen),
  input_metadata = meta_df,
  output = "data/processed/03_maaslin3_out",
  fixed_effects = "When",
  reference = c("When,Early"),
  min_prevalence = 0.10,
  normalization = "TSS",
  transform = "LOG",
  median_comparison_abundance = TRUE,
  median_comparison_prevalence = FALSE,
  max_significance = 0.05,
  plot_summary_plot = FALSE,
  plot_associations = FALSE,
  cores = 1,
  verbosity = "ERROR",
  warn_prevalence = FALSE,
  bypass_small_group_warning = TRUE
)

m3_ab <- as.data.table(m3_out$fit_data_abundance$results)[metadata == "When" & value == "Late"]
m3_pr <- as.data.table(m3_out$fit_data_prevalence$results)[metadata == "When" & value == "Late"]

m3_clean <- merge(
  m3_ab[, .(Genus = feature, MaAsLin3_Coef = round(coef, 4),
            MaAsLin3_Q_Abund = qval_individual, MaAsLin3_Q_Joint = qval_joint)],
  m3_pr[, .(Genus = feature, MaAsLin3_Q_Prev = qval_individual)],
  by = "Genus", all = TRUE
)
m3_clean[, MaAsLin3_Sig := (!is.na(MaAsLin3_Q_Joint) & MaAsLin3_Q_Joint < 0.05) |
                           (!is.na(MaAsLin3_Q_Abund) & MaAsLin3_Q_Abund < 0.05)]

# 4. Method 3: microeco LEfSe --------------------------------------------------
cat("[METHOD 3] 运行 microeco::trans_diff (LEfSe)...\n")
mt_gen <- file2meco::phyloseq2meco(ps_gen)
mt_gen$sample_table$When <- factor(mt_gen$sample_table$When, levels = c("Early", "Late"))
mt_gen$tidy_dataset()

set.seed(SEED)
lefse_obj <- microeco::trans_diff$new(
  dataset = mt_gen,
  method = "lefse",
  group = "When",
  taxa_level = "Genus",
  alpha = 0.05,
  lefse_subgroup = NULL
)
lefse_raw <- as.data.table(lefse_obj$res_diff)
lefse_raw[, Genus := sub(".*\\|g__", "", Taxa)]
lefse_raw[, Genus := sub("^g__", "", Genus)]
lefse_clean <- lefse_raw[, .(
  Genus,
  LEfSe_Enriched_Group = Group,
  LEfSe_LDA = round(LDA, 4),
  LEfSe_Q = P.adj,
  LEfSe_Sig = (P.adj < 0.05 & LDA >= 2.0)
)]
# 若同一属名对应多行则取最高 LDA 行
lefse_clean <- lefse_clean[order(-LEfSe_LDA), .SD[1], by = Genus]

# 5. Method 4: 传统 TSS Wilcoxon BH FDR 基准对照 -------------------------------
tss_gen <- sweep(otu_gen, 1L, rowSums(otu_gen), "/")
wilcox_dt <- rbindlist(lapply(keep_gen, function(g) {
  v_e <- tss_gen[meta_df$When == "Early", g]
  v_l <- tss_gen[meta_df$When == "Late", g]
  pv <- suppressWarnings(wilcox.test(v_l, v_e)$p.value)
  l2fc <- log2((mean(v_l) + 1e-5) / (mean(v_e) + 1e-5))
  data.table(Genus = g, Mean_Rel_Early = round(mean(v_e), 5), Mean_Rel_Late = round(mean(v_l), 5),
             TSS_Log2FC = round(l2fc, 4), TSS_P = pv)
}))
wilcox_dt[, TSS_FDR := p.adjust(TSS_P, method = "BH")]
wilcox_dt[, TSS_Wilcoxon_Sig := (TSS_FDR < 0.05)]

# 6. 合并交叉验证全表并导出 ----------------------------------------------------
cons_dt <- Reduce(function(x, y) merge(x, y, by = "Genus", all = TRUE),
                  list(wilcox_dt, anc_clean, m3_clean, lefse_clean))
for (col in c("ANCOMBC2_Sig", "MaAsLin3_Sig", "LEfSe_Sig", "TSS_Wilcoxon_Sig")) {
  cons_dt[is.na(get(col)), (col) := FALSE]
}
cons_dt[, N_Modern_Calls := as.integer(ANCOMBC2_Sig) + as.integer(MaAsLin3_Sig) + as.integer(LEfSe_Sig)]
cons_dt[, Consensus_2plus := (N_Modern_Calls >= 2L)]
cons_dt <- cons_dt[order(-N_Modern_Calls, -abs(ANCOMBC2_LFC))]

fwrite(cons_dt, "repro_results/03_daa_consensus_genera.tsv", sep = "\t")

n_anc  <- cons_dt[ANCOMBC2_Sig == TRUE, .N]
n_m3   <- cons_dt[MaAsLin3_Sig == TRUE, .N]
n_lef  <- cons_dt[LEfSe_Sig == TRUE, .N]
n_tss  <- cons_dt[TSS_Wilcoxon_Sig == TRUE, .N]
n_cons <- cons_dt[Consensus_2plus == TRUE, .N]

daa_summary <- data.table(
  Metric = c(
    "Total_Tested_Genera",
    "ANCOMBC2_Sig_Genera",
    "MaAsLin3_Sig_Genera",
    "LEfSe_Sig_Genera",
    "Naive_TSS_Wilcoxon_Sig_Genera",
    "Consensus_2plus_Genera"
  ),
  Value = c(
    as.character(nrow(cons_dt)),
    as.character(n_anc),
    as.character(n_m3),
    as.character(n_lef),
    as.character(n_tss),
    as.character(n_cons)
  )
)
fwrite(daa_summary, "repro_results/03_daa_method_comparison.tsv", sep = "\t")
cat(sprintf("[DAA DONE] Tested=%d | ANCOM-BC2=%d | MaAsLin3=%d | LEfSe=%d | TSS=%d | Consensus(>=2)=%d\n",
            nrow(cons_dt), n_anc, n_m3, n_lef, n_tss, n_cons))

# 7. 绘制交叉验证双联图 --------------------------------------------------------
cons_dt[, Sort_Effect := fifelse(!is.na(ANCOMBC2_LFC), abs(ANCOMBC2_LFC), abs(MaAsLin3_Coef))]
cons_dt[is.na(Sort_Effect), Sort_Effect := 0]
top_plot_dt <- cons_dt[order(-N_Modern_Calls, -Sort_Effect)][1:min(18L, .N)]
top_plot_dt[, Genus_Label := factor(Genus, levels = rev(Genus))]

mat_long <- melt(
  top_plot_dt,
  id.vars = "Genus_Label",
  measure.vars = c("ANCOMBC2_Sig", "MaAsLin3_Sig", "LEfSe_Sig", "TSS_Wilcoxon_Sig"),
  variable.name = "Method", value.name = "Detected"
)
mat_long[, Method := factor(
  Method,
  levels = c("ANCOMBC2_Sig", "MaAsLin3_Sig", "LEfSe_Sig", "TSS_Wilcoxon_Sig"),
  labels = c("ANCOM-BC2", "MaAsLin 3", "LEfSe", "TSS Wilcoxon")
)]

p_heat <- ggplot(mat_long, aes(x = Method, y = Genus_Label, fill = Detected)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = ifelse(Detected, "Sig", "")), size = 3.0, color = "white", fontface = "bold", family = "sans") +
  scale_fill_manual(values = c("FALSE" = "#ECF0F1", "TRUE" = "#1F618D"), guide = "none") +
  labs(
    title = "A. Multi-Method Differential Genus Concordance",
    subtitle = sprintf("Consensus (>=2 modern tools): %d genera", n_cons),
    x = NULL, y = "Genus (Prevalence >= 10%)"
  ) +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), plot.title = element_text(face = "bold"))

lfc_df <- top_plot_dt[!is.na(ANCOMBC2_LFC)]
lfc_df[, Direction := ifelse(ANCOMBC2_LFC >= 0, "Enriched in Late", "Enriched in Early")]
na_df <- top_plot_dt[is.na(ANCOMBC2_LFC)]

p_lfc <- ggplot(lfc_df, aes(x = ANCOMBC2_LFC, y = Genus_Label, color = Direction, shape = Consensus_2plus)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7F8C8D") +
  geom_errorbar(aes(xmin = ANCOMBC2_LFC - ANCOMBC2_SE, xmax = ANCOMBC2_LFC + ANCOMBC2_SE),
                orientation = "y", width = 0.28, linewidth = 0.7) +
  geom_point(size = 3.4) +
  geom_text(data = na_df, aes(x = 0, y = Genus_Label, label = "Filtered in ANCOM-BC2 (structural zero / sparse)"),
            inherit.aes = FALSE, size = 2.8, color = "#7F8C8D", fontface = "italic", family = "sans") +
  scale_y_discrete(limits = levels(top_plot_dt$Genus_Label)) +
  scale_color_manual(values = c("Enriched in Late" = "#D55E00", "Enriched in Early" = "#0072B2")) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), labels = c("Single/None", "Consensus >=2")) +
  labs(
    title = "B. ANCOM-BC2 Bias-Corrected Log Fold Change (+/- SE)",
    subtitle = "Unrarefied genus counts with sampling-fraction offset",
    x = "ANCOM-BC2 Log Fold Change (Late vs Early)",
    y = NULL, color = "Direction", shape = "Validation"
  ) +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  "repro_results/03_daa_consensus_and_effect_sizes.png",
  p_heat + p_lfc + plot_layout(widths = c(1, 1.35)),
  width = 13.5, height = 6.8, dpi = 300, bg = "white"
)
