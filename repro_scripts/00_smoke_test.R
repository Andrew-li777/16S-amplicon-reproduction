# ==============================================================================
# 00_smoke_test.R — 16S 扩增子 10 秒确定性冒烟自检脚本
#
# 核心目标 (REPRO-D):
#   1. 验证路径解析、文件读写与图形设备全线畅通，不依赖外网下载与大型外部参考库。
#   2. 构造含"成分闭合效应 (Compositional Closure Artifact)"与"试剂盒污染 (Kitome)"的
#      确定性合成 16S 计数矩阵 (24 样本 × 40 ASV)：
#      - 断言 1: 传统 TSS 相对丰度检验对 35 个真实不变背景菌产生大量假阳性；
#      - 断言 2: CLR 中心对数比变换大幅压低背景假阳性并精准检出 3 个真实变化菌；
#      - 断言 3: 阴性对照 (Neg_Ctrl) 精准捕获 2 个试剂盒背景污染 ASV。
#   3. 输出诊断双联图 (00_smoke_16s_compositional_qc.png) 与定量摘要 TSV。
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("repro_results", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

# 1. 构建确定性合成 16S 样本元数据与计数矩阵 --------------------------------
SEED <- 20261125L
set.seed(SEED)

n_ctrl <- 10L
n_dis  <- 10L
n_neg  <- 4L
n_samp <- n_ctrl + n_dis + n_neg
n_asv  <- 40L

sample_ids <- c(
  sprintf("Ctrl_%02d", seq_len(n_ctrl)),
  sprintf("Dis_%02d", seq_len(n_dis)),
  sprintf("Neg_%02d", seq_len(n_neg))
)
groups <- c(rep("Control", n_ctrl), rep("Disease", n_dis), rep("Neg_Ctrl", n_neg))
dna_conc <- c(
  rnorm(n_ctrl, mean = 18.0, sd = 2.5),
  rnorm(n_dis,  mean = 18.5, sd = 2.5),
  rnorm(n_neg,  mean = 0.35, sd = 0.08)
)
dna_conc <- pmax(dna_conc, 0.1)

asv_ids <- c(
  "ASV01_Bloom",
  "ASV02_TrueUp",
  "ASV03_TrueDown",
  "ASV04_Kitome1",
  "ASV05_Kitome2",
  sprintf("ASV%02d_Invariant", 6:n_asv)
)
truth_class <- c(
  "True_Bloom",
  "True_Differential",
  "True_Differential",
  "Kitome_Contaminant",
  "Kitome_Contaminant",
  rep("Invariant_Background", n_asv - 5L)
)

# 绝对丰度期望矩阵 (ASV × Sample)
abs_mat <- matrix(0, nrow = n_asv, ncol = n_samp, dimnames = list(asv_ids, sample_ids))
# 35 个真实不变背景菌：Control 与 Disease 保持完全相同的绝对菌量 (~800 counts)
for (i in 6:n_asv) {
  base_mu <- runif(1, 650, 1150)
  abs_mat[i, groups == "Control"]  <- rnbinom(n_ctrl, mu = base_mu, size = 45)
  abs_mat[i, groups == "Disease"]  <- rnbinom(n_dis,  mu = base_mu, size = 45)
  abs_mat[i, groups == "Neg_Ctrl"] <- rbinom(n_neg, size = 2, prob = 0.1)
}
# 1 个高丰度爆发菌 (ASV01_Bloom)：引发 Disease 组强烈成分闭合挤压效应
abs_mat[1, groups == "Control"]  <- rnbinom(n_ctrl, mu = 1500,  size = 40)
abs_mat[1, groups == "Disease"]  <- rnbinom(n_dis,  mu = 16500, size = 40)
abs_mat[1, groups == "Neg_Ctrl"] <- rbinom(n_neg, size = 3, prob = 0.1)

# 2 个中度真实差异菌
abs_mat[2, groups == "Control"]  <- rnbinom(n_ctrl, mu = 320,  size = 40)
abs_mat[2, groups == "Disease"]  <- rnbinom(n_dis,  mu = 1900, size = 40)
abs_mat[3, groups == "Control"]  <- rnbinom(n_ctrl, mu = 1400, size = 40)
abs_mat[3, groups == "Disease"]  <- rnbinom(n_dis,  mu = 220,  size = 40)

# 2 个试剂盒背景污染菌：固定微量 DNA 混入，在低生物量 Neg_Ctrl 中占主导
abs_mat[4, ] <- rnbinom(n_samp, mu = 260, size = 35)
abs_mat[5, ] <- rnbinom(n_samp, mu = 210, size = 35)

# 模拟常规测序文库标准化抽样（生物样本 ~25,000 reads，阴性对照 ~1,500 reads）
set.seed(SEED)
counts_mat <- abs_mat
for (j in seq_len(n_samp)) {
  target_depth <- if (groups[j] == "Neg_Ctrl") 1500L else 25000L
  probs <- abs_mat[, j] / sum(abs_mat[, j])
  counts_mat[, j] <- as.vector(rmultinom(1L, size = target_depth, prob = probs))
}

# 2. 对比分析：朴素 TSS 相对丰度 vs CLR 中心对数比变换 -------------------------
bio_idx <- which(groups %in% c("Control", "Disease"))
bio_counts <- counts_mat[, bio_idx]
bio_groups <- groups[bio_idx]
neg_idx <- which(groups == "Neg_Ctrl")

# 相对丰度 (TSS)
tss_all <- sweep(counts_mat, 2L, colSums(counts_mat), "/")
tss_bio <- tss_all[, bio_idx]

# CLR 变换（以稳定伪计数 +0.5 取对数并按样本中心化）
log_bio <- log2(bio_counts + 0.5)
clr_bio <- sweep(log_bio, 2L, colMeans(log_bio), "-")

eval_list <- lapply(seq_len(n_asv), function(i) {
  c_tss <- tss_bio[i, bio_groups == "Control"]
  d_tss <- tss_bio[i, bio_groups == "Disease"]
  p_tss <- suppressWarnings(wilcox.test(d_tss, c_tss)$p.value)
  log2fc_tss <- log2((mean(d_tss) + 1e-6) / (mean(c_tss) + 1e-6))

  c_clr <- clr_bio[i, bio_groups == "Control"]
  d_clr <- clr_bio[i, bio_groups == "Disease"]
  p_clr <- suppressWarnings(wilcox.test(d_clr, c_clr)$p.value)
  diff_clr <- mean(d_clr) - mean(c_clr)

  neg_mean_tss <- mean(tss_all[i, neg_idx])
  bio_mean_tss <- mean(tss_bio[i, ])
  neg_bio_ratio <- (neg_mean_tss + 1e-6) / (bio_mean_tss + 1e-6)

  data.table(
    ASV = asv_ids[i],
    Truth_Class = truth_class[i],
    Mean_TSS_Ctrl = round(mean(c_tss), 5),
    Mean_TSS_Dis = round(mean(d_tss), 5),
    Log2FC_TSS = round(log2fc_tss, 4),
    P_TSS = p_tss,
    Mean_CLR_Diff = round(diff_clr, 4),
    P_CLR = p_clr,
    Neg_Mean_TSS = round(neg_mean_tss, 5),
    Neg_to_Bio_Ratio = round(neg_bio_ratio, 2)
  )
})

audit_dt <- rbindlist(eval_list)
audit_dt[, FDR_TSS := p.adjust(P_TSS, method = "BH")]
audit_dt[, FDR_CLR := p.adjust(P_CLR, method = "BH")]
audit_dt[, Kitome_Flag := Neg_to_Bio_Ratio >= 5.0 & Neg_Mean_TSS >= 0.05]

# 3. 定量断言防线 ---------------------------------------------------------------
naive_fp <- audit_dt[Truth_Class == "Invariant_Background" & FDR_TSS < 0.05, .N]
clr_fp   <- audit_dt[Truth_Class == "Invariant_Background" & FDR_CLR < 0.05, .N]
true_rec <- audit_dt[Truth_Class %in% c("True_Bloom", "True_Differential") & FDR_CLR < 0.05, .N]
kitome_n <- audit_dt[Truth_Class == "Kitome_Contaminant" & Kitome_Flag == TRUE, .N]
kitome_fp <- audit_dt[Truth_Class != "Kitome_Contaminant" & Kitome_Flag == TRUE, .N]

stopifnot(
  naive_fp >= 15L,
  clr_fp <= 3L,
  clr_fp < naive_fp,
  true_rec == 3L,
  kitome_n == 2L,
  kitome_fp == 0L
)

# 4. 导出定量审计表与汇总指标 ---------------------------------------------------
fwrite(audit_dt, "data/processed/smoke_16s_feature_audit.tsv", sep = "\t")

summary_dt <- data.table(
  Metric = c(
    "Smoke_Samples_Total",
    "Smoke_ASVs_Total",
    "Naive_TSS_False_Positives",
    "CLR_False_Positives",
    "True_Diff_Recovered_CLR",
    "Kitome_Contaminants_Flagged",
    "Kitome_False_Positives",
    "Random_Seed",
    "Smoke_Status"
  ),
  Value = c(
    as.character(n_samp),
    as.character(n_asv),
    as.character(naive_fp),
    as.character(clr_fp),
    as.character(true_rec),
    as.character(kitome_n),
    as.character(kitome_fp),
    as.character(SEED),
    "PASS"
  )
)
fwrite(summary_dt, "repro_results/00_smoke_verification_summary.tsv", sep = "\t")

# 5. 绘制双联诊断质控图 ---------------------------------------------------------
pal_truth <- c(
  "True_Bloom"           = "#D55E00",
  "True_Differential"    = "#0072B2",
  "Invariant_Background" = "#7F8C8D",
  "Kitome_Contaminant"   = "#CC79A7"
)

p_comp <- ggplot(
  audit_dt[Truth_Class != "Kitome_Contaminant"],
  aes(x = Log2FC_TSS, y = Mean_CLR_Diff, color = Truth_Class)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#95A5A6") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#95A5A6") +
  geom_point(aes(size = ifelse(Truth_Class == "Invariant_Background", 2.4, 4.0)), alpha = 0.88) +
  scale_size_identity() +
  scale_color_manual(values = pal_truth) +
  annotate(
    "label", x = -2.95, y = 2.75, hjust = 0, vjust = 1, size = 3.2, family = "sans",
    label = sprintf("Naive TSS FP (Invariant): %d / 35\nCLR FP (Invariant): %d / 35\nTrue Diff Recovered (CLR): %d / 3",
                    naive_fp, clr_fp, true_rec),
    fill = "#F8F9FA", color = "#2C3E50"
  ) +
  labs(
    title = "A. Compositional Closure Artifact vs CLR Correction",
    subtitle = "Single blooming taxon shifts TSS proportions of 35 truly invariant taxa",
    x = "Naive Relative Abundance log2FC (Disease / Control)",
    y = "CLR Mean Difference (Disease - Control)",
    color = "Ground Truth"
  ) +
  theme_bw(base_family = "sans", base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

kit_df <- data.table(
  Sample = sample_ids,
  Group = factor(groups, levels = c("Control", "Disease", "Neg_Ctrl")),
  DNA_Conc = dna_conc,
  Kitome_RelAbund = colSums(tss_all[4:5, ]) * 100
)

p_kit <- ggplot(kit_df, aes(x = DNA_Conc, y = Kitome_RelAbund, color = Group, shape = Group)) +
  geom_point(size = 3.6, alpha = 0.9) +
  scale_color_manual(values = c("Control" = "#009E73", "Disease" = "#E69F00", "Neg_Ctrl" = "#CC79A7")) +
  labs(
    title = "B. Negative Control Kitome Dilution Signature",
    subtitle = "Reagent contaminants scale inversely with sample input DNA",
    x = "Extracted DNA Concentration (ng/uL)",
    y = "Kitome Relative Abundance (%)",
    color = "Sample Type",
    shape = "Sample Type"
  ) +
  theme_bw(base_family = "sans", base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

p_comb <- p_comp + p_kit +
  plot_annotation(
    title = "16S Amplicon Smoke Verification: Compositionality & Decontamination Defense",
    theme = theme(plot.title = element_text(face = "bold", size = 13, family = "sans"))
  )

ggsave(
  "repro_results/00_smoke_16s_compositional_qc.png",
  p_comb, width = 12, height = 5.5, dpi = 300, bg = "white"
)

# 6. 写入基础环境版本快照（在运行 04 后会补充完整 16S 包版本） ---------------------
bioc_ver <- tryCatch(as.character(BiocManager::version()), error = function(e) "Unknown")
core_pkgs <- c("dada2", "phyloseq", "microeco", "vegan", "ANCOMBC", "maaslin3",
               "decontam", "ggplot2", "data.table", "patchwork")
ver_rows <- rbindlist(c(
  list(data.table(Package = "R", Version = paste0(R.version$major, ".", R.version$minor), Repository = "Base R")),
  lapply(core_pkgs, function(p) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) "Not Installed")
    repo <- if (p %in% c("dada2", "phyloseq", "ANCOMBC", "maaslin3", "decontam")) paste0("Bioconductor ", bioc_ver) else "CRAN"
    data.table(Package = p, Version = v, Repository = repo)
  })
))
fwrite(ver_rows, "repro_results/environment_versions.tsv", sep = "\t")

cat(sprintf(
  "[SMOKE PASS] Naive_TSS_FP=%d | CLR_FP=%d | True_Recovered=%d/3 | Kitome_Flagged=%d/2\n",
  naive_fp, clr_fp, true_rec, kitome_n
))
