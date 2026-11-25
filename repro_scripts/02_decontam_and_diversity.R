# ==============================================================================
# 02_decontam_and_diversity.R — Track A-1 阴性对照去污染与群落生态多样性统计防线
#
# 分析内容:
#   Part A (低生物量去污基准): 使用 decontam::isContaminant (Frequency + Prevalence
#          联合模型) 对 MUClite 口腔样本 (539 生物样本 + 30 阴性对照) 进行背景菌识别，
#          绘制 DNA 浓度反比斜率诊断图与阴性对照流行率富集图。
#   Part B (小鼠肠道 16S 多样性防线): 读取 Track B 上游产物 01_ps_miseq.rds (19 样本):
#     1. 绘制门 (Phylum) 与 Top 15 属 (Genus) 相对丰度堆叠图；
#     2. 绘制稀释曲线并在最小深度处抽平，计算 Observed / Chao1 / Shannon / Faith's PD
#        及 Early vs Late Wilcoxon 秩和检验；
#     3. 计算 Bray-Curtis / Weighted UniFrac / Unweighted UniFrac / Aitchison (CLR) 四大
#        距离矩阵，同步运行 PERMANOVA (adonis2, by="margin") 与组内离散度齐性检验
#        (betadisper + permutest) —— 严格防止离散度伪装为组间中心偏移假阳性。
# ==============================================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(decontam)
  library(vegan)
  library(picante)
  library(ape)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
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
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("repro_results", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Part A: decontam 低生物量去污染基准 (MUClite: 539 True + 30 Control)
# ==============================================================================
cat("[PART A] 加载 decontam 内置 MUClite 数据集执行联合去污染检验...\n")
ext_dir <- system.file("extdata", package = "decontam")
ps_muc <- readRDS(file.path(ext_dir, "MUClite.rds"))
ps_muc <- prune_samples(sample_sums(ps_muc) > 0, ps_muc)
sample_data(ps_muc)$is.neg <- sample_data(ps_muc)$Sample_or_Control == "Control Sample"
sample_data(ps_muc)$conc_pos <- pmax(as.numeric(sample_data(ps_muc)$quant_reading), 1)

contam_comb <- isContaminant(
  ps_muc, method = "combined",
  conc = "conc_pos", neg = "is.neg", threshold = 0.1
)
contam_dt <- as.data.table(contam_comb, keep.rownames = "ASV")

# 补充计算每个 ASV 在生物样本 vs 阴性对照中的流行率 (Prevalence)
ps_muc_pa <- transform_sample_counts(ps_muc, function(x) ifelse(x > 0, 1, 0))
pa_neg  <- prune_samples(sample_data(ps_muc_pa)$is.neg == TRUE, ps_muc_pa)
pa_true <- prune_samples(sample_data(ps_muc_pa)$is.neg == FALSE, ps_muc_pa)
contam_dt[, Prev_Neg := taxa_sums(pa_neg) / nsamples(pa_neg)]
contam_dt[, Prev_True := taxa_sums(pa_true) / nsamples(pa_true)]
contam_dt[, Total_Reads := taxa_sums(ps_muc)[ASV]]

n_contam_flagged <- contam_dt[contaminant == TRUE, .N]
stopifnot(n_contam_flagged >= 1L)

fwrite(contam_dt[order(p)], "repro_results/02_decontam_all_taxa.tsv", sep = "\t")

decontam_summary <- data.table(
  Metric = c(
    "MUClite_Total_Samples",
    "MUClite_True_Samples",
    "MUClite_Negative_Controls",
    "MUClite_Total_Taxa",
    "Decontam_Flagged_Contaminants",
    "Decontam_Threshold"
  ),
  Value = c(
    as.character(nsamples(ps_muc)),
    as.character(nsamples(pa_true)),
    as.character(nsamples(pa_neg)),
    as.character(ntaxa(ps_muc)),
    as.character(n_contam_flagged),
    "0.10"
  )
)
fwrite(decontam_summary, "repro_results/02_decontam_summary.tsv", sep = "\t")
cat(sprintf("[PART A] decontam 识别出 %d / %d 个污染特征 (30 阴性对照 vs 539 生物样本)\n",
            n_contam_flagged, ntaxa(ps_muc)))

# 绘制 decontam 双联图: Panel A = 阴性 vs 生物样本流行率对比; Panel B = Top 污染 ASV 浓度反比斜率
p_dec_prev <- ggplot(contam_dt, aes(x = Prev_True, y = Prev_Neg, color = contaminant)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#95A5A6") +
  geom_point(aes(size = log10(Total_Reads + 1)), alpha = 0.78) +
  scale_color_manual(values = c("FALSE" = "#34495E", "TRUE" = "#E74C3C"),
                     labels = c("True Bio Taxa", sprintf("Contaminant (n=%d)", n_contam_flagged))) +
  scale_size_continuous(range = c(1.4, 4.2), guide = "none") +
  labs(
    title = "A. Prevalence in Negative Controls vs True Samples",
    subtitle = "Combined Frequency + Prevalence Test (MUClite, n=569)",
    x = "Prevalence in True Oral Samples (n=539)",
    y = "Prevalence in Negative Controls (n=30)",
    color = "Decontam Call"
  ) +
  theme_bw(base_family = "sans", base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

top_contam_ids <- contam_dt[contaminant == TRUE][order(p)][1:min(4L, .N), ASV]
ps_muc_rel <- transform_sample_counts(ps_muc, function(x) x / sum(x))
muc_meta <- data.frame(sample_data(ps_muc_rel))
freq_df <- rbindlist(lapply(top_contam_ids, function(aid) {
  rel_vals <- as.numeric(otu_table(ps_muc_rel)[, aid])
  if (taxa_are_rows(ps_muc_rel)) rel_vals <- as.numeric(otu_table(ps_muc_rel)[aid, ])
  tax_lbl <- as.character(tax_table(ps_muc_rel)[aid, "Genus"])
  if (is.na(tax_lbl) || tax_lbl == "") tax_lbl <- aid else tax_lbl <- paste0(aid, " (", tax_lbl, ")")
  data.table(
    Taxon = tax_lbl,
    DNA_Conc = muc_meta$conc_pos,
    RelAbund = rel_vals,
    Is_Neg = ifelse(muc_meta$is.neg, "Negative Control", "True Sample")
  )
}))
freq_df <- freq_df[RelAbund > 0]

p_dec_freq <- ggplot(freq_df, aes(x = DNA_Conc, y = RelAbund)) +
  geom_point(aes(color = Is_Neg), alpha = 0.65, size = 1.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#C0392B", linewidth = 0.8) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("True Sample" = "#2C3E50", "Negative Control" = "#D55E00")) +
  facet_wrap(~Taxon, scales = "free_y", ncol = 2) +
  labs(
    title = "B. Inverse Frequency-DNA Concentration Scaling of Top Contaminants",
    subtitle = "Reagent contaminants display log-log negative slope against sample DNA",
    x = "Sample DNA Quantification (log10)",
    y = "Relative Abundance (log10)",
    color = "Sample Type"
  ) +
  theme_bw(base_family = "sans", base_size = 10) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  "repro_results/02_decontam_contaminant_qc.png",
  p_dec_prev + p_dec_freq + plot_layout(widths = c(1, 1.35)),
  width = 13.5, height = 5.6, dpi = 300, bg = "white"
)

# ==============================================================================
# Part B: 小鼠肠道 MiSeq 队列组成、α/β 多样性与 PERMANOVA/betadisper 防线
# ==============================================================================
cat("[PART B] 读取 data/processed/01_ps_miseq.rds 执行生态多样性全分析...\n")
ps <- readRDS("data/processed/01_ps_miseq.rds")
stopifnot(nsamples(ps) == 19L)

meta_df <- data.frame(sample_data(ps))
meta_df$SampleID <- rownames(meta_df)
meta_df <- meta_df[order(meta_df$Day, meta_df$SampleID), ]

# 1. 物种组成堆叠柱状图 (Phylum & Top 12 Genus)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_phy <- tax_glom(ps_rel, taxrank = "Phylum", NArm = FALSE)
phy_melt <- data.table(psmelt(ps_phy))
phy_melt[is.na(Phylum) | Phylum == "", Phylum := "Unclassified"]
phy_top <- phy_melt[, .(MeanAbund = mean(Abundance)), by = Phylum][order(-MeanAbund)][1:min(6L, .N), Phylum]
phy_melt[, Taxon_Plot := ifelse(Phylum %in% phy_top, Phylum, "Other")]
phy_agg <- phy_melt[, .(Abundance = sum(Abundance)), by = .(Sample, When, Day, Taxon_Plot)]
phy_agg[, Sample := factor(Sample, levels = meta_df$SampleID)]

ps_gen <- tax_glom(ps_rel, taxrank = "Genus", NArm = FALSE)
gen_melt <- data.table(psmelt(ps_gen))
gen_melt[is.na(Genus) | Genus == "", Genus := paste0("Unclassified_", Family)]
gen_top <- gen_melt[, .(MeanAbund = mean(Abundance)), by = Genus][order(-MeanAbund)][1:min(12L, .N), Genus]
gen_melt[, Taxon_Plot := ifelse(Genus %in% gen_top, Genus, "Other")]
gen_agg <- gen_melt[, .(Abundance = sum(Abundance)), by = .(Sample, When, Day, Taxon_Plot)]
gen_agg[, Sample := factor(Sample, levels = meta_df$SampleID)]

p_phy <- ggplot(phy_agg, aes(x = Sample, y = Abundance * 100, fill = Taxon_Plot)) +
  geom_col(width = 0.85) +
  facet_grid(~When, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "A. Phylum-Level Composition (Early vs Late Post-Weaning)", x = NULL, y = "Relative Abundance (%)", fill = "Phylum") +
  theme_bw(base_family = "sans", base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

pal_gen <- colorRampPalette(brewer.pal(12, "Paired"))(length(unique(gen_agg$Taxon_Plot)))
p_gen <- ggplot(gen_agg, aes(x = Sample, y = Abundance * 100, fill = Taxon_Plot)) +
  geom_col(width = 0.85) +
  facet_grid(~When, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_gen) +
  labs(title = "B. Top 12 Genera Relative Abundance", x = "Mouse Fecal Sample (Ordered by Day)", y = "Relative Abundance (%)", fill = "Genus") +
  theme_bw(base_family = "sans", base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold"))

ggsave("repro_results/02_taxonomic_composition_bar.png", p_phy / p_gen, width = 12, height = 8.5, dpi = 300, bg = "white")

# 2. 稀释曲线 + 抽平后 Alpha 多样性 (含 Faith's PD)
otu_raw <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) otu_raw <- t(otu_raw)
min_depth <- min(rowSums(otu_raw))

rc_list <- suppressWarnings(vegan::rarecurve(otu_raw, step = 200, tidy = TRUE))
rc_dt <- as.data.table(rc_list)
rc_dt[, When := meta_df[as.character(Site), "When"]]

p_rare <- ggplot(rc_dt, aes(x = Sample, y = Species, group = Site, color = When)) +
  geom_line(linewidth = 0.85, alpha = 0.85) +
  geom_vline(xintercept = min_depth, linetype = "dashed", color = "#2C3E50") +
  annotate("text", x = min_depth + 120, y = max(rc_dt$Species) * 0.15,
           label = sprintf("Min depth = %d", min_depth), hjust = 0, size = 3.2, family = "sans") +
  scale_color_manual(values = c("Early" = "#0072B2", "Late" = "#D55E00")) +
  labs(
    title = "A. ASV Rarefaction Curves",
    subtitle = "Sequencing saturation reached prior to minimum library depth",
    x = "Sequencing Depth (Reads)", y = "Observed ASVs", color = "Stage"
  ) +
  theme_bw(base_family = "sans", base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

set.seed(SEED)
ps_rare <- rarefy_even_depth(ps, sample.size = min_depth, rngseed = SEED, replace = FALSE, verbose = FALSE)
saveRDS(ps_rare, "data/processed/02_ps_miseq_rarefied.rds")

alpha_est <- estimate_richness(ps_rare, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
otu_rare_mat <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_rare_mat <- t(otu_rare_mat)
pd_res <- picante::pd(otu_rare_mat, phy_tree(ps_rare), include.root = TRUE)

alpha_df <- data.table(
  SampleID = rownames(alpha_est),
  When = meta_df[rownames(alpha_est), "When"],
  Day = meta_df[rownames(alpha_est), "Day"],
  Observed = alpha_est$Observed,
  Chao1 = alpha_est$Chao1,
  Shannon = alpha_est$Shannon,
  Faith_PD = pd_res[rownames(alpha_est), "PD"]
)
fwrite(alpha_df, "repro_results/02_alpha_diversity_per_sample.tsv", sep = "\t")

alpha_long <- melt(alpha_df, id.vars = c("SampleID", "When", "Day"),
                   measure.vars = c("Observed", "Chao1", "Shannon", "Faith_PD"),
                   variable.name = "Index", value.name = "Value")
pvals_alpha <- alpha_long[, .(
  P_Wilcox = suppressWarnings(wilcox.test(Value ~ When)$p.value)
), by = Index]
pvals_alpha[, Label := sprintf("Wilcoxon p = %.3f", P_Wilcox)]

p_alpha <- ggplot(alpha_long, aes(x = When, y = Value, fill = When)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.14, size = 2.1, alpha = 0.85) +
  geom_text(data = pvals_alpha, aes(x = 1.5, y = Inf, label = Label),
            vjust = 1.6, size = 3.1, inherit.aes = FALSE, family = "sans") +
  facet_wrap(~Index, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Early" = "#0072B2", "Late" = "#D55E00")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(
    title = "B. Rarefied Alpha Diversity (Early vs Late)",
    subtitle = sprintf("Rarefied at even depth = %d reads/sample", min_depth),
    x = "Time Stage", y = "Diversity Metric", fill = "Stage"
  ) +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  "repro_results/02_alpha_diversity_and_rarefaction.png",
  p_rare + p_alpha + plot_layout(widths = c(1, 1.25)),
  width = 13, height = 5.8, dpi = 300, bg = "white"
)

# 3. 四大 Beta 多样性距离 + PERMANOVA (adonis2) + 离散度齐性校验 (betadisper) ---
dist_bray <- phyloseq::distance(ps_rare, method = "bray")
dist_wuni <- phyloseq::distance(ps_rare, method = "wunifrac")
dist_uuni <- phyloseq::distance(ps_rare, method = "unifrac")

# Aitchison 距离：基于未抽平计数矩阵做 +0.5 伪计数 CLR 变换后算欧氏距离
log_raw <- log(otu_raw + 0.5)
clr_raw <- sweep(log_raw, 1L, rowMeans(log_raw), "-")
dist_ait <- dist(clr_raw, method = "euclidean")

dist_list <- list(
  "Bray-Curtis (Rarefied)"        = dist_bray,
  "Weighted UniFrac (Rarefied)"   = dist_wuni,
  "Unweighted UniFrac (Rarefied)" = dist_uuni,
  "Aitchison CLR (Unrarefied)"    = dist_ait
)

beta_stats <- list()
pcoa_plots <- list()

for (dname in names(dist_list)) {
  dmat <- dist_list[[dname]]
  samp_order <- attr(dmat, "Labels")
  meta_sub <- data.frame(sample_data(ps))[samp_order, , drop = FALSE]

  set.seed(SEED)
  ado <- vegan::adonis2(dmat ~ When, data = meta_sub, permutations = 999L, by = "margin")
  set.seed(SEED)
  bdisp <- vegan::betadisper(dmat, meta_sub$When)
  set.seed(SEED)
  bperm <- vegan::permutest(bdisp, permutations = 999L)

  r2_val <- as.numeric(ado$R2[1])
  f_val  <- as.numeric(ado$F[1])
  p_ado  <- as.numeric(ado$`Pr(>F)`[1])
  p_disp <- as.numeric(bperm$tab$`Pr(>F)`[1])

  beta_stats[[dname]] <- data.table(
    Distance_Metric = dname,
    PERMANOVA_F = round(f_val, 4),
    PERMANOVA_R2 = round(r2_val, 4),
    PERMANOVA_P = p_ado,
    Betadisper_P = p_disp,
    Dispersion_Homogeneous = (p_disp >= 0.05)
  )

  pcoa_res <- ape::pcoa(dmat)
  var_exp <- round(pcoa_res$values$Relative_eig[1:2] * 100, 1)
  coords_df <- data.table(
    SampleID = rownames(pcoa_res$vectors),
    Axis1 = pcoa_res$vectors[, 1],
    Axis2 = pcoa_res$vectors[, 2],
    When = meta_sub$When,
    Day = meta_sub$Day
  )

  ann_txt <- sprintf("PERMANOVA R2=%.3f, p=%.3f\nbetadisper p=%.3f", r2_val, p_ado, p_disp)
  pcoa_plots[[dname]] <- ggplot(coords_df, aes(x = Axis1, y = Axis2, color = When)) +
    stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.7) +
    geom_point(size = 3.2, alpha = 0.9) +
    scale_color_manual(values = c("Early" = "#0072B2", "Late" = "#D55E00")) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.12,
             label = ann_txt, size = 3.0, fill = "#F8F9FA", family = "sans") +
    labs(
      title = dname,
      x = sprintf("PCoA1 (%.1f%%)", var_exp[1]),
      y = sprintf("PCoA2 (%.1f%%)", var_exp[2]),
      color = "Time Stage"
    ) +
    theme_bw(base_family = "sans", base_size = 10.5) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
}

beta_dt <- rbindlist(beta_stats)
fwrite(beta_dt, "repro_results/02_permanova_betadisper_table.tsv", sep = "\t")

p_beta_all <- wrap_plots(pcoa_plots, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "Beta Diversity Ordination & Statistical Defense (PERMANOVA + betadisper)",
    subtitle = "Evaluating community centroid shifts alongside multivariate dispersion homogeneity",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, family = "sans"),
      legend.position = "bottom"
    )
  ) & theme(legend.position = "bottom")

ggsave("repro_results/02_beta_pcoa_4metrics_defense.png", p_beta_all, width = 12.5, height = 9.8, dpi = 300, bg = "white")

# 4. 导出定量摘要指标 ----------------------------------------------------------
bray_row <- beta_dt[Distance_Metric == "Bray-Curtis (Rarefied)"]
wuni_row <- beta_dt[Distance_Metric == "Weighted UniFrac (Rarefied)"]
ait_row  <- beta_dt[Distance_Metric == "Aitchison CLR (Unrarefied)"]

div_summary <- data.table(
  Metric = c(
    "MiSeq_Final_Samples",
    "MiSeq_Clean_ASVs",
    "Rarefaction_Depth",
    "Median_Shannon_Early",
    "Median_Shannon_Late",
    "Shannon_Wilcoxon_P",
    "FaithPD_Wilcoxon_P",
    "Bray_PERMANOVA_R2",
    "Bray_PERMANOVA_P",
    "Bray_Betadisper_P",
    "WeightedUniFrac_PERMANOVA_R2",
    "Aitchison_PERMANOVA_R2"
  ),
  Value = c(
    as.character(nsamples(ps)),
    as.character(ntaxa(ps)),
    as.character(min_depth),
    sprintf("%.4f", median(alpha_df[When == "Early", Shannon])),
    sprintf("%.4f", median(alpha_df[When == "Late", Shannon])),
    sprintf("%.5f", pvals_alpha[Index == "Shannon", P_Wilcox]),
    sprintf("%.5f", pvals_alpha[Index == "Faith_PD", P_Wilcox]),
    sprintf("%.4f", bray_row$PERMANOVA_R2),
    sprintf("%.5f", bray_row$PERMANOVA_P),
    sprintf("%.5f", bray_row$Betadisper_P),
    sprintf("%.4f", wuni_row$PERMANOVA_R2),
    sprintf("%.4f", ait_row$PERMANOVA_R2)
  )
)
fwrite(div_summary, "repro_results/02_diversity_stats_summary.tsv", sep = "\t")
cat("[DONE] 02_decontam_and_diversity.R 执行完成。\n")
