# ==============================================================================
# 04_microeco_network_and_func.R — Track A-3 共现网络、Keystone 识别与生态功能/标志物
#
# 核心内容 (对齐 Nature Protocols 2025/2026 microeco 工作流):
#   1. 基于 microeco 内置 90 样本湿地 16S 队列构建 Spearman 共现网络与模块划分;
#   2. 计算节点模块内连通度 (Zi) 与模块间连通度 (Pi), 识别 Keystone 物种
#      (Module hubs / Connectors / Network hubs);
#   3. db-RDA 限制性排序解析理化环境因子 (pH/Conductivity/TN/TOC) 驱动机制;
#   4. FAPROTAX 原核生物群落生态功能预测;
#   5. Random Forest 随机森林多分类器构建与 Top 关键属生物标志物评估;
#   6. 在真实子进程末尾汇出全部 R/Bioconductor/Python 版本锁表与 sessionInfo。
# ==============================================================================

suppressPackageStartupMessages({
  library(microeco)
  library(igraph)
  library(randomForest)
  library(vegan)
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

SEED <- 20261125L
dir.create("repro_results", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

# 1. 加载 microeco 内置 90 样本 16S 队列与理化元数据 ----------------------------
cat("[STEP 1] 加载 microeco::dataset (90 samples x 404 features) 与 env_data_16S...\n")
data(dataset, package = "microeco")
data(env_data_16S, package = "microeco")
mt <- clone(dataset)
mt$tidy_dataset()
mt$cal_abund()

n_samples <- nrow(mt$sample_table)
n_asvs <- nrow(mt$otu_table)
stopifnot(n_samples == 90L)

# 2. 共现网络构建、模块划分与 Zi-Pi Keystone 识别 ------------------------------
cat("[STEP 2] 构建 Spearman 共现网络并计算 Zi-Pi 拓扑角色...\n")
set.seed(SEED)
tn <- trans_network$new(dataset = mt, cor_method = "spearman", filter_thres = 0.001)
tn$cal_network(network_method = "COR", COR_p_thres = 0.01, COR_cut = 0.55)
set.seed(SEED)
tn$cal_module(method = "cluster_fast_greedy")
tn$cal_network_attr()
tn$get_node_table(node_roles = TRUE)

net_attr <- tn$res_network_attr
node_dt <- as.data.table(tn$res_node_table)
fwrite(node_dt, "repro_results/04_network_nodes_zipi.tsv", sep = "\t")

n_nodes <- as.integer(round(net_attr["Vertex", 1]))
n_edges <- as.integer(round(net_attr["Edge", 1]))
modularity_val <- as.numeric(net_attr["Modularity", 1])
n_modules <- length(unique(node_dt$module))
n_mod_hubs <- node_dt[taxa_roles == "Module hubs", .N]
n_connectors <- node_dt[taxa_roles == "Connectors", .N]
n_net_hubs <- node_dt[taxa_roles == "Network hubs", .N]

p_zipi <- tn$plot_taxa_roles(use_type = 1, add_label = TRUE, label_text_size = 3.0) +
  labs(
    title = "A. Topological Keystone Roles (Zi-Pi Analysis)",
    subtitle = sprintf("Nodes=%d | Edges=%d | Modularity=%.3f | Hubs+Connectors=%d",
                       n_nodes, n_edges, modularity_val, n_mod_hubs + n_connectors + n_net_hubs)
  ) +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"))

# 绘制主模块门水平组成与节点度分布
top_mods <- node_dt[, .N, by = module][order(-N)][1:min(6L, .N), module]
mod_sub <- node_dt[module %in% top_mods]
mod_sub[, module := factor(module, levels = top_mods)]
mod_sub[is.na(Phylum) | Phylum == "", Phylum := "Unclassified"]
mod_sub[, Phylum_Clean := sub("^p__", "", Phylum)]
top_phys <- mod_sub[, .N, by = Phylum_Clean][order(-N)][1:min(6L, .N), Phylum_Clean]
mod_sub[, Phylum_Plot := ifelse(Phylum_Clean %in% top_phys, Phylum_Clean, "Other")]

p_mod <- ggplot(mod_sub, aes(x = module, fill = Phylum_Plot)) +
  geom_bar(width = 0.72) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "B. Top Co-Occurrence Modules by Phylum",
    subtitle = sprintf("Total modules=%d (cluster_fast_greedy); COR |r|>=0.55, FDR<0.01", n_modules),
    x = "Network Module", y = "Node Count (ASVs)", fill = "Phylum"
  ) +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

ggsave(
  "repro_results/04_microeco_network_and_zipi.png",
  p_zipi + p_mod + plot_layout(widths = c(1.2, 1)),
  width = 13.5, height = 5.8, dpi = 300, bg = "white"
)

# 3. 环境因子约束排序 (db-RDA) + FAPROTAX 功能 + Random Forest 标志物 -----------
cat("[STEP 3] 执行 db-RDA 环境关联、FAPROTAX 功能预测与 Random Forest 分类...\n")
t_func <- trans_func$new(mt)
t_func$cal_func(prok_database = "FAPROTAX")
t_func$cal_func_FR(abundance_weighted = TRUE)
n_funcs <- ncol(t_func$res_func_FR)
p_func <- t_func$plot_func_FR(add_facet = TRUE) +
  labs(
    title = "A. FAPROTAX Functional Groups (Abundance-Weighted Across 90 Samples)",
    x = "Samples (ordered across CW / IW / TW wetland groups)"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", family = "sans", size = 11)
  )

t_env <- trans_env$new(dataset = mt, add_data = env_data_16S)
set.seed(SEED)
t_env$cal_ordination(method = "dbRDA", use_measure = "bray")
t_env$trans_ordination(show_taxa = 8, adjust_arrow_length = TRUE, max_perc_env = 1.4, max_perc_tax = 1.4)
p_dbrda <- t_env$plot_ordination(plot_color = "Group") +
  labs(title = "B. db-RDA Ordination (Bray-Curtis + Soil Factors)") +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))

# 属水平相对丰度矩阵构建 Random Forest 多分类模型 (CW / IW / TW)
raw_gen_mat <- mt$taxa_abund$Genus
raw_gen_names <- sub("^.*\\|g__", "", rownames(raw_gen_mat))
valid_gen <- nzchar(trimws(raw_gen_names)) &
  !grepl("^(unclassified|uncultured|unknown|NA|g__)$|Incertae_Sedis|Incertae Sedis", raw_gen_names, ignore.case = TRUE)
raw_gen_mat <- raw_gen_mat[valid_gen, , drop = FALSE]
clean_names <- raw_gen_names[valid_gen]
# 合并同名属丰度后再转置为 样本 × 属 矩阵
gen_summed <- rowsum(as.matrix(raw_gen_mat), group = clean_names, reorder = FALSE)
gen_mat <- t(gen_summed)
y_grp <- factor(mt$sample_table[rownames(gen_mat), "Group"])

set.seed(SEED)
rf_fit <- randomForest::randomForest(x = gen_mat, y = y_grp, ntree = 500L, importance = TRUE)
oob_acc <- 1 - tail(rf_fit$err.rate[, "OOB"], 1)

imp_dt <- as.data.table(randomForest::importance(rf_fit), keep.rownames = "Genus")
imp_dt <- imp_dt[order(-MeanDecreaseGini)]
fwrite(imp_dt, "repro_results/04_rf_top_biomarkers.tsv", sep = "\t")

top_rf <- imp_dt[1:min(12L, .N)]
top_rf[, Genus_Label := factor(Genus, levels = rev(Genus))]
p_rf <- ggplot(top_rf, aes(x = MeanDecreaseGini, y = Genus_Label)) +
  geom_segment(aes(x = 0, xend = MeanDecreaseGini, y = Genus_Label, yend = Genus_Label),
               color = "#2C3E50", linewidth = 0.8) +
  geom_point(color = "#009E73", size = 3.3) +
  labs(
    title = "C. Random Forest Genus Biomarkers",
    subtitle = sprintf("500 trees across 3 wetland groups | OOB Accuracy = %.2f%%", oob_acc * 100),
    x = "Mean Decrease Gini", y = NULL
  ) +
  theme_bw(base_family = "sans", base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"))

p_bottom <- p_dbrda + p_rf + plot_layout(widths = c(1.15, 1))
ggsave(
  "repro_results/04_microeco_dbrda_faprotax_rf.png",
  p_func / p_bottom + plot_layout(heights = c(1.25, 1)),
  width = 14, height = 11.5, dpi = 300, bg = "white"
)

# 4. 导出定量摘要指标 ----------------------------------------------------------
summary_04 <- data.table(
  Metric = c(
    "Microeco_Samples",
    "Microeco_ASVs",
    "Network_Nodes",
    "Network_Edges",
    "Network_Modules",
    "Network_Modularity",
    "Keystone_Module_Hubs",
    "Keystone_Connectors",
    "Keystone_Network_Hubs",
    "FAPROTAX_Functions",
    "RF_OOB_Accuracy"
  ),
  Value = c(
    as.character(n_samples),
    as.character(n_asvs),
    as.character(n_nodes),
    as.character(n_edges),
    as.character(n_modules),
    sprintf("%.4f", modularity_val),
    as.character(n_mod_hubs),
    as.character(n_connectors),
    as.character(n_net_hubs),
    as.character(n_funcs),
    sprintf("%.4f", oob_acc)
  )
)
fwrite(summary_04, "repro_results/04_microeco_summary.tsv", sep = "\t")

# 5. 由真实子进程导出全部依赖精准版本锁表与 sessionInfo() -----------------------
bioc_ver <- tryCatch(as.character(BiocManager::version()), error = function(e) "Unknown")
bioc_label <- paste0("Bioconductor ", bioc_ver)

bioc_pkgs <- c("dada2", "phyloseq", "Biostrings", "DECIPHER", "ANCOMBC", "microbiome",
               "maaslin3", "Maaslin2", "ALDEx2", "decontam", "ComplexHeatmap")
cran_pkgs <- c("microeco", "file2meco", "vegan", "ape", "phangorn", "picante",
               "GUniFrac", "igraph", "randomForest", "data.table", "ggplot2",
               "patchwork", "ggrepel", "RColorBrewer")

all_pkgs <- c(cran_pkgs, bioc_pkgs)
ver_rows <- rbindlist(c(
  list(data.table(Package = "R", Version = paste0(R.version$major, ".", R.version$minor), Repository = "Base R")),
  lapply(all_pkgs, function(p) {
    v <- tryCatch(as.character(packageVersion(p)), error = function(e) "Not Installed")
    repo <- if (p %in% bioc_pkgs) bioc_label else "CRAN"
    data.table(Package = p, Version = v, Repository = repo)
  })
))
fwrite(ver_rows, "repro_results/environment_versions.tsv", sep = "\t")

capture.output(sessionInfo(), file = "logs/sessionInfo.txt")
cat(sprintf("[DONE] 04_microeco_network_and_func.R 完成 | RF OOB=%.2f%% | Nodes=%d, Edges=%d\n",
            oob_acc * 100, n_nodes, n_edges))
