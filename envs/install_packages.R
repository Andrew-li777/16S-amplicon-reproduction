# ==============================================================================
# 16S rRNA 扩增子复现项目：R 环境最小依赖检查与安装脚本
#
# 设计说明:
#   - 采用 "核心双轨组 (core) + 扩展组 (qiime2_bridge / network_extra)" 分层管理，
#     与本项目 R 4.6.0 / Bioconductor 3.23 基线完全对齐。
#   - maaslin3 若在当前 Bioconductor 镜像拉取受阻，自动回退至 biobakery/maaslin3 GitHub 或 Maaslin2。
#
# 用法:
#   Rscript envs/install_packages.R --dry-run   # 仅检查当前缺失包，不执行安装
#   Rscript envs/install_packages.R             # 自动安装缺失的 CRAN & Bioconductor 包
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args

profile <- "core"
if ("--profile" %in% args) {
  idx <- which(args == "--profile")
  if (idx < length(args)) profile <- args[idx + 1]
} else if (any(grepl("^--profile=", args))) {
  profile <- sub("^--profile=", "", args[grep("^--profile=", args)][1])
}

cran_core <- c(
  "BiocManager", "remotes", "yaml", "jsonlite",
  "data.table", "dplyr", "tidyr", "ggplot2", "patchwork",
  "ggrepel", "pheatmap", "RColorBrewer", "scales",
  "vegan", "ape", "phangorn", "picante", "GUniFrac",
  "igraph", "randomForest", "microeco", "file2meco"
)

bioc_core <- c(
  "dada2", "phyloseq", "Biostrings", "DECIPHER",
  "ANCOMBC", "microbiome", "maaslin3", "Maaslin2", "ALDEx2",
  "decontam", "ComplexHeatmap"
)

profile_extras <- switch(
  profile,
  "network_extra" = list(
    cran = c("rgexf", "ggalluvial", "SpiecEasi"),
    bioc = c("mia", "miaViz")
  ),
  list(cran = character(0), bioc = character(0))
)

target_cran <- unique(c(cran_core, profile_extras$cran))
target_bioc <- unique(c(bioc_core, profile_extras$bioc))

installed <- rownames(installed.packages())
missing_cran <- setdiff(target_cran, installed)
missing_bioc <- setdiff(target_bioc, installed)

cat("=== 16S 扩增子复现 R 环境依赖检查 (Profile:", profile, ") ===\n")
cat("R 版本:", R.version.string, "\n")
if ("BiocManager" %in% installed) {
  cat("Bioconductor 版本:", as.character(BiocManager::version()), "\n")
}
cat("CRAN 目标:", length(target_cran), "个 | 已装:", length(target_cran) - length(missing_cran),
    "| 缺失:", length(missing_cran), "\n")
if (length(missing_cran) > 0) cat("  缺失 CRAN:", paste(missing_cran, collapse = ", "), "\n")

cat("Bioc 目标:", length(target_bioc), "个 | 已装:", length(target_bioc) - length(missing_bioc),
    "| 缺失:", length(missing_bioc), "\n")
if (length(missing_bioc) > 0) cat("  缺失 Bioc:", paste(missing_bioc, collapse = ", "), "\n")

if (dry_run) {
  cat("\n[--dry-run] 仅完成检查，未执行任何包安装。\n")
  quit(status = 0)
}

if (length(missing_cran) > 0) {
  cat("\n>>> 安装缺失 CRAN 包...\n")
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (length(missing_bioc) > 0) {
  cat("\n>>> 安装缺失 Bioconductor 包...\n")
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

# maaslin3 回退保护（若镜像缺包则尝试 GitHub）
if (!("maaslin3" %in% rownames(installed.packages()))) {
  cat("\n>>> 尝试从 GitHub biobakery/maaslin3 回退安装...\n")
  tryCatch({
    remotes::install_github("biobakery/maaslin3", upgrade = "never")
  }, error = function(e) {
    cat("[WARN] maaslin3 GitHub 回退失败，后续脚本将使用 Maaslin2 作为降级替代:", conditionMessage(e), "\n")
  })
}

installed_after <- rownames(installed.packages())
still_missing <- setdiff(c(target_cran, target_bioc), installed_after)
cat("\n=== 依赖安装流程结束 | 剩余缺失:", length(still_missing), "个 ===\n")
if (length(still_missing) > 0) {
  cat("剩余缺失列表:", paste(still_missing, collapse = ", "), "\n")
}
