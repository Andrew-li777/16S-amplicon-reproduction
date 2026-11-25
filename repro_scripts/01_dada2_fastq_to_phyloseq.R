# ==============================================================================
# 01_dada2_fastq_to_phyloseq.R — Track B 上游端到端 16S FASTQ 处理流水线
#
# 流程阶段:
#   1. 自动拉取与校验 Mothur MiSeq_SOP 双端 FASTQ + HMP Mock + SILVA v138.2 训练集
#   2. 绘制 R1/R2 质量分布并执行 Windows 兼容 filterAndTrim (truncLen = 240/160)
#   3. 学习测序错误率模型 -> DADA2 单碱基降噪 -> 双端拼接 -> 一致性去嵌合体
#   4. HMP Mock 阳性对照 20 菌株真值核验 (零假阳性检验)
#   5. SILVA v138.2 门-属-种七级注释 + 过滤叶绿体/线粒体/未定门
#   6. DECIPHER 多序列比对 + phangorn GTR+G+I 系统发育建树 + 中点定根 (Midpoint Root)
#   7. 封装并缓存标准 phyloseq 对象供下游 quick/full 模式共用
# ==============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
  library(Biostrings)
  library(DECIPHER)
  library(phangorn)
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

options(timeout = max(600L, getOption("timeout")))
SEED <- 20261125L
IS_WINDOWS <- (.Platform$OS.type == "windows")

dir.create("data/raw/silva_v138.2", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/filtered_fastq", recursive = TRUE, showWarnings = FALSE)
dir.create("repro_results", recursive = TRUE, showWarnings = FALSE)

safe_download <- function(urls, dest_path, min_bytes = 1e6) {
  if (file.exists(dest_path) && file.info(dest_path)$size >= min_bytes) {
    cat(sprintf("[CACHE HIT] %s (%.2f MB)\n", dest_path, file.info(dest_path)$size / 1048576))
    return(invisible(dest_path))
  }
  for (u in urls) {
    cat(sprintf("[DOWNLOAD] %s -> %s\n", u, dest_path))
    ok <- tryCatch({
      status <- download.file(u, destfile = dest_path, mode = "wb", quiet = FALSE)
      identical(as.integer(status), 0L) && file.exists(dest_path) && file.info(dest_path)$size >= min_bytes
    }, error = function(e) {
      cat("[WARN] 下载失败:", conditionMessage(e), "\n")
      FALSE
    })
    if (ok) {
      cat(sprintf("[OK] 已缓存 %s (%.2f MB)\n", dest_path, file.info(dest_path)$size / 1048576))
      return(invisible(dest_path))
    }
    if (file.exists(dest_path)) unlink(dest_path)
  }
  stop("所有备用下载源均失败: ", dest_path)
}

# 1. 获取并校验 MiSeq_SOP 原始数据与 SILVA v138.2 参考数据库 ----------------------
miseq_zip <- "data/raw/miseqsopdata.zip"
safe_download(
  urls = c("https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"),
  dest_path = miseq_zip,
  min_bytes = 5e6
)

existing_r1 <- list.files("data/raw", pattern = "_R1_001\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
if (length(existing_r1) < 20L) {
  cat("[UNZIP] 解压 miseqsopdata.zip 至 data/raw/ ...\n")
  unzip(miseq_zip, exdir = "data/raw")
  existing_r1 <- list.files("data/raw", pattern = "_R1_001\\.fastq(\\.gz)?$", recursive = TRUE, full.names = TRUE)
}
stopifnot(length(existing_r1) >= 20L)
miseq_dir <- dirname(existing_r1[1])
cat("[DATA] MiSeq_SOP FASTQ 目录:", miseq_dir, "| 样本对数:", length(existing_r1), "\n")

silva_genus_fa <- "data/raw/silva_v138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz"
silva_species_fa <- "data/raw/silva_v138.2/silva_v138.2_assignSpecies.fa.gz"

safe_download(
  urls = c(
    "https://zenodo.org/api/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz/content",
    "https://zenodo.org/api/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz/content"
  ),
  dest_path = silva_genus_fa,
  min_bytes = 5e7
)

safe_download(
  urls = c(
    "https://zenodo.org/api/records/14169026/files/silva_v138.2_assignSpecies.fa.gz/content",
    "https://zenodo.org/api/records/4587955/files/silva_species_assignment_v138.1.fa.gz/content"
  ),
  dest_path = silva_species_fa,
  min_bytes = 2e7
)

# 2. 质量分布诊断与双端截断过滤 (Windows 下 filterAndTrim 必须单线程) -----------
fnFs <- sort(list.files(miseq_dir, pattern = "_R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(miseq_dir, pattern = "_R2_001\\.fastq(\\.gz)?$", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
stopifnot(length(fnFs) == length(fnRs), "Mock" %in% sample.names)

p_qF <- plotQualityProfile(fnFs[1:2]) +
  geom_vline(xintercept = 240, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  labs(title = "Forward Reads (R1) Quality & truncLen=240") +
  theme_bw(base_family = "sans")
p_qR <- plotQualityProfile(fnRs[1:2]) +
  geom_vline(xintercept = 160, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  labs(title = "Reverse Reads (R2) Quality & truncLen=160") +
  theme_bw(base_family = "sans")
ggsave(
  "repro_results/01_fastq_quality_profiles.png",
  p_qF / p_qR, width = 10, height = 8, dpi = 300, bg = "white"
)

filtFs <- file.path("data/processed/filtered_fastq", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("data/processed/filtered_fastq", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

cat("[STEP 2] 执行 filterAndTrim (multithread =", !IS_WINDOWS, ")...\n")
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2),
  truncQ = 2, rm.phix = TRUE, compress = TRUE,
  multithread = !IS_WINDOWS
)

# 3. 测序错误率学习、DADA2 降噪、双端拼接与去嵌合体 ----------------------------
cat("[STEP 3] 学习 R1/R2 测序错误率模型...\n")
set.seed(SEED)
errF <- learnErrors(filtFs, multithread = TRUE, verbose = FALSE)
set.seed(SEED)
errR <- learnErrors(filtRs, multithread = TRUE, verbose = FALSE)

p_errF <- plotErrors(errF, nominalQ = TRUE) +
  labs(title = "A. Forward (R1) Error Model Fit") +
  theme_bw(base_family = "sans", base_size = 9)
p_errR <- plotErrors(errR, nominalQ = TRUE) +
  labs(title = "B. Reverse (R2) Error Model Fit") +
  theme_bw(base_family = "sans", base_size = 9)
ggsave(
  "repro_results/01_dada2_error_models.png",
  p_errF | p_errR, width = 14, height = 6.5, dpi = 300, bg = "white"
)

cat("[STEP 4] 执行 DADA2 单碱基降噪与双端拼接...\n")
dadaFs <- dada(filtFs, err = errF, multithread = TRUE, verbose = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, verbose = FALSE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 12, verbose = FALSE)

seqtab <- makeSequenceTable(mergers)
cat(sprintf("[INFO] 原始拼接 ASV 表维度: %d 样本 × %d 序列\n", nrow(seqtab), ncol(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = FALSE)
nonchim_ratio <- sum(seqtab.nochim) / sum(seqtab)
cat(sprintf("[INFO] 去嵌合后 ASV 表维度: %d 样本 × %d ASV | 读段保留比=%.4f\n",
            nrow(seqtab.nochim), ncol(seqtab.nochim), nonchim_ratio))

getN <- function(x) sum(getUniques(x))
track_dt <- data.table(
  Sample = sample.names,
  input = out[, 1],
  filtered = out[, 2],
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)
track_dt[, pct_retained := round(nonchim / input * 100, 2)]
fwrite(track_dt, "repro_results/01_dada2_read_tracking.tsv", sep = "\t")

track_long <- melt(
  track_dt,
  id.vars = "Sample",
  measure.vars = c("input", "filtered", "denoisedF", "merged", "nonchim"),
  variable.name = "Stage", value.name = "Reads",
  variable.factor = FALSE
)
track_long[, Reads := as.numeric(Reads)]
track_long[, Stage := factor(Stage, levels = c("input", "filtered", "denoisedF", "merged", "nonchim"))]
p_track <- ggplot(track_long, aes(x = Stage, y = Reads, group = Sample, color = Sample == "Mock")) +
  geom_line(alpha = 0.75, linewidth = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "#2C3E50", "TRUE" = "#D55E00"), labels = c("Mouse Gut (n=19)", "HMP Mock Control")) +
  labs(
    title = "DADA2 Pipeline Read Retention Across 20 Paired-End MiSeq Samples",
    subtitle = sprintf("Overall non-chimeric read retention: %.2f%% | Total unique ASVs: %d",
                       nonchim_ratio * 100, ncol(seqtab.nochim)),
    x = "Processing Stage", y = "Read Count", color = "Sample Type"
  ) +
  theme_bw(base_family = "sans", base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(
  "repro_results/01_dada2_read_retention.png",
  p_track, width = 9.5, height = 5.5, dpi = 300, bg = "white"
)

# 4. HMP Mock 阳性对照 20 菌株真值核验 -----------------------------------------
mock_fasta <- list.files(miseq_dir, pattern = "HMP_MOCK.*\\.fasta$", full.names = TRUE)[1]
stopifnot(!is.na(mock_fasta), file.exists(mock_fasta))
mock_ref <- getSequences(mock_fasta)
unqs_mock <- seqtab.nochim["Mock", ]
unqs_mock <- sort(unqs_mock[unqs_mock > 0], decreasing = TRUE)
mock_matches <- sum(sapply(names(unqs_mock), function(x) any(grepl(x, mock_ref, fixed = TRUE))))
mock_fp <- length(unqs_mock) - mock_matches

mock_summary <- data.table(
  Metric = c(
    "Mock_Reference_Strains",
    "Mock_Inferred_ASVs",
    "Mock_Exact_Matches",
    "Mock_False_Positives",
    "Mock_Total_Reads",
    "Nonchim_Read_Ratio",
    "Total_ASVs_Before_Filter"
  ),
  Value = c(
    as.character(length(mock_ref)),
    as.character(length(unqs_mock)),
    as.character(mock_matches),
    as.character(mock_fp),
    as.character(sum(unqs_mock)),
    sprintf("%.5f", nonchim_ratio),
    as.character(ncol(seqtab.nochim))
  )
)
fwrite(mock_summary, "repro_results/01_mock_accuracy_summary.tsv", sep = "\t")
cat(sprintf("[MOCK CHECK] 检出 ASV=%d | 精确匹配参考菌=%d | 假阳性=%d\n",
            length(unqs_mock), mock_matches, mock_fp))

# 5. SILVA v138.2 注释、非细菌剔除与有根进化树构建 -------------------------------
cat("[STEP 5] 运行 SILVA v138.2 物种注释 (assignTaxonomy + addSpecies)...\n")
set.seed(SEED)
taxa <- assignTaxonomy(seqtab.nochim, silva_genus_fa, multithread = TRUE)
set.seed(SEED)
taxa <- addSpecies(taxa, silva_species_fa)

# 剔除 Mock 对照保留 19 个小鼠肠道样本
mouse_samples <- setdiff(rownames(seqtab.nochim), "Mock")
seqtab_mouse <- seqtab.nochim[mouse_samples, , drop = FALSE]
keep_nonzero <- colSums(seqtab_mouse) > 0
seqtab_mouse <- seqtab_mouse[, keep_nonzero, drop = FALSE]
taxa_mouse <- taxa[keep_nonzero, , drop = FALSE]

n_asv_mouse_raw <- ncol(seqtab_mouse)
is_valid_tax <- !(is.na(taxa_mouse[, "Phylum"])) &
  (taxa_mouse[, "Kingdom"] %in% "Bacteria") &
  !(taxa_mouse[, "Order"] %in% "Chloroplast") &
  !(taxa_mouse[, "Family"] %in% "Mitochondria")

seqtab_clean <- seqtab_mouse[, is_valid_tax, drop = FALSE]
taxa_clean <- taxa_mouse[is_valid_tax, , drop = FALSE]
cat(sprintf("[TAX FILTER] 小鼠肠道 ASV: %d -> %d (过滤 %d 条未定门/叶绿体/线粒体)\n",
            n_asv_mouse_raw, ncol(seqtab_clean), n_asv_mouse_raw - ncol(seqtab_clean)))

# 解析实验元信息 (Subject / Gender / Day / When)
samdf <- data.frame(
  SampleID = mouse_samples,
  Subject = sapply(strsplit(mouse_samples, "D"), `[`, 1),
  Gender = substr(mouse_samples, 1, 1),
  Day = as.integer(sapply(strsplit(mouse_samples, "D"), `[`, 2)),
  stringsAsFactors = FALSE
)
samdf$When <- factor(ifelse(samdf$Day > 100, "Late", "Early"), levels = c("Early", "Late"))
rownames(samdf) <- mouse_samples

# DECIPHER 多序列比对 + phangorn GTR+G+I 建树 + Midpoint 定根
cat("[STEP 6] 构建 DECIPHER + phangorn GTR+G+I 中点定根系统发育树...\n")
seqs_clean <- colnames(seqtab_clean)
names(seqs_clean) <- seqs_clean
alignment <- AlignSeqs(DNAStringSet(seqs_clean), anchor = NA, processors = 2, verbose = FALSE)
phang_align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang_align)
treeNJ <- NJ(dm)
fit <- pml(treeNJ, data = phang_align)
fitGTR <- update(fit, k = 4, inv = 0.2)
set.seed(SEED)
fitGTR <- optim.pml(
  fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
  rearrangement = "NNI", control = pml.control(trace = 0)
)
rooted_tree <- phangorn::midpoint(fitGTR$tree)

# 组装标准 phyloseq 对象，保留 DNA 原序列到 refseq 并赋予可读 ASV ID
ps <- phyloseq(
  otu_table(seqtab_clean, taxa_are_rows = FALSE),
  sample_data(samdf),
  tax_table(taxa_clean),
  phy_tree(rooted_tree)
)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- sprintf("ASV%03d", seq_len(ntaxa(ps)))

saveRDS(seqtab.nochim, "data/processed/01_dada2_seqtab_nochim.rds")
saveRDS(ps, "data/processed/01_ps_miseq.rds")

cat(sprintf("[DONE] 已保存 data/processed/01_ps_miseq.rds (%d 样本 × %d ASV)\n",
            nsamples(ps), ntaxa(ps)))
