# 16S 清洗中间产物与检查点说明 (`data/processed/`)

本目录保存由 `repro_scripts/` 各阶段自动生成的中间检查点对象（Checkpoints），既支持 `--mode quick` 跳过上游耗时步骤直接复现全部下游统计，又防止临时明细文件污染轻量可审阅的 `repro_results/` 目录。

---

## 1. 目录文件清单与结构契约

| 文件 / 子目录 | 生成脚本 | 实测大小 / 规模 | 内部对象结构与下游使用契约 |
|---|---|---|---|
| `smoke_16s_feature_audit.tsv` | [00_smoke_test.R](../../repro_scripts/00_smoke_test.R) | `6,570` bytes (`40` 行特征) | 记录 40 条合成 ASV 的真实分类（`True_Bloom` / `True_Differential` / `Invariant_Background` / `Kitome_Contaminant`）及 TSS vs. CLR / 阴性富集审计指标（白名单纳入 Git） |
| `filtered_fastq/` | [01_dada2_fastq_to_phyloseq.R](../../repro_scripts/01_dada2_fastq_to_phyloseq.R) | `40` 个 `.fastq.gz`，合计 `22,770,799` bytes (`21.72 MiB`) | 经 `filterAndTrim(truncLen=c(240,160), maxEE=c(2,2))` 过滤的 20 对样本高质量双端读段 |
| `01_dada2_seqtab_nochim.rds` | [01_dada2_fastq_to_phyloseq.R](../../repro_scripts/01_dada2_fastq_to_phyloseq.R) | `11,372` bytes (`20 × 232` 矩阵) | 包含 19 例小鼠样本 + 1 例 `Mock` 对照的全量去嵌合 ASV 计数矩阵（行名样本、列名原始 DNA 序列） |
| `01_ps_miseq.rds` | [01_dada2_fastq_to_phyloseq.R](../../repro_scripts/01_dada2_fastq_to_phyloseq.R) | `20,043` bytes (`19 × 214` `phyloseq`) | **全流水线核心检查点**：剔除 `Mock`、叶绿体、线粒体及未定门 ASV 后的**未抽平原始计数表** + SILVA v138.2 七级分类 + `refseq` DNA + `phangorn::midpoint()` 定根 GTR+G+I 进化树（`quick` 模式前置必填输入） |
| `02_ps_miseq_rarefied.rds` | [02_decontam_and_diversity.R](../../repro_scripts/02_decontam_and_diversity.R) | `20,042` bytes (`19 × 214` `phyloseq`) | 按最小样本深度 `2,516` reads（`rngseed = 20261125`）均匀抽平后的 `phyloseq` 对象，专供 Alpha/Beta 经典生态距离计算 |
| `03_maaslin3_out/` | [03_daa_ancombc2_maaslin3_lefse.R](../../repro_scripts/03_daa_ancombc2_maaslin3_lefse.R) | `9` 个文件，合计 `59,078` bytes | `maaslin3::maaslin3` 原生写入的明细结果（`all_results.tsv`、`significant_results.tsv`、特征转换表及日志） |

---

## 2. 上下游边界契约警示

- **严禁向 `03_daa_ancombc2_maaslin3_lefse.R` 传入 `02_ps_miseq_rarefied.rds`**：`ANCOM-BC2` 与 `MaAsLin 3` 均需基于**未抽平原始计数表** (`01_ps_miseq.rds`) 估计文库抽样分数偏差与流行率模型。
- 若需从零重新生成本目录全部产物，在仓库根目录运行 `python repro_scripts/run_all.py --mode full`。
