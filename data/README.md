# 16S 扩增子复现工程数据总台账 (`data/`)

本目录严格践行**“原始只读输入 (`raw/`) 与可再生中间产物 (`processed/`) 物理隔离”**的生信工程准则，共托管 Track B 上游原始 FASTQ 测序队列、SILVA 参考库及 Track A 下游多算法中间容器对象。

---

## 1. 数据资产全景汇总

| 层次 | 路径 | 来源 / 生成方式 | 规模与实测体量 | Git 跟踪规则 |
|---|---|---|---|---|
| **原始测序输入** | [`raw/`](./raw/README.md) | Mothur `miseqsopdata.zip` + [Zenodo 14169026](https://zenodo.org/records/14169026) SILVA v138.2 | `45` 个 `MiSeq_SOP` 文件 (`163.09 MiB`) + `zip` (`35.17 MiB`) + SILVA 库 (`200.17 MiB`)，总计 `~398.4 MiB` | 大文件全屏蔽，仅放行 `README.md` 与 `.gitkeep` |
| **过滤后读段** | [`processed/filtered_fastq/`](./processed/README.md) | `01_dada2_fastq_to_phyloseq.R` (`filterAndTrim`) | `40` 个压缩配对 FASTQ（`20` 样本 × R1/R2），合计 `22,770,799` bytes (`21.72 MiB`) | 全屏蔽（可由 `raw/` 重新生成） |
| **核心中间对象** | [`processed/01_ps_miseq.rds`](./processed/README.md) | `01_dada2_fastq_to_phyloseq.R` 封装 | `20,043` bytes (`19.57 KiB`)；`19` 样本 × `214` ASV（含 SILVA 注释、DNA 序列与中点定根 GTR+G+I 进化树） | 屏蔽 `.rds` 二进制，由 `--mode full` 重建并支撑 `--mode quick` 秒级复用 |
| **抽平对象** | [`processed/02_ps_miseq_rarefied.rds`](./processed/README.md) | `02_decontam_and_diversity.R` (`rarefy_even_depth`) | `20,042` bytes (`19.57 KiB`)；统一抽平至最小深度 `2,516` reads/样本 | 屏蔽 `.rds` 二进制 |
| **MaAsLin 3 产出** | [`processed/03_maaslin3_out/`](./processed/README.md) | `03_daa_ancombc2_maaslin3_lefse.R` | `9` 个明细 TSV/Log，合计 `59,078` bytes (`57.69 KiB`) | 屏蔽临时目录，精选交叉汇总写入 `repro_results/` |
| **冒烟审计表** | [`processed/smoke_16s_feature_audit.tsv`](./processed/README.md) | `00_smoke_test.R` | `6,570` bytes (`40` ASV 合成真值 vs TSS/CLR/Kitome 审计) | 白名单放行纳入 Git |

---

## 2. 实验设计与样本元数据契约

1. **小鼠肠道时间序列队列（`MiSeq_SOP`，19 例雌性 C57BL/6 小鼠粪便 V4 测序）**：
   - **`Early` 断奶早期组（`n = 9`，`Day 0`–`Day 9`）**：`F3D0, F3D1, F3D2, F3D3, F3D5, F3D6, F3D7, F3D8, F3D9`（个体差异与定植波动大）。
   - **`Late` 断奶晚期稳定组（`n = 10`，`Day 141`–`Day 150`）**：`F3D141`–`F3D150`。
2. **阳性标准品对照（`Mock`，1 例）**：
   - 配套 `HMP_MOCK.v35.fasta`（含 `32` 条已知菌株 V3–V5 参考序列片段，对应 20 种标准微生物群落），在 `01` 阶段精准用于检验 DADA2 降噪假阳性率（实测检出 `20/20` 精确匹配、`0` 假阳性）。
3. **内置低生物量去污队列（`decontam::MUClite.rds`）** 与 **湿地生态网络队列（`microeco::dataset` + `env_data_16S`）**：
   - 随相应 CRAN/Bioconductor R 包发行并在运行时自动加载，详见 [raw/README.md](./raw/README.md)。
