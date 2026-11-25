# 16S 原始只读输入数据说明 (`data/raw/`)

本目录存放 **Track B 上游流水线** 自动拉取并缓存的原始测序压缩包、解压后的双端 FASTQ、Mock 阳性对照参考序列以及 SILVA v138.2 分类参考库。**所有分析脚本严禁原地覆写本目录下的任何文件**；大体积测序与参考库均已被 [`.gitignore`](../../.gitignore) 屏蔽，仅本说明文档纳入版本控制。

---

## 1. 原始文件台账与实测体量

| 目录 / 文件 | 来源与 Accession / URL | 规模与实测大小 | 生物学内容与用途 |
|---|---|---|---|
| `miseqsopdata.zip` | [Mothur MiSeq SOP](https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip) | `36,879,679` bytes (`35.17 MiB`) | Mothur 经典小鼠肠道 16S V4 双端测序教学压缩包 |
| `MiSeq_SOP/` | 由 `miseqsopdata.zip` 解压生成 | `45` 个文件，合计 `171,010,210` bytes (`163.09 MiB`) | 包含 19 只小鼠粪便时间序列配对 FASTQ（`F3D0`–`F3D150`）、1 份 `Mock_S280` 阳性对照双端 FASTQ、20 株参考序列 `HMP_MOCK.v35.fasta`（`17,599` bytes）及原始设计表 |
| `silva_v138.2/silva_nr99_v138.2_toGenus_trainset.fa.gz` | [Zenodo 14169026](https://zenodo.org/records/14169026) | `139,996,892` bytes (`133.51 MiB`) | DADA2 格式化 SILVA nr99 v138.2 界–属（Kingdom–Genus）贝叶斯训练集 |
| `silva_v138.2/silva_v138.2_assignSpecies.fa.gz` | [Zenodo 14169026](https://zenodo.org/records/14169026) | `69,899,323` bytes (`66.66 MiB`) | DADA2 格式化 SILVA v138.2 精确 100% 匹配种水平（Species）补充参考库 |

---

## 2. 独立内置参考数据集（运行时内存加载）

除上述落盘至 `data/raw/` 的大体积原始文件外，下游 **Track A** 脚本还在运行时直接加载以下两个随 R 包发行的标准内置参考数据集：
1. **`MUClite.rds`**（随 Bioconductor `decontam 1.32.0` 发行，`system.file("extdata", "MUClite.rds", package = "decontam")`）：`569` 例人类口腔黏膜炎拭子样本（含 **`539` 例真实样本** + **`30` 例阴性空白对照**）× **`1,951` 个微生物特征**，附带 `quant_reading` DNA 定量浓度列，用于 [02_decontam_and_diversity.R](../../repro_scripts/02_decontam_and_diversity.R) 验证频率 + 流行率联合去污染模型。
2. **`dataset` & `env_data_16S`**（随 CRAN `microeco 2.4.3` 发行）：`90` 例中国北方湿地土壤 16S 样本（`CW`/`IW`/`TW` 三组各 30 例）× **`404` 个特征** + 配套理化因子（pH、Conductivity、Temperature、Precipitation、TOC、TN、NH4、NO3、Altitude 等），用于 [04_microeco_network_and_func.R](../../repro_scripts/04_microeco_network_and_func.R)。

---

## 3. 自动获取与完整性校验约定

首次运行 `--mode full` 时，[01_dada2_fastq_to_phyloseq.R](../../repro_scripts/01_dada2_fastq_to_phyloseq.R) 会自动执行：
- 二进制模式（`mode = "wb"`, `timeout = 600`）下载并校验最小字节数（`miseqsopdata.zip >= 5 MB`、`toGenus >= 50 MB`、`assignSpecies >= 20 MB`）；
- 若 Zenodo `14169026`（SILVA v138.2）临时不可达，自动切换至备用源 [Zenodo 4587955](https://zenodo.org/records/4587955)（SILVA v138.1）；
- 已存在的有效缓存将触发 `[CACHE HIT]` 跳过重复下载。
