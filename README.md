# scCRISPRra

**S**ingle-**C**ell CRISPR Screening **R**eadout **A**nalysis

单细胞 CRISPR 筛选表型分析工具包（v0.2.0）。输入为同时含有
**sgRNA 捕获计数**（CRISPR Guide assay）与**连续表型评分**（metadata 一列）的
Seurat 对象，输出基因水平的显著性（p_high / p_low / FDR）与可视化。
算法本身 virus-agnostic：表型可以是病毒 RNA 载量、感染荧光强度、
任意基因表达量等任何连续评分。

## 两种并列的 gene 级评价方案

包提供两条从 sgRNA 到基因的完整管线，**输入相同、方向语义相同、
可视化函数通用**，可互为正交验证：

| | 方案 A：shape 分解法 | 方案 B：直接法 RRA |
|---|---|---|
| sgRNA 级统计量 | 携带细胞表型分布 vs 全局分布的 2 阶正交多项式系数（`linear`/`quadratic`） | 携带细胞 normalized rank 的 RRA 聚集检验（`p_high`/`p_low`）+ 平均 rank（`mean_rank`，AUC 型） |
| gene 级聚合 | `run_rra_pipeline(score_cols = c("linear","quadratic"))` | `run_direct_rra_pipeline()`（`mean_rank` 经 RRA 聚合） |
| 优点 | 保留分布形状信息（linear 抓整体偏移、quadratic 抓形状偏移） | 不做分布参数化、计算便宜、对窗口/截断选择更稳健 |
| 主入口 | `summarize_sgrna_shapes()` → `run_rra_pipeline()` | `direct_sgrna_rra()` → `run_direct_rra_pipeline()` |

实践中两方案 top50 重叠基因的方向一致率通常 > 93%（真实数据三窗口验证），
共识基因可信度最高；`mean_rank` 单独作为 score 时与 shape `linear` 的
gene 级 signed Spearman 约 0.54–0.66。

## 分析流程

```
Seurat 对象（sgRNA count assay + 连续表型 metadata 列）
   │
   ├─ ① validate_seurat()         输入预检（assay / 表型列 / sample 列）
   │
   ├─ ② call_sgrna()              Poisson-Gaussian 混合模型 guide calling
   │                               （区分 ambient 空转 UMI 与真实转导，
   │                                 参考 Replogle et al. 2020）
   │
   ├─ ③A summarize_sgrna_shapes() 每 sgRNA 携带细胞表型分布 vs 全局分布
   │   │                          → 2 阶正交多项式系数（linear / quadratic）
   │   └─ run_rra_pipeline()      join library → alpha-RRA → gene 级 p/FDR/score
   │
   ├─ ③B direct_sgrna_rra()       窗口内细胞 normalized rank → "sgRNA 当基因、
   │   │                          细胞当重复"复用 alpha_rra → sgRNA 级 p/FDR/mean_rank
   │   └─ run_direct_rra_pipeline()  mean_rank → alpha-RRA → gene 级 p/FDR/score
   │
   └─ ④ plot_rra_volcano() 等     火山图（两方案通用）/ shape 分布 / 单 sgRNA 直方图
```

## 安装

要求 R >= 4.2.0、Seurat >= 5.0.0；依赖 Matrix、ggplot2、ggrepel、gridExtra、dplyr、tidyr。

```r
## 方式一：开发中直接加载（无需 install；需要 pkgload 包）
pkgload::load_all("<scCRISPRra 包路径>", export_all = TRUE, attach = TRUE)

## 方式二：常规安装
devtools::install("<scCRISPRra 包路径>")
library(scCRISPRra)
```

## 输入要求

### 1. Seurat 对象（核心输入）

| 组件 | 要求 | 参数 |
|------|------|------|
| CRISPR Guide assay | `counts` 为 sgRNA × cell 矩阵，rownames = sgRNA id（cellranger feature reference 中 sgRNA 应填入 `id` 字段） | `guide_assay`，默认 `"CRISPR Guide"` |
| 表型评分列 | metadata 中的**连续数值**列（如 `rsv_sum = log2(病毒 UMI + 1)`） | `phenotype_col`，必传 |
| sample 列 | 可选，metadata 中标识批次/样本的列 | `sample_col` |

注意：Seurat 会把 assay 名中的空格替换为点，因此实际常用 `"CRISPR_Guide"`。
用 `validate_seurat()` 可提前校验以上三项。

从 count 矩阵构建最小可用对象：

```r
seu <- CreateSeuratObject(counts = guide_mat, assay = "CRISPR_Guide")
seu$rsv_sum <- phenotype[colnames(seu)]   # 命名数值向量
```

### 2. sgRNA library（gene 聚合需要）

三列 data.frame：`sgrna_id` / `sequence` / `gene_name`。
可用 `load_sgrna_library(path, sep, col_idx)` 读取，按**列位置**（默认 1/2/3）取列，
兼容 MAGeCK `library.txt`（tab 分隔）与普通 csv。

关键：`sgrna_id` 必须与 guide assay 的 rownames **完全一致**——
`run_rra_pipeline()` 与 `run_direct_rra_pipeline()` 按其内连接，无交集会直接报错。

## 输出说明

### 方案 A：sgRNA 级表（`summarize_sgrna_shapes()` 返回）

每行一个 sgRNA：

| 列 | 含义 |
|------|------|
| `sgrna` | sgRNA id |
| `linear` | 1 阶正交多项式系数：携带细胞表型分布相对全局的整体偏移（正值 = 偏向高值端） |
| `quadratic` | 2 阶系数：分布形状/宽度相对全局的偏移 |
| `freq` | 通过 calling 且落在表型区间内的携带 cell 数 |

另有 attributes 可用 `attr()` 提取：`breaks`、`global_hist`、`phenotype_range`、
`n_cells_in_range`（与 shape 分解所用的 bin 定义完全一致，供下游可视化复用）。

### 方案 B：sgRNA 级表（`direct_sgrna_rra()` 返回）

每行一个 sgRNA：

| 列 | 含义 |
|------|------|
| `sgrna` | sgRNA id |
| `n_cells` | 窗口内携带 cell 数（`min_cells_per_sgrna` 过滤依据） |
| `cells_low` / `cells_high` | 聚集在表型低值端 / 高值端（ratio_rank ≤ alpha）的 cell 数 |
| `p_high` / `p_low` | RRA 聚集检验 p 值（方向语义见下节） |
| `fdr_high` / `fdr_low` | BH 校正后 FDR |
| `mean_rank` | 携带细胞平均 normalized rank（AUC 型统计量，表型越低值越小） |

### gene 级表（两方案结构一致）

`run_rra_pipeline()` 返回 list（按 score 列命名）；`run_direct_rra_pipeline()`
返回 list（`sgrna_level` + `gene_level`）。gene 级表每行一个 gene：
`gene`、`sgRNA_count`、`sgRNA_high_in_alpha`、`sgRNA_low_in_alpha`、
`p_high`、`p_low`、`fdr_high`、`fdr_low`、`score`（该 gene 全部 sgRNA 的
score 均值；方案 B 中为 mean_rank 均值）、`rank_high`、`rank_low`、
`p_two`（双侧）、`alpha_high/low`、`genename_high/low`
（显著基因标注，供火山图 label）。

### ⚠ p_high / p_low 方向语义（易混淆，务必读）

两方案的 rank 均按表型（或 shape 系数）**升序**计算，因此方向语义完全一致：

| 显著端 | rank 低端含义 | 敲除细胞的表型 | 表型 = 病毒载量时的生物学解读 |
|--------|--------------|--------------|------------------------------|
| `p_high` 显著 | sgRNA/系数聚集在 rank 低端 | 富集于表型**低值**端 | 敲除后感染载量降低 → 候选**促病毒宿主因子**（受体 / 入胞相关） |
| `p_low` 显著 | 聚集在 rank 高端 | 富集于表型**高值**端 | 敲除后载量升高 → 候选**限制因子** |

`high/low` 命名沿袭 MAGeCK 惯例，指代 **ratio_rank 的端点**而非表型高低；
火山图标题中的 "high-end / low-end enrichment" 同样对应 p 名。
若自行调用 `alpha_rra()` 并传入降序 rank，则两端语义互换。

## 快速上手

```r
library(scCRISPRra)
set.seed(123)

## 1. 输入：Seurat（guide assay + 表型列）与 library
seu   <- CreateSeuratObject(counts = guide_mat, assay = "CRISPR_Guide")
seu$rsv_sum <- phenotype[colnames(seu)]
lib   <- data.frame(sgrna_id = ..., sequence = ..., gene_name = ...)

## 2. 校验 → guide calling
validate_seurat(seu, phenotype_col = "rsv_sum", guide_assay = "CRISPR_Guide")
bool_mat       <- call_sgrna(seu, guide_assay = "CRISPR_Guide",
                             method = "poisson_gaussian",
                             min_umi = 3L, max_guides_per_cell = 3L)

## 3A. 方案 A：shape 分解 → RRA
sgrna_summary  <- summarize_sgrna_shapes(seu, bool_mat,
                                         phenotype_col = "rsv_sum",
                                         guide_assay = "CRISPR_Guide")
rra_out        <- run_rra_pipeline(sgrna_summary,
                                   score_cols = c("linear", "quadratic"),
                                   library_df = lib)
plot_rra_volcano(rra_out[["linear"]], direction = "high")

## 3B. 方案 B：直接法 RRA（同一 bool_mat，无需 shape 分解）
sgl            <- direct_sgrna_rra(seu, bool_mat, phenotype_col = "rsv_sum")
direct_out     <- run_direct_rra_pipeline(sgl, library_df = lib)
plot_rra_volcano(direct_out$gene_level, direction = "high")
```

完整可运行的端到端示例见 `demo/end_to_end_demo.R`（自包含 mock 数据，两方案
并列演示并交叉验证模拟效应，`Rscript demo/end_to_end_demo.R` 直接运行），
详细讲解见 `vignettes/scCRISPRra.Rmd`。真实数据（cellranger + sgRNA 比对 +
病毒 UMI）的完整驱动脚本见项目根目录 `04.run_scCRISPRra.R`。

## 关键参数与注意事项

- **`min_cells_per_sgrna`**（默认 `20`）：每 sgRNA 至少携带多少 cell 才保留。
  非负整数；与 bin 数（breaks 长度）**无关**，不会随 bin 数量自动调整。
  两方案中含义一致（方案 B 中为窗口内携带细胞数）。
- **`phenotype_range` vs `phenotype_quantile`**：前者为固定有效区间
  （如病毒载量 log2 尺度的 `c(9, 13)`），提供时优先于后者（默认 q01/q99 分位数截断）。
  固定区间有先验依据时优先用 `phenotype_range`。两方案支持同一套窗口参数。
- **`weight_by_umi`**（默认 `TRUE`，仅方案 A）：多 sgRNA 细胞按 guide UMI 占比
  软去卷积。注意：若传入的 Seurat 对象中 guide assay 并非 sgRNA-UMI 结构
  （如 gene-union 表达矩阵），必须设 `FALSE`，否则系数会静默变为 NA。
  方案 B 不使用 UMI 权重（RRA 基于 rank，天然免权重）。
- **`max_guides_per_cell`**：`1L` 严格丢弃全部 multi-sgRNA 细胞；
  `Inf` 全保留；中间值保留 UMI 最高的前 N 个。
- **`min_umi`**（默认 `3L`）：guide UMI 阈值，低于此值不判 positive（cellranger 推荐）。
- **`alpha`**（方案 B sgRNA 级，默认 `0.25`）：ratio_rank 聚集阈值。sgRNA 级的
  "重复"是细胞（数百个），远多于 gene 级的 sgRNA 数（2-10 个），故比 gene 级的
  `rra_alpha`（默认 0.1）适当放宽。

## 函数参考

| 分类 | 函数 | 说明 |
|------|------|------|
| 输入校验 | `validate_seurat` | 检查 guide assay / 表型列 / sample 列 |
| Guide calling | `call_sgrna` | `method`：`poisson_gaussian`（推荐）/ `per_cell` / `simple` / `none` |
| 方案 A：shape 分解 | `compute_breaks` | 表型截断区间内生成等距 breaks |
| | `compute_global_hist` | 全局归一化频率分布 |
| | `decompose_shape` | 单 sgRNA 分布偏移的 2 阶正交多项式分解 |
| | `summarize_sgrna_shapes` | 全部 sgRNA 的 linear/quadratic 系数总表（主入口） |
| 方案 B：直接法 | `direct_sgrna_rra` | sgRNA 级 RRA（细胞当重复）+ mean_rank（主入口） |
| | `run_direct_rra_pipeline` | mean_rank 经 RRA 聚合到 gene（主入口） |
| RRA（共用） | `alpha_rra` / `alpha_rra_v2` | MAGeCK 风格 RRA（v2 带 n_eff 收缩） |
| | `load_sgrna_library` | 读 library（mageck library.txt 风格） |
| | `run_rra_pipeline` | join library → score 列逐列 RRA（方案 A 主入口） |
| 可视化（两方案通用） | `plot_rra_volcano` / `plot_rra_volcano_both` | 单侧 / 双侧火山图 |
| | `plot_shape_score_distribution` | shape 系数分布 + 标记目标 sgRNA |
| | `plot_sgrna_distribution` | 单 sgRNA 携带细胞的表型直方图 |
| | `plot_sgrna_bin_offset` | 各 bin 频率偏移条形图 |

## 测试

```r
devtools::test("<scCRISPRra 包路径>")   # testthat：calling / shape / rra / direct / validate
```

## 算法参考

- Replogle, J. et al. (2020) *Nature Biotechnology* — direct capture Perturb-seq 与
  Poisson-Gaussian guide calling 思路
- Li, W. et al. (2014) *Genome Biology* — MAGeCK RRA 聚合算法

## License

MIT
