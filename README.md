# scCRISPRra

**S**ingle-**C**ell CRISPR Screening **R**eadout **A**nalysis

单细胞 CRISPR 筛选表型分析工具包（v0.1.0）。输入为同时含有
**sgRNA 捕获计数**（CRISPR Guide assay）与**连续表型评分**（metadata 一列）的
Seurat 对象，输出基因水平的显著性（p_high / p_low / FDR）与可视化。
算法本身 virus-agnostic：表型可以是病毒 RNA 载量、感染荧光强度、
任意基因表达量等任何连续评分。

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
   ├─ ③ summarize_sgrna_shapes()  每 sgRNA 携带细胞的表型分布 vs 全局分布
   │                               → 正交多项式 shape 系数（linear / quadratic）
   │                               + 扩展指标（mean_delta / 分位数偏移 / Wasserstein）
   │
   ├─ ④ run_rra_pipeline()        join sgRNA library → MAGeCK 风格 alpha-RRA
   │                               → gene 级 p_high / p_low / FDR / score
   │
   └─ ⑤ plot_rra_volcano() 等     火山图 / shape 系数分布 / 单 sgRNA 表型直方图
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
`run_rra_pipeline()` 按其内连接，无交集会直接报错。

## 输出说明

### sgRNA 级表（`summarize_sgrna_shapes()` 返回）

每行一个 sgRNA：

| 列 | 含义 |
|------|------|
| `sgrna` | sgRNA id |
| `linear` / `quadratic` | 正交多项式 shape 系数：相对全局分布的偏移（正值 = 偏向表型高值端） |
| `freq` | 通过 calling 且落在表型区间内的携带 cell 数 |
| `mean_delta` | 携带 cell 表型均值 − 全局均值 |
| `frac_low` / `frac_high` | 低于全局 q10 / 高于全局 q90 的比例 |
| `q10_dev` … `q90_dev` | 各分位数相对全局的偏移 |
| `wasserstein_dist` | 与全局分布的 Wasserstein-1 距离 |

另有 attributes 可用 `attr()` 提取：`breaks`、`global_hist`、`phenotype_range`、
`n_cells_in_range`（与 shape 分解所用的 bin 定义完全一致，供下游可视化复用）。

### gene 级表（`run_rra_pipeline()` 返回，list 按 score 列命名）

每行一个 gene：`gene`、`sgRNA_count`、`sgRNA_high_in_alpha`、`sgRNA_low_in_alpha`、
`p_high`、`p_low`、`fdr_high`、`fdr_low`、`score`（该 gene 全部 sgRNA 系数均值）、
`rank_high`、`rank_low`、`p_two`（双侧）、`alpha_high/low`、`genename_high/low`
（显著基因标注，供火山图 label）。

### ⚠ p_high / p_low 方向语义（易混淆，务必读）

`run_rra_pipeline()` 内 rank 按 shape 系数**升序**计算：

| 显著端 | sgRNA 聚集位置 | linear 系数 | 敲除细胞的表型 | 表型 = 病毒载量时的生物学解读 |
|--------|--------------|------------|--------------|------------------------------|
| `p_high` 显著 | rank 低端（系数小） | 偏负 | 富集于表型**低值**端 | 敲除后感染载量降低 → 候选**促病毒宿主因子**（受体 / 入胞相关） |
| `p_low` 显著 | rank 高端（系数大） | 偏正 | 富集于表型**高值**端 | 敲除后载量升高 → 候选**限制因子** |

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

## 2. 校验 → guide calling → shape 分解 → RRA
validate_seurat(seu, phenotype_col = "rsv_sum", guide_assay = "CRISPR_Guide")
bool_mat       <- call_sgrna(seu, guide_assay = "CRISPR_Guide",
                             method = "poisson_gaussian",
                             min_umi = 3L, max_guides_per_cell = 3L)
sgrna_summary  <- summarize_sgrna_shapes(seu, bool_mat,
                                         phenotype_col = "rsv_sum",
                                         guide_assay = "CRISPR_Guide")
rra_out        <- run_rra_pipeline(sgrna_summary,
                                  score_cols = c("linear", "quadratic"),
                                  library_df = lib)

## 3. 火山图
plot_rra_volcano(rra_out[["linear"]], direction = "high")
```

完整可运行的端到端示例见 `demo/end_to_end_demo.R`（自包含 mock 数据，
`Rscript demo/end_to_end_demo.R` 直接运行），详细讲解见
`vignettes/scCRISPRra.Rmd`。真实数据（cellranger + sgRNA 比对 + 病毒 UMI）
的完整驱动脚本见项目根目录 `04.run_scCRISPRna.R`。

## 关键参数与注意事项

- **`min_cells_per_sgrna`**（默认 `20`）：每 sgRNA 至少携带多少 cell 才保留。
  非负整数；与 bin 数（breaks 长度）**无关**，不会随 bin 数量自动调整。
- **`phenotype_range` vs `phenotype_quantile`**：前者为固定有效区间
  （如病毒载量 log2 尺度的 `c(9, 13)`），提供时优先于后者（默认 q01/q99 分位数截断）。
  固定区间有先验依据时优先用 `phenotype_range`。
- **`weight_by_umi`**（默认 `TRUE`）：多 sgRNA 细胞按 guide UMI 占比软去卷积。
  注意：若传入的 Seurat 对象中 guide assay 并非 sgRNA-UMI 结构
  （如 gene-union 表达矩阵），必须设 `FALSE`，否则系数会静默变为 NA。
- **`max_guides_per_cell`**：`1L` 严格丢弃全部 multi-sgRNA 细胞；
  `Inf` 全保留；中间值保留 UMI 最高的前 N 个。
- **`min_umi`**（默认 `3L`）：guide UMI 阈值，低于此值不判 positive（cellranger 推荐）。

## 函数参考

| 分类 | 函数 | 说明 |
|------|------|------|
| 输入校验 | `validate_seurat` | 检查 guide assay / 表型列 / sample 列 |
| Guide calling | `call_sgrna` | `method`：`poisson_gaussian`（推荐）/ `per_cell` / `simple` / `none` |
| Shape 分解 | `compute_breaks` | 表型截断区间内生成等距 breaks |
| | `compute_global_hist` | 全局归一化频率分布 |
| | `decompose_shape` | 单 sgRNA 分布偏移的正交多项式分解 |
| | `summarize_sgrna_shapes` | 全部 sgRNA 的 shape 系数总表（主入口） |
| RRA | `alpha_rra` / `alpha_rra_v2` | MAGeCK 风格 RRA（v2 带 n_eff 收缩） |
| | `load_sgrna_library` | 读 library（mageck library.txt 风格） |
| | `run_rra_pipeline` | join library → 分数列逐列 RRA（主入口） |
| 可视化 | `plot_rra_volcano` / `plot_rra_volcano_both` | 单侧 / 双侧火山图 |
| | `plot_shape_score_distribution` | shape 系数分布 + 标记目标 sgRNA |
| | `plot_sgrna_distribution` | 单 sgRNA 携带细胞的表型直方图 |
| | `plot_sgrna_bin_offset` | 各 bin 频率偏移条形图 |

## 测试

```r
devtools::test("<scCRISPRra 包路径>")   # testthat：calling / shape / rra / validate
```

## 算法参考

- Replogle, J. et al. (2020) *Nature Biotechnology* — direct capture Perturb-seq 与
  Poisson-Gaussian guide calling 思路
- Li, W. et al. (2014) *Genome Biology* — MAGeCK RRA 聚合算法

## License

MIT
