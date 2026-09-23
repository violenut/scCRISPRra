## ============================================================================
## scCRISPRra 端到端演示脚本（自包含 mock 数据，无需任何外部输入）
## ============================================================================
##
## 流程与真实数据驱动脚本 04.run_scCRISPRra.R 完全一致：
##   构建 Seurat → validate_seurat → call_sgrna → summarize_sgrna_shapes
##   → run_rra_pipeline → 输出 CSV + 火山图
## 区别：本脚本所有输入均由模拟数据生成，任何机器上可直接运行：
##   Rscript end_to_end_demo.R
##
## ──────────────── 包的最小输入（本 demo 实际演示的两项）─────────────────────
## (1) Seurat 对象
##     - guide assay（demo 中命名 "CRISPR_Guide"；包默认 "CRISPR Guide"）：
##       counts = sgRNA × cell 矩阵，rownames = sgRNA id，数值为 guide UMI
##     - metadata 表型列（demo 中 "rsv_sum"）：连续数值，如 log2(病毒 UMI + 1)
##     - metadata sample 列（可选；demo 中为 G1/G2/G3）
## (2) sgRNA library：三列 data.frame（sgrna_id / sequence / gene_name）；
##     真实项目中用 load_sgrna_library(path, sep) 读取（按列位置 1/2/3 取列，
##     兼容 MAGeCK library.txt）。sgrna_id 必须与 guide assay rownames 一致。
##
## ──────────────── 真实项目的输入文件（对照 04.run_scCRISPRra.R）──────────────
## (1) cellranger raw_feature_bc_matrix/：每批一个目录
##     （matrix.mtx / barcodes.tsv / features.tsv），用 Read10X() 读取；
##     多批合并时给 barcode 加样品前缀（如 G1_）避免冲突
## (2) sgRNA matched-seq csv：列 1 = matched_seq（CBD 16nt + 接头 + sgRNA 序列），
##     列 2 = matched_length，列 3 = matched_read；
##     CBD 与 sgRNA 序列从 matched_seq 拆出，折叠成 sgRNA × cell count 矩阵
## (3) 病毒基因 UMI 表 tsv：列 = gene / cell / count；
##     各病毒基因 UMI 求和后 log2(+1) 得表型列 rsv_sum
## (4) sgRNA library csv：三列 id / sgRNA(20nt 序列) / Gene
##
## ──────────────── 输出（写入 <当前工作目录>/scCRISPRra_demo_output/）──────────
## sgrna_summary.csv          sgRNA 级 shape 系数表
##                            （sgrna/linear/quadratic/freq/mean_delta/frac_*...）
## rra_linear.csv             gene 级 RRA 表
## rra_quadratic.csv          （p_high/p_low/fdr_high/fdr_low/score/rank_*/p_two...）
## volcano_linear_p_high.pdf  p_high 侧火山图（KO→表型低值端富集的基因）
## volcano_linear_p_low.pdf   p_low 侧火山图（KO→表型高值端富集的基因）
## volcano_linear_both.pdf    双侧拼图
## shape_score_dist_linear.pdf  linear 系数分布（虚线标记 Gene01 的一个 sgRNA）
## sgrna_dist_selected.pdf    Gene01 / Gene02 各一条 sgRNA 的携带细胞表型直方图
## bin_offset_selected.pdf    同上 sgRNA 各 bin 的频率偏移条形图
##
## ──────────────── 方向语义（run_rra_pipeline 内 rank 按系数升序）──────────────
##   p_high 显著 → 该 gene 的 sgRNA 聚集在系数低端 → 敲除细胞富集于表型低值端
##                （表型 = 病毒载量时 = 候选促病毒宿主因子）
##   p_low  显著 → 敲除细胞富集于表型高值端（= 候选限制因子）
## ============================================================================

set.seed(123)

## ──────────────────── Step 0: 加载包 ────────────────────
## 优先用已安装的 scCRISPRra；未安装时按脚本位置（<pkg>/demo/ 的上一级）
## pkgload::load_all 加载开发版；也可在 source 前自定义变量 PKG_ROOT
if (!require(scCRISPRra, quietly = TRUE)) {
  candidates <- character(0)
  if (exists("PKG_ROOT", inherits = TRUE)) candidates <- c(candidates, PKG_ROOT)
  arg_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg_file) == 1) {
    p <- normalizePath(sub("^--file=", "", arg_file), mustWork = FALSE)
    candidates <- c(candidates, dirname(dirname(p)))   # <pkg>/demo/x.R → <pkg>
  }
  pkg_root <- candidates[file.exists(file.path(candidates, "DESCRIPTION"))][1]
  if (is.na(pkg_root)) {
    stop("未找到 scCRISPRra 包：请 install 后 library(scCRISPRra)，",
         "或用 Rscript 运行本脚本，或先定义 PKG_ROOT <- '<包绝对路径>'")
  }
  pkgload::load_all(pkg_root, export_all = TRUE, attach = TRUE)
}
library(ggplot2)   # ggsave 及绘图函数依赖
out_dir <- file.path(getwd(), "scCRISPRra_demo_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ──────────────────── Step 1: 模拟输入数据 ────────────────────
## 模拟 3 个 sample × 1000 cell；60 条 sgRNA 打 12 个基因（每基因 5 条）。
## 其中 Gene01 = 模拟"促病毒宿主因子"（敲除 → 载量下降），
##      Gene02 = 模拟"限制因子"（敲除 → 载量上升），其余为无效应对照。
## 表型 rsv_sum 模拟双峰感染分布（低载量峰 ~9.2 / 高载量峰 ~12.2）。
n_per_sample <- 1000
samples <- c("G1", "G2", "G3")
n_cells <- n_per_sample * length(samples)
cell_ids <- unlist(lapply(samples, function(g)
  paste0(g, "_", sprintf("bc%05d", seq_len(n_per_sample)), "-1")))
sample_of <- rep(samples, each = n_per_sample)

n_guides <- 60
guide_ids <- paste0("sg", sprintf("%02d", seq_len(n_guides)))
genes <- paste0("Gene", sprintf("%02d", 1:12))
gene_map <- rep(genes, each = n_guides / length(genes))  # sg01-05→Gene01 ...

## 每细胞分配 1 条主 sgRNA；8% 细胞额外携带 1 条（模拟 doublet）
g1_idx <- sample(n_guides, n_cells, replace = TRUE)
is_doublet <- runif(n_cells) < 0.08
g2_idx <- sample(n_guides, n_cells, replace = TRUE)

## guide UMI count 矩阵（sgRNA × cell）：
## 全局 ambient ~ Poisson(0.15)（模拟空转/环境 UMI）
## 携带 guide    ~ Poisson(0.15) + Poisson(30)（模拟真实转导）
counts <- matrix(rpois(n_guides * n_cells, 0.15),
                 nrow = n_guides, ncol = n_cells,
                 dimnames = list(guide_ids, cell_ids))
counts[cbind(g1_idx, seq_len(n_cells))] <-
  counts[cbind(g1_idx, seq_len(n_cells))] + rpois(n_cells, 30)
db <- which(is_doublet)
counts[cbind(g2_idx[db], db)] <-
  counts[cbind(g2_idx[db], db)] + rpois(length(db), 30)

## 表型：双峰基底 + 基因敲除效应
base <- ifelse(runif(n_cells) < 0.6,
               rnorm(n_cells, 9.2, 0.6),
               rnorm(n_cells, 12.2, 0.8))
shift_of_gene <- setNames(c(-1.2, 1.2, rep(0, 10)), genes)  # Gene01/Gene02 有效应
shift <- shift_of_gene[gene_map[g1_idx]] +
  ifelse(is_doublet, shift_of_gene[gene_map[g2_idx]], 0)
phenotype <- setNames(base + shift, cell_ids)

## ---- 按"包的最小输入"组装：Seurat + guide assay + 表型列 + sample 列 ----
seu <- CreateSeuratObject(
  counts = counts,
  assay = "CRISPR_Guide",
  meta.data = data.frame(sample = sample_of, rsv_sum = phenotype,
                         row.names = cell_ids))
cat("Seurat:", ncol(seu), "cell ×", nrow(seu), "sgRNA；assays:",
    paste(Assays(seu), collapse = ", "), "\n")

## sgRNA library：三列 sgrna_id / sequence / gene_name
## （真实项目：lib <- load_sgrna_library("CROP_sglibrary.csv", sep = ",")）
library_df <- data.frame(
  sgrna_id  = guide_ids,
  sequence  = replicate(n_guides, paste0(sample(c("A","C","G","T"), 20,
                                                replace = TRUE), collapse = "")),
  gene_name = gene_map,
  stringsAsFactors = FALSE
)

## ──────────────────── Step 2: 校验 + guide calling ────────────────────
validate_seurat(seu, phenotype_col = "rsv_sum",
                guide_assay = "CRISPR_Guide", sample_col = "sample")
cat("validate_seurat: OK\n")

## Poisson-Gaussian per-guide calling；每细胞最多保留 3 条 sgRNA
bool_mat <- call_sgrna(seu, guide_assay = "CRISPR_Guide",
                       method = "poisson_gaussian",
                       min_umi = 3L, max_guides_per_cell = 3L, seed = 123)

guides_per_cell <- colSums(bool_mat)
cat("calling 阳性总数:", sum(bool_mat),
    "；携带 >=1 sgRNA 的细胞:", sum(guides_per_cell > 0), "/", n_cells, "\n")
cat("每细胞 sgRNA 数分布:\n"); print(table(guides_per_cell))

## ──────────────────── Step 3: shape 分解 ────────────────────
## phenotype_range = NULL → 用 q01/q99 分位数截断；
## 真实项目可传固定区间（如 RSV 载量 log2 尺度的 c(9, 13)），
## 提供时忽略分位数、breaks 直接从区间端点生成。
## demo 用 break_step = 0.25（bin 少、小样本更平滑）；真实项目用 0.01。
sgrna_summary <- summarize_sgrna_shapes(seu, bool_mat,
                                       phenotype_col = "rsv_sum",
                                       phenotype_range = NULL,
                                       break_step = 0.25,
                                       poly_degree = 2L,
                                       guide_assay = "CRISPR_Guide")
cat("shape 分解：", nrow(sgrna_summary), "条 sgRNA 通过 min_cells_per_sgrna=20 过滤\n")

## 与 shape 分解完全一致的 bin / 区间信息（attributes，供可视化复用）
used_range <- attr(sgrna_summary, "phenotype_range")
global_hist <- attr(sgrna_summary, "global_hist")
cat(sprintf("表型区间 [%.2f, %.2f]，%d bins，区间内 %d 细胞\n",
            used_range[1], used_range[2], nrow(global_hist),
            attr(sgrna_summary, "n_cells_in_range")))

## ──────────────────── Step 4: RRA（sgRNA → gene 聚合）────────────────────
rra_out <- run_rra_pipeline(sgrna_summary,
                            score_cols = c("linear", "quadratic"),
                            library_df = library_df,
                            rra_alpha = 0.1, rra_variant = "v1")

cat("\n===== RRA 结果（linear，p_high 侧 top5：敲除→表型低值端）=====\n")
print(head(rra_out$linear[order(rra_out$linear$p_high),
                          c("gene","sgRNA_count","p_high","fdr_high","score")]))
cat("\n===== RRA 结果（linear，p_low 侧 top5：敲除→表型高值端）=====\n")
print(head(rra_out$linear[order(rra_out$linear$p_low),
                          c("gene","sgRNA_count","p_low","fdr_low","score")]))

## 模拟效应验证：Gene01 应 p_high 显著；Gene02 应 p_low 显著
cat("\n===== 模拟效应验证 =====\n")
print(rra_out$linear[rra_out$linear$gene %in% c("Gene01","Gene02"),
                     c("gene","p_high","p_low","fdr_high","fdr_low","score")])
cat("[预期] Gene01（模拟促病毒因子，KO→低载量）→ p_high 显著；\n")
cat("[预期] Gene02（模拟限制因子，  KO→高载量）→ p_low 显著\n")

## ──────────────────── Step 5: 输出表格 ────────────────────
write.csv(sgrna_summary, file.path(out_dir, "sgrna_summary.csv"), row.names = FALSE)
for (col in names(rra_out)) {
  write.csv(rra_out[[col]], file.path(out_dir, paste0("rra_", col, ".csv")),
            row.names = FALSE)
}
cat("\n已写出: sgrna_summary.csv, rra_linear.csv, rra_quadratic.csv\n")

## ──────────────────── Step 6: 可视化 ────────────────────
## 单图失败不影响其余（与 04.run_scCRISPRra.R 一致的写法）
save_pdf <- function(name, p, w, h) {
  tryCatch({
    ggsave(file.path(out_dir, name), p, width = w, height = h, device = "pdf")
    cat("wrote", name, "\n")
  }, error = function(e) message("  [skip] ", name, ": ", conditionMessage(e)))
}

## 6.1 火山图（linear 的 high/low 两侧 + 双侧拼图）
save_pdf("volcano_linear_p_high.pdf",
         plot_rra_volcano(rra_out$linear, direction = "high", ylim_max = 6), 8, 6)
save_pdf("volcano_linear_p_low.pdf",
         plot_rra_volcano(rra_out$linear, direction = "low", ylim_max = 6), 8, 6)
save_pdf("volcano_linear_both.pdf",
         plot_rra_volcano_both(rra_out$linear, ylim_max = 6), 12, 6)

## 6.2 shape 系数分布（标记 Gene01 的一条 sgRNA）
g_gene01 <- guide_ids[which(gene_map == "Gene01")[1]]
g_gene02 <- guide_ids[which(gene_map == "Gene02")[1]]
save_pdf("shape_score_dist_linear.pdf",
         plot_shape_score_distribution(sgrna_summary, score_col = "linear",
                                       target_sgrna = g_gene01,
                                       binwidth = 0.005), 7, 5)

## 6.3 示例 sgRNA 的携带细胞表型分布 + 各 bin 频率偏移
save_pdf("sgrna_dist_selected.pdf",
         plot_sgrna_distribution(phenotype, bool_mat,
                                 target_sgrnas = c(g_gene01, g_gene02),
                                 binwidth = 0.4), 8, 5)
save_pdf("bin_offset_selected.pdf",
         plot_sgrna_bin_offset(global_hist, phenotype, bool_mat,
                               target_sgrnas = c(g_gene01, g_gene02)), 10, 5)

cat("\n[完成] 输出目录:", out_dir, "\n")
cat("表: sgrna_summary.csv / rra_linear.csv / rra_quadratic.csv\n")
cat("图: volcano_linear_{p_high,p_low,both}.pdf, shape_score_dist_linear.pdf,",
    "sgrna_dist_selected.pdf, bin_offset_selected.pdf\n")
