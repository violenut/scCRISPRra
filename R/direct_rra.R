#' 直接法 sgRNA 级 RRA（Layer 1：细胞当重复、sgRNA 当基因）
#'
#' 跳过多项式 shape 分解，把 MAGeCK 风格 RRA 直接下推到 sgRNA→细胞层：
#' 在表型区间内对每个细胞算全局 normalized rank（升序：表型低 → rank 小），
#' 然后把"sgRNA 当 gene、携带细胞当重复"复用 \code{\link{alpha_rra}}，
#' 检验每个 sgRNA 携带细胞的 rank 是否聚集在某一端。
#'
#' 与 shape 法的语义对应（两者 rank 均为表型升序）：
#' \itemize{
#'   \item \code{p_high} 显著 = 携带细胞聚集在表型低值端 = 敲除后载量降低 =
#'         候选促病毒宿主因子（对应 shape 法 p_high 显著 = linear 系数偏负）
#'   \item \code{p_low} 显著 = 携带细胞聚集在表型高值端 = 候选限制因子
#' }
#'
#' 输出中的 \code{mean_rank} 是携带细胞平均 normalized rank（AUC 型统计量，
#' 表型越低值越小），作为 Layer 2 基因聚合的 score 列。
#'
#' @param seu Seurat 对象（含表型 metadata 列）或直接的命名数值表型向量
#' @param sgrna_calls boolean matrix（sgRNA × cell），来自 \code{\link{call_sgrna}}
#' @param phenotype_col metadata 中连续表型评分列名（仅在 \code{seu} 为 Seurat
#'   对象时使用），默认 \code{NULL}
#' @param phenotype_range 长度 2 数值向量 \code{c(low, high)} 或 \code{NULL}（默认）。
#'   提供时只保留表型落在 \code{[low, high]} 的 cell
#' @param phenotype_quantile 长度 2 数值向量，默认 \code{c(0.01, 0.99)}。
#'   仅在 \code{phenotype_range = NULL} 时生效
#' @param min_cells_per_sgrna 每 sgRNA 窗口内至少携带多少 cell 才保留，默认 20
#' @param alpha alpha_rra 的 ratio_rank 阈值。sgRNA 级重复数（细胞数）远多于
#'   gene 级（sgRNA 数），适当放宽到 0.25，默认 \code{0.25}
#' @param fdr_method p 值校正方法，默认 \code{"BH"}
#'
#' @return data.frame，每行一个 sgRNA，列：
#'   \code{sgrna}（sgRNA id）、\code{n_cells}（窗口内携带 cell 数）、
#'   \code{cells_low}（聚集在表型低值端 ratio_rank<=alpha 的 cell 数）、
#'   \code{cells_high}（聚集在表型高值端的 cell 数）、
#'   \code{p_high, p_low, fdr_high, fdr_low}、
#'   \code{mean_rank}（携带细胞平均 normalized rank）
#' @export
#'
#' @examples
#' \dontrun{
#' calls <- call_sgrna(seu)
#' sgl <- direct_sgrna_rra(seu, calls, phenotype_col = "virus_load",
#'                         phenotype_range = c(9, 13))
#' }
direct_sgrna_rra <- function(seu,
                             sgrna_calls,
                             phenotype_col = NULL,
                             phenotype_range = NULL,
                             phenotype_quantile = c(0.01, 0.99),
                             min_cells_per_sgrna = 20L,
                             alpha = 0.25,
                             fdr_method = "BH") {
  ## 接受 Seurat 对象（内部提取表型）或直接的命名数值表型向量
  if (inherits(seu, "Seurat")) {
    if (is.null(phenotype_col)) {
      stop("传入 Seurat 对象时必须指定 phenotype_col")
    }
    phenotype <- .get_phenotype_vector(seu, phenotype_col)
  } else {
    phenotype <- seu
  }
  stopifnot("'phenotype' 必须是命名数值向量" =
              is.numeric(phenotype) && !is.null(names(phenotype)))
  stopifnot("'sgrna_calls' 必须是 logical matrix" =
              is.logical(sgrna_calls) && is.matrix(sgrna_calls))
  stopifnot("'min_cells_per_sgrna' 必须是非负整数" =
              is.numeric(min_cells_per_sgrna) && length(min_cells_per_sgrna) == 1 &&
              min_cells_per_sgrna >= 0 &&
              min_cells_per_sgrna == round(min_cells_per_sgrna))
  if (!is.null(phenotype_range)) {
    stopifnot("'phenotype_range' 必须是长度 2 的数值 c(low, high) 且 low < high" =
                is.numeric(phenotype_range) && length(phenotype_range) == 2 &&
                phenotype_range[1] < phenotype_range[2])
  }

  ## 1. 窗口截断（与 summarize_sgrna_shapes 同一套规则）
  if (!is.null(phenotype_range)) {
    lo <- phenotype_range[1]
    hi <- phenotype_range[2]
  } else {
    lo <- stats::quantile(phenotype, probs = phenotype_quantile[1])
    hi <- stats::quantile(phenotype, probs = phenotype_quantile[2])
  }
  pheno_in <- phenotype[phenotype >= lo & phenotype <= hi]
  n_w <- length(pheno_in)
  if (n_w < 2) {
    stop(sprintf("表型区间 [%s, %s] 内仅 %d 个 cell，无法计算",
                 format(lo), format(hi), n_w))
  }

  ## 2. 窗口内 normalized rank（升序：表型低值端 → 小 rank）
  r_win <- setNames(rank(pheno_in) / (n_w + 1), names(pheno_in))

  ## 3. 展开 positive (sgRNA, cell) 对，仅保留窗口内细胞
  idx <- which(sgrna_calls, arr.ind = TRUE)
  gid <- rownames(sgrna_calls)[idx[, "row"]]
  cid <- colnames(sgrna_calls)[idx[, "col"]]
  keep <- cid %in% names(pheno_in)
  df <- data.frame(gene = gid[keep], sgrna = cid[keep],
                   rank = unname(r_win[cid[keep]]),
                   stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    stop("表型区间内没有任何 positive call，请检查窗口或 calling 结果")
  }

  ## 4. min_cells_per_sgrna 过滤
  n_tab <- table(df$gene)
  ok_ids <- names(n_tab)[n_tab >= min_cells_per_sgrna]
  df <- df[df$gene %in% ok_ids, ]
  if (nrow(df) == 0) {
    stop(sprintf("所有 sgRNA 窗口内携带细胞数均 < %d，请降低 min_cells_per_sgrna",
                 min_cells_per_sgrna))
  }

  ## 5. sgRNA 级 RRA：sgRNA 当 gene、细胞当重复
  rr <- alpha_rra(df, N_total = n_w, alpha = alpha, fdr_method = fdr_method)
  out <- data.frame(
    sgrna      = rr$gene,
    n_cells    = rr$sgRNA_count,
    cells_low  = rr$sgRNA_high_in_alpha,   # ratio_rank 小端 = 表型低值端
    cells_high = rr$sgRNA_low_in_alpha,    # ratio_rank 大端 = 表型高值端
    p_high     = rr$p_high,
    p_low      = rr$p_low,
    fdr_high   = rr$fdr_high,
    fdr_low    = rr$fdr_low,
    stringsAsFactors = FALSE
  )

  ## 6. mean_rank（AUC 型统计量）：窗口内携带细胞平均 normalized rank
  mr <- tapply(df$rank, df$gene, mean)
  out$mean_rank <- unname(mr[out$sgrna])

  out <- out[order(out$sgrna), ]
  rownames(out) <- NULL
  out
}


#' 直接法 gene 级整合（Layer 2：mean_rank 经 RRA 聚合到基因）
#'
#' 把 \code{\link{direct_sgrna_rra}} 输出的 sgRNA 级表交给
#' \code{\link{run_rra_pipeline}}，以 \code{mean_rank} 为 score 列做
#' gene 级 RRA（join library 在此发生）。mean_rank 升序 rank 与 shape 法
#' linear 升序语义一致：p_high 显著 = 促病毒方向，p_low 显著 = 限制因子方向。
#'
#' @param sgrna_level data.frame（来自 \code{\link{direct_sgrna_rra}}）
#' @param library_df data.frame（来自 \code{\link{load_sgrna_library}}），
#'   含 \code{sgrna_id, gene_name}
#' @param rra_alpha RRA alpha 阈值，默认 0.1
#' @param rra_variant \code{"v1"} 或 \code{"v2"}
#' @param fdr_method p 值校正方法，默认 \code{"BH"}
#' @param sig_threshold 显著性阈值（火山图标注），默认 0.01
#'
#' @return list：
#'   \code{sgrna_level}（原样传入的 sgRNA 级表）、
#'   \code{gene_level}（gene 级 RRA 表，含 score/rank_high/rank_low/p_two/
#'   alpha_*/genename_* 可视化辅助列，可直接喂 \code{\link{plot_rra_volcano}}）
#' @export
#'
#' @examples
#' \dontrun{
#' sgl <- direct_sgrna_rra(seu, calls, phenotype_col = "virus_load")
#' lib <- load_sgrna_library("library.txt", sep = "\t")
#' dir_out <- run_direct_rra_pipeline(sgl, library_df = lib)
#' plot_rra_volcano(dir_out$gene_level, direction = "high")
#' }
run_direct_rra_pipeline <- function(sgrna_level,
                                    library_df,
                                    rra_alpha = 0.1,
                                    rra_variant = c("v1", "v2"),
                                    fdr_method = "BH",
                                    sig_threshold = 0.01) {
  stopifnot("'sgrna_level' 必须是含 sgrna, mean_rank 列的 data.frame" =
              is.data.frame(sgrna_level) &&
              all(c("sgrna", "mean_rank") %in% colnames(sgrna_level)))
  rra_variant <- match.arg(rra_variant)

  gene_level <- run_rra_pipeline(
    sgrna_level[, c("sgrna", "mean_rank"), drop = FALSE],
    score_cols = "mean_rank",
    library_df = library_df,
    rra_alpha = rra_alpha,
    rra_variant = rra_variant,
    fdr_method = fdr_method,
    sig_threshold = sig_threshold
  )[["mean_rank"]]

  list(sgrna_level = sgrna_level, gene_level = gene_level)
}
