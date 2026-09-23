#' 计算表型分布的 breaks（基于分位数 + 步长）
#'
#' 根据 \code{phenotype_quantile} 截断两端 outlier，然后在
#' \code{[q01, q99]} 区间按 \code{break_step} 步长生成等距 breaks，
#' 保留到 0.01 精度。
#'
#' @param phenotype 命名数值向量，names = cell barcode
#' @param phenotype_quantile 长度 2 数值向量，默认 \code{c(0.01, 0.99)}
#' @param break_step 步长，默认 \code{0.01}
#'
#' @return 数值向量（breaks）
#' @export
#'
#' @examples
#' pheno <- setNames(rnorm(100), paste0("cell", 1:100))
#' compute_breaks(pheno)
compute_breaks <- function(phenotype,
                           phenotype_quantile = c(0.01, 0.99),
                           break_step = 0.01) {
  stopifnot("'phenotype' 必须是命名数值向量" =
              is.numeric(phenotype) && !is.null(names(phenotype)))
  stopifnot("'phenotype_quantile' 长度必须为 2 且在 [0,1] 范围" =
              length(phenotype_quantile) == 2 &&
              all(phenotype_quantile >= 0) && all(phenotype_quantile <= 1))
  stopifnot("'break_step' 必须是正数" =
              is.numeric(break_step) && length(break_step) == 1 && break_step > 0)

  q01 <- stats::quantile(phenotype, probs = phenotype_quantile[1])
  q99 <- stats::quantile(phenotype, probs = phenotype_quantile[2])
  breaks <- round(seq(q01, q99, by = break_step), 2)
  breaks
}


#' 计算全局归一化频率分布
#'
#' 在给定 breaks 上对表型值做 \code{hist()}，返回每个 bin 的归一化频率
#' （\code{freq / sum(freq)}）。
#'
#' @param phenotype 命名数值向量
#' @param breaks 数值向量（来自 \code{\link{compute_breaks}}）
#'
#' @return data.frame，列：\code{start, end, freq}（freq 已归一化）
#' @export
#'
#' @examples
#' pheno <- setNames(rnorm(100), paste0("cell", 1:100))
#' brks <- compute_breaks(pheno)
#' global_hist <- compute_global_hist(pheno, brks)
compute_global_hist <- function(phenotype, breaks) {
  stopifnot("'phenotype' 必须是命名数值向量" =
              is.numeric(phenotype) && !is.null(names(phenotype)))
  stopifnot("'breaks' 必须是单调递增数值向量" =
              is.numeric(breaks) && length(breaks) >= 2 &&
              all(diff(breaks) > 0))

  # 只保留落在 breaks 范围内的 phenotype，避免 hist 报 "some 'x' not counted"
  in_range <- phenotype >= min(breaks) & phenotype <= max(breaks)
  pheno_in <- phenotype[in_range]

  if (length(pheno_in) == 0) {
    return(data.frame(start = head(breaks, -1),
                      end = tail(breaks, -1),
                      freq = rep(0, length(breaks) - 1)))
  }

  h <- graphics::hist(pheno_in, breaks = breaks, plot = FALSE,
                      include.lowest = TRUE, right = TRUE)
  freq <- h$counts / sum(h$counts)
  data.frame(start = head(h$breaks, -1),
             end = tail(h$breaks, -1),
             freq = freq)
}


#' 对单个 sgRNA 的表型分布做正交多项式分解（固定 2 阶）
#'
#' 把该 sgRNA 携带的 cell 的表型分布与全局分布相减得到偏移，
#' 再用 \code{stats::poly()} 的 2 阶正交多项式基投影得到系数：
#' \code{linear}（整体偏移方向）与 \code{quadratic}（分布形状/宽度偏移）。
#'
#' @param per_sgrna_phenotype 命名数值向量，该 sgRNA 携带 cell 的表型
#' @param global_hist data.frame（来自 \code{\link{compute_global_hist}}）
#' @param weights 每个 cell 的权重（如 UMI 占比软去卷积权重），长度需与
#'   \code{per_sgrna_phenotype} 一致。默认 \code{NULL} 即等权（每 cell 计 1）
#'
#' @return 命名数值向量 \code{c(linear, quadratic)}
#' @export
#'
#' @examples
#' pheno <- setNames(rnorm(100), paste0("cell", 1:100))
#' brks <- compute_breaks(pheno)
#' global_hist <- compute_global_hist(pheno, brks)
#' sub_pheno <- pheno[1:30]
#' shape <- decompose_shape(sub_pheno, global_hist)
decompose_shape <- function(per_sgrna_phenotype, global_hist, weights = NULL) {
  stopifnot("'per_sgrna_phenotype' 必须是数值向量" =
              is.numeric(per_sgrna_phenotype))
  stopifnot("'global_hist' 必须是 compute_global_hist 输出的 data.frame" =
              is.data.frame(global_hist) &&
              all(c("start", "end", "freq") %in% colnames(global_hist)))
  if (!is.null(weights)) {
    stopifnot("'weights' 必须是数值向量且长度与 per_sgrna_phenotype 一致" =
                is.numeric(weights) &&
                length(weights) == length(per_sgrna_phenotype))
  }

  degree <- 2L
  # 去掉 NA（避免 cut/split 传播 NA level）
  valid_idx <- !is.na(per_sgrna_phenotype)
  per_sgrna_phenotype <- per_sgrna_phenotype[valid_idx]
  if (!is.null(weights)) weights <- weights[valid_idx]

  breaks <- c(global_hist$start[1], global_hist$end)
  # 只保留落在 breaks 范围内的 phenotype
  in_range <- per_sgrna_phenotype >= min(breaks) & per_sgrna_phenotype <= max(breaks)
  pheno_in <- per_sgrna_phenotype[in_range]
  w_in <- if (is.null(weights)) rep(1, length(per_sgrna_phenotype)) else weights
  w_in <- w_in[in_range]

  if (length(pheno_in) == 0) {
    shape <- rep(0, degree)
  } else {
    # 加权计数：cut 与 hist(include.lowest=TRUE, right=TRUE) 语义一致；
    # split(drop=FALSE) 保留所有 bin level，空 bin 权重和为 0
    bin <- cut(pheno_in, breaks = breaks, include.lowest = TRUE, right = TRUE)
    freq_num <- as.numeric(vapply(split(w_in, bin), sum, numeric(1)))
    s <- sum(freq_num)
    if (s == 0) {
      shape <- rep(0, degree)
    } else {
      freq_offset <- freq_num / s - global_hist$freq
      freq_offset[is.nan(freq_offset) | is.infinite(freq_offset)] <- 0
      poly_basis <- stats::poly((global_hist$start + global_hist$end) / 2,
                                degree = degree)
      # freq_offset 长度 = nrow(global_hist) = nrow(poly_basis)
      # 用 crossprod 实现 t(freq_offset) %*% poly_basis（1 x degree）
      shape <- as.numeric(crossprod(freq_offset, poly_basis))
    }
  }
  names(shape) <- c("linear", "quadratic")
  shape
}


#' 对所有 sgRNA 计算 shape 分解系数（linear + quadratic）
#'
#' shape 方案的主函数：先用 \code{phenotype_quantile} 或 \code{phenotype_range}
#' 截断表型区间，在截断后的数据上算 breaks 与全局 hist，然后对每个 sgRNA 在其
#' 携带 cell 上做 2 阶正交多项式分解，输出 sgRNA-level 系数表。
#'
#' 注意：本函数不依赖 sgRNA library，只用 sgRNA id 做分组。
#' library join 推迟到 \code{\link{run_rra_pipeline}}。
#'
#' @param seu Seurat 对象（含表型 metadata 列）或直接的命名数值表型向量。
#'   传入 Seurat 时内部从 \code{phenotype_col} 提取表型向量
#' @param sgrna_calls boolean matrix（sgRNA × cell），rownames = sgRNA id，
#'   colnames = cell barcode。来自 \code{\link{call_sgrna}}
#' @param phenotype_col metadata 中连续表型评分列名（仅在 \code{seu} 为 Seurat
#'   对象时使用），默认 \code{NULL}
#' @param phenotype_quantile 长度 2 数值向量，默认 \code{c(0.01, 0.99)}。
#'   仅在 \code{phenotype_range = NULL} 时生效：用下/上分位数截断表型 outlier
#' @param phenotype_range 长度 2 数值向量 \code{c(low, high)} 或 \code{NULL}（默认）。
#'   提供时忽略 \code{phenotype_quantile}：只保留表型落在 \code{[low, high]} 的
#'   cell，且 breaks 直接定义为 \code{seq(low, high, by = break_step)}。
#'   适用于有先验有效区间（如病毒载量 log2 scale 的 \code{c(9, 13)}）的场景
#' @param break_step 步长，默认 \code{0.01}
#' @param min_cells_per_sgrna 每 sgRNA 至少携带多少 cell 才保留在输出中。
#'   默认 20（与 bin 数解耦，bin 稀疏度不影响 rank，详见 vignette）
#' @param guide_assay CRISPR Guide assay 名，默认 \code{"CRISPR Guide"}。
#'   仅在 \code{weight_by_umi = TRUE} 且 \code{seu} 为 Seurat 对象时使用
#' @param weight_by_umi 是否用 UMI 占比对多 sgRNA cell 做软去卷积加权，默认
#'   \code{TRUE}。多 sgRNA cell 的每个 sgRNA 按 \code{UMI / 该 cell positive
#'   sgRNA UMI 总和} 软分配，单一 sgRNA cell 权重为 1
#'
#' @return data.frame，每行一个 sgRNA，列：
#'   \code{sgrna}（sgRNA id）、\code{linear}、\code{quadratic}（shape 系数）、
#'   \code{freq}（携带 cell 数）
#' @export
#'
#' @examples
#' \dontrun{
#' calls <- call_sgrna(seu)
#' sgrna_summary <- summarize_sgrna_shapes(seu, calls, phenotype_col = "virus_load")
#' }
summarize_sgrna_shapes <- function(seu,
                                    sgrna_calls,
                                    phenotype_col = NULL,
                                    phenotype_quantile = c(0.01, 0.99),
                                    phenotype_range = NULL,
                                    break_step = 0.01,
                                    min_cells_per_sgrna = 20L,
                                    guide_assay = "CRISPR Guide",
                                    weight_by_umi = TRUE) {
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

  # 1. 确定截断区间：固定区间优先，否则用分位数
  if (!is.null(phenotype_range)) {
    lo <- phenotype_range[1]
    hi <- phenotype_range[2]
  } else {
    lo <- stats::quantile(phenotype, probs = phenotype_quantile[1])
    hi <- stats::quantile(phenotype, probs = phenotype_quantile[2])
  }
  keep_cells <- names(phenotype)[phenotype >= lo & phenotype <= hi]
  phenotype_filt <- phenotype[keep_cells]

  # 2. breaks + 全局 hist
  breaks <- round(seq(lo, hi, by = break_step), 2)
  # 确保 max(breaks) 覆盖 hi（避免 seq 末点 < hi 时丢边界 cell）
  while (max(breaks) < hi) {
    breaks <- c(breaks, max(breaks) + break_step)
  }
  # 进一步把 phenotype_filt 限定到 breaks 范围（浮点安全）
  phenotype_filt <- phenotype_filt[
    phenotype_filt >= min(breaks) & phenotype_filt <= max(breaks)
  ]
  global_hist <- compute_global_hist(phenotype_filt, breaks)

  # 3. UMI 占比权重（软去卷积：多 sgRNA 细胞按 UMI 占比软分配）
  if (weight_by_umi && inherits(seu, "Seurat")) {
    guide_count <- .get_guide_matrix(seu, guide_assay)      # sgRNA × cell
    common_cells <- intersect(colnames(guide_count), colnames(sgrna_calls))
    common_guides <- intersect(rownames(guide_count), rownames(sgrna_calls))
    guide_count <- guide_count[common_guides, common_cells, drop = FALSE]
    calls_sp <- Matrix::Matrix(
      sgrna_calls[common_guides, common_cells, drop = FALSE], sparse = TRUE
    )
    pos_umi <- guide_count * calls_sp                     # 仅 positive 位置保留 UMI
    cell_total <- Matrix::colSums(pos_umi)                # 每细胞 positive UMI 总和
    cell_total[cell_total <= 0] <- 1                      # 避免除零
  } else {
    cell_total <- NULL
  }

  # 4. per-sgrna shape 分解
  guide_ids <- rownames(sgrna_calls)
  if (is.null(guide_ids)) {
    stop("'sgrna_calls' 必须有 rownames（sgRNA id）")
  }

  result_list <- lapply(guide_ids, function(gid) {
    pos_cells <- colnames(sgrna_calls)[sgrna_calls[gid, ]]
    pos_cells <- intersect(pos_cells, names(phenotype_filt))
    if (length(pos_cells) < 1) {
      return(data.frame(sgrna = gid, linear = NA_real_,
                        quadratic = NA_real_, freq = 0L,
                        stringsAsFactors = FALSE))
    }
    pheno_sub <- phenotype_filt[pos_cells]
    w <- if (is.null(cell_total)) NULL else
         as.numeric(guide_count[gid, pos_cells] / cell_total[pos_cells])
    shape <- decompose_shape(pheno_sub, global_hist, weights = w)
    data.frame(sgrna = gid,
               linear = unname(shape["linear"]),
               quadratic = unname(shape["quadratic"]),
               freq = length(pos_cells),
               stringsAsFactors = FALSE)
  })

  result <- do.call(rbind, result_list)
  result <- result[, c("sgrna", "linear", "quadratic", "freq")]

  # 5. min_cells_per_sgrna 过滤
  result <- result[result$freq >= min_cells_per_sgrna, ]
  rownames(result) <- NULL

  # 6. 附带 breaks / global_hist / 使用的区间，供下游（如 bin-offset 图）复用，
  #    保证可视化与 shape 分解用完全一致的 bin 定义
  attr(result, "breaks") <- breaks
  attr(result, "global_hist") <- global_hist
  attr(result, "phenotype_range") <- c(min(breaks), max(breaks))
  attr(result, "n_cells_in_range") <- length(phenotype_filt)
  result
}
