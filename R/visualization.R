#' RRA 火山图（high 或 low 侧）
#'
#' \code{direction = "high"} 画 p_high 侧（检验 sgRNA 聚集在 shape 系数低端），
#' \code{direction = "low"} 画 p_low 侧（聚集在高端）。
#' x 轴为 score（gene 的 shape 系数均值），y 轴为 -log10(FDR)，
#' 显著 gene 的 label 自动通过 \code{ggrepel} 避免重叠。
#'
#' @param rra_result data.frame（来自 \code{\link{run_rra_pipeline}} 的某个元素）
#' @param direction \code{"high"} 或 \code{"low"}
#' @param sig_threshold 显著性阈值（与 \code{run_rra_pipeline} 的 \code{sig_threshold} 一致），
#'   默认 0.01；用于决定 alpha 透明度与 label 显示
#' @param ylim_max y 轴上限，默认 1.2
#' @param midpoint \code{scale_color_gradient2} 的 midpoint，默认 10
#' @param y_value y 轴取值：\code{"fdr"}（默认，BH 校正 FDR）或 \code{"p"}
#'   （未校正 p 值）
#'
#' @return ggplot 对象
#' @export
#'
#' @examples
#' \dontrun{
#' rra_out <- run_rra_pipeline(sgrna_summary, library_df = lib)
#' plot_rra_volcano(rra_out[["linear"]], direction = "high")
#' }
plot_rra_volcano <- function(rra_result,
                             direction = c("high", "low"),
                             sig_threshold = 0.01,
                             ylim_max = 4,
                             midpoint = 10,
                             y_value = c("fdr", "p")) {
  direction <- match.arg(direction)
  y_value <- match.arg(y_value)

  y_col <- if (y_value == "fdr") {
    if (direction == "high") "fdr_high" else "fdr_low"
  } else {
    if (direction == "high") "p_high" else "p_low"
  }
  ylab <- paste0("-log10(", if (y_value == "fdr") "FDR" else "p", "_",
                 direction, ")")
  rank_col <- if (direction == "high") "rank_high" else "rank_low"
  alpha_col <- if (direction == "high") "alpha_high" else "alpha_low"
  genename_col <- if (direction == "high") "genename_high" else "genename_low"

  p <- ggplot2::ggplot(rra_result,
                       ggplot2::aes(.data[["score"]],
                                    -log10(.data[[y_col]]),
                                    colour = log2(.data[[rank_col]]))) +
    ggplot2::geom_point(alpha = rra_result[[alpha_col]]) +
    ggplot2::scale_color_gradient2(low = "#4575B4", mid = "yellow",
                                   high = "#D73027", midpoint = midpoint) +
    ggrepel::geom_text_repel(
      data = rra_result,
      ggplot2::aes(.data[["score"]], -log10(.data[[y_col]]),
                   label = .data[[genename_col]]),
      size = 3, box.padding = grid::unit(0.5, "lines"),
      point.padding = grid::unit(0.8, "lines"),
      segment.color = "black", show.legend = FALSE) +
    ggplot2::coord_cartesian(ylim = c(0, ylim_max)) +
    ggplot2::theme_bw()

  if (direction == "high") {
    p <- p + ggplot2::labs(x = "score", y = ylab,
                           title = "high-end enrichment")
  } else {
    p <- p + ggplot2::scale_color_gradient2(high = "#4575B4",
                                            mid = "yellow",
                                            low = "#D73027",
                                            midpoint = midpoint) +
      ggplot2::labs(x = "score", y = ylab,
                    title = "low-end enrichment")
  }
  p
}


#' 同时画 high/low 两侧的火山图
#'
#' 调用 \code{\link{plot_rra_volcano}} 两次（high 和 low），
#' 用 \code{gridExtra::grid.arrange} 拼成左右两图。
#'
#' @inheritParams plot_rra_volcano
#' @return gtable 对象（grid.arrange 输出）
#' @export
#'
#' @examples
#' \dontrun{
#' plot_rra_volcano_both(rra_out[["linear"]])
#' }
plot_rra_volcano_both <- function(rra_result, sig_threshold = 0.01,
                                  ylim_max = 1.2, midpoint = 10,
                                  y_value = c("fdr", "p")) {
  p1 <- plot_rra_volcano(rra_result, direction = "high",
                         sig_threshold = sig_threshold,
                         ylim_max = ylim_max, midpoint = midpoint,
                         y_value = y_value)
  p2 <- plot_rra_volcano(rra_result, direction = "low",
                         sig_threshold = sig_threshold,
                         ylim_max = ylim_max, midpoint = midpoint,
                         y_value = y_value)
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}


#' shape 系数分布 + 标记目标 sgRNA
#'
#' 画所有 sgRNA 的某个 shape 系数（如 linear）的直方图，
#' 并用虚线标出指定 sgRNA 的位置。
#'
#' @param sgrna_summary data.frame（来自 \code{\link{summarize_sgrna_shapes}}）
#' @param score_col 系数列名，默认 \code{"linear"}
#' @param target_sgrna 单个 sgRNA id（必须是 \code{sgrna_summary$sgrna} 中的值）
#' @param binwidth 直方图 bin 宽，默认 0.002
#'
#' @return ggplot 对象
#' @export
#'
#' @examples
#' \dontrun{
#' plot_shape_score_distribution(sgrna_summary,
#'                               target_sgrna = "g1_1")
#' }
plot_shape_score_distribution <- function(sgrna_summary,
                                          score_col = "linear",
                                          target_sgrna = NULL,
                                          binwidth = 0.002) {
  stopifnot("'score_col' 必须存在于 sgrna_summary" =
              score_col %in% colnames(sgrna_summary))

  target_val <- if (!is.null(target_sgrna)) {
    stopifnot("'target_sgrna' 必须在 sgrna_summary$sgrna 中" =
              target_sgrna %in% sgrna_summary$sgrna)
    sgrna_summary[[score_col]][sgrna_summary$sgrna == target_sgrna]
  } else NA_real_

  p <- ggplot2::ggplot(sgrna_summary, ggplot2::aes(.data[[score_col]])) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            binwidth = binwidth) +
    ggplot2::theme_minimal()

  if (!is.na(target_val)) {
    p <- p + ggplot2::geom_vline(xintercept = target_val,
                                 color = "blue",
                                 linetype = "dotted",
                                 linewidth = 0.8)
  }
  p
}


#' 单个 sgRNA 携带 cell 的表型分布直方图
#'
#' 对指定 sgRNA 的携带 cell 子集画表型直方图，并叠加 kernel density 曲线。
#' 多 sgRNA 可同时画，自动 facet。
#'
#' 注：原脚本（03.rra_for_rsv_v7.R）叠加双高斯密度曲线（依赖 Mclust），
#' 因 Mclust 已移除，此处改为单条 kernel density 曲线。
#'
#' @param phenotype 命名数值向量
#' @param sgrna_calls boolean matrix（sgRNA × cell）
#' @param target_sgrnas 字符向量，要画的 sgRNA id 列表
#' @param binwidth 直方图 bin 宽，默认 0.4
#'
#' @return ggplot 对象（facet by sgrna）
#' @export
#'
#' @examples
#' \dontrun{
#' plot_sgrna_distribution(phenotype, calls,
#'                         target_sgrnas = c("g1_1", "g1_2"))
#' }
plot_sgrna_distribution <- function(phenotype,
                                    sgrna_calls,
                                    target_sgrnas,
                                    binwidth = 0.4) {
  stopifnot("'phenotype' 必须是命名数值向量" =
              is.numeric(phenotype) && !is.null(names(phenotype)))
  stopifnot("'sgrna_calls' 必须是 logical matrix" =
              is.logical(sgrna_calls) && is.matrix(sgrna_calls))
  stopifnot("'target_sgrnas' 必须非空且都在 sgrna_calls 的 rownames 中" =
              length(target_sgrnas) > 0 &&
              all(target_sgrnas %in% rownames(sgrna_calls)))

  # 收集每个 sgRNA 的携带 cell 的表型值
  df_list <- lapply(target_sgrnas, function(gid) {
    pos_cells <- colnames(sgrna_calls)[sgrna_calls[gid, ]]
    pos_cells <- intersect(pos_cells, names(phenotype))
    if (length(pos_cells) == 0) return(NULL)
    data.frame(phenotype = phenotype[pos_cells], sgrna = gid)
  })
  df <- do.call(rbind, df_list)
  if (is.null(df) || nrow(df) == 0) {
    stop("目标 sgRNA 在指定 cell 集合中无 positive call")
  }

  ggplot2::ggplot(df, ggplot2::aes(.data[["phenotype"]])) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            binwidth = binwidth) +
    ggplot2::geom_density(color = "red3", linewidth = 1) +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~ sgrna)
}


#' per-sgrna 各 bin 的频率偏移条形图
#'
#' 对每个目标 sgRNA，画每个 bin 上
#' \code{(sgRNA 频率 - 全局频率)} 的条形图，可看出该 sgRNA
#' 相对全局分布的偏离方向。
#'
#' @param global_hist data.frame（来自 \code{\link{compute_global_hist}}）
#' @param phenotype 命名数值向量
#' @param sgrna_calls boolean matrix（sgRNA × cell）
#' @param target_sgrnas 字符向量
#' @param normalize_by_sqrt 是否除以 \code{sqrt(global_freq)}（与原脚本 L352 一致），
#'   默认 FALSE
#'
#' @return ggplot 对象（facet by sgrna）
#' @export
#'
#' @examples
#' \dontrun{
#' brks <- compute_breaks(phenotype)
#' gh <- compute_global_hist(phenotype, brks)
#' plot_sgrna_bin_offset(gh, phenotype, calls,
#'                       target_sgrnas = c("g1_1"))
#' }
plot_sgrna_bin_offset <- function(global_hist,
                                  phenotype,
                                  sgrna_calls,
                                  target_sgrnas,
                                  normalize_by_sqrt = FALSE) {
  stopifnot("'global_hist' 含 start/end/freq 列" =
              all(c("start", "end", "freq") %in% colnames(global_hist)))
  stopifnot("'target_sgrnas' 必须在 sgrna_calls 的 rownames 中" =
              all(target_sgrnas %in% rownames(sgrna_calls)))

  breaks <- c(global_hist$start[1], global_hist$end)

  df_list <- lapply(target_sgrnas, function(gid) {
    pos_cells <- colnames(sgrna_calls)[sgrna_calls[gid, ]]
    pos_cells <- intersect(pos_cells, names(phenotype))
    if (length(pos_cells) == 0) return(NULL)
    ## 仅统计落在 breaks 范围内的 cell（与 shape 分解阶段的 q01/q99 截断一致），
    ## 否则 hist 对越界值报 "有些 x 值没有被计数"
    x <- phenotype[pos_cells]
    x <- x[x >= min(breaks) & x <= max(breaks)]
    if (length(x) == 0) return(NULL)
    h <- graphics::hist(x, breaks = breaks, plot = FALSE)
    s <- sum(h$counts)
    if (s == 0) return(NULL)
    freq_offset <- h$counts / s - global_hist$freq
    if (normalize_by_sqrt) {
      freq_offset <- freq_offset / sqrt(global_hist$freq)
    }
    data.frame(start = head(h$breaks, -1),
               freq = freq_offset,
               sgrna = gid)
  })
  df <- do.call(rbind, df_list)
  if (is.null(df) || nrow(df) == 0) {
    stop("目标 sgRNA 无足够数据")
  }

  ggplot2::ggplot(df, ggplot2::aes(.data[["start"]], .data[["freq"]])) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_grid(~ sgrna) +
    ggplot2::theme_bw()
}
