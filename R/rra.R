#' MAGeCK 风格的 RRA 算法（v1）
#'
#' 基于原版 alpha_rra（mageck_demo.R）：对每个 gene，把其 sgRNA 的 rank
#' 归一化为 ratio_rank = rank/(N_total+1)，再用 Beta 分布检验
#' 是否聚集在分布的低端（p_low）或高端（p_high）。
#'
#' p_high：检验 sgRNA 是否聚集在 ratio_rank 小端（对应 shape 系数低）。
#' p_low：检验 sgRNA 是否聚集在 ratio_rank 大端（对应 shape 系数高）。
#'
#' 命名上 \code{p_high/p_low} 沿用 MAGeCK 习惯指代 rank 的两端，
#' 具体语义取决于调用者传入的 rank 方向（升序或降序）。
#'
#' @param df data.frame，必须含列：\code{gene, sgrna, rank}
#' @param N_total 总 sgRNA 数（用于 ratio_rank 归一化），默认 \code{max(df$rank)}
#' @param alpha 阈值参数，默认 0.05
#' @param fdr_method p 值校正方法，默认 \code{"BH"}
#'
#' @return data.frame，列：
#'   \code{gene, sgRNA_count, sgRNA_high_in_alpha, sgRNA_low_in_alpha,
#'         p_high, p_low, fdr_high, fdr_low}
#' @export
#'
#' @examples
#' df <- data.frame(
#'   gene = c("g1", "g1", "g2", "g2"),
#'   sgrna = c("s1", "s2", "s3", "s4"),
#'   rank = c(1, 2, 50, 100)
#' )
#' alpha_rra(df, alpha = 0.1)
alpha_rra <- function(df, N_total = NULL, alpha = 0.05, fdr_method = "BH") {
  stopifnot("'df' 必须包含列 gene, sgrna, rank" =
              all(c("gene", "sgrna", "rank") %in% names(df)))

  if (is.null(N_total)) {
    N_total <- max(df$rank, na.rm = TRUE)
  }

  df <- df[!is.na(df$rank), ]
  df <- df[order(df$gene, df$rank), ]

  df$ratio_rank <- df$rank / (N_total + 1)

  df_split <- split(df, df$gene)

  result_list <- lapply(df_split, function(data_tmp) {
    n <- nrow(data_tmp)

    x_high <- sort(data_tmp$ratio_rank)
    x_low  <- sort(1 - data_tmp$ratio_rank)

    k <- seq_len(n)

    p_high_all <- stats::pbeta(x_high, shape1 = k, shape2 = n + 1 - k)
    p_low_all  <- stats::pbeta(x_low,  shape1 = k, shape2 = n + 1 - k)

    j_high <- sum(x_high <= alpha)
    j_low  <- sum(x_low  <= alpha)

    p_high <- if (j_high > 0) min(p_high_all[seq_len(j_high)]) else 1
    p_low  <- if (j_low  > 0) min(p_low_all[seq_len(j_low)])   else 1

    data.frame(
      gene = data_tmp$gene[1],
      sgRNA_count = n,
      sgRNA_high_in_alpha = j_high,
      sgRNA_low_in_alpha = j_low,
      p_high = p_high,
      p_low = p_low,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, result_list)
  rownames(result_df) <- NULL

  result_df$fdr_high <- stats::p.adjust(result_df$p_high, method = fdr_method)
  result_df$fdr_low  <- stats::p.adjust(result_df$p_low,  method = fdr_method)

  result_df <- result_df[order(result_df$fdr_low, result_df$fdr_high), ]
  rownames(result_df) <- NULL

  result_df
}


#' MAGeCK 风格的 RRA 算法（v2，带 n_eff 收缩）
#'
#' 与 \code{\link{alpha_rra}} 的区别：用 \code{n_eff = max(j_high, min(2, n))}
#' 替代原始 n 作为 Beta 分布的 shape2 参数，对 sgRNA 数较少的 gene
#' 做 shrinkage，缓解小样本下 p 值过乐观的问题。
#'
#' @inheritParams alpha_rra
#'
#' @return 同 \code{\link{alpha_rra}}，附加列 \code{n_eff}
#' @export
#'
#' @examples
#' df <- data.frame(
#'   gene = c("g1", "g1", "g2", "g2"),
#'   sgrna = c("s1", "s2", "s3", "s4"),
#'   rank = c(1, 2, 50, 100)
#' )
#' alpha_rra_v2(df, alpha = 0.1)
alpha_rra_v2 <- function(df, N_total = NULL, alpha = 0.05, fdr_method = "BH") {
  stopifnot("'df' 必须包含列 gene, sgrna, rank" =
              all(c("gene", "sgrna", "rank") %in% names(df)))

  if (is.null(N_total)) {
    N_total <- max(df$rank, na.rm = TRUE)
  }

  df <- df[!is.na(df$rank), ]
  df <- df[order(df$gene, df$rank), ]

  df$ratio_rank <- df$rank / (N_total + 1)

  df_split <- split(df, df$gene)

  result_list <- lapply(df_split, function(data_tmp) {
    n <- nrow(data_tmp)

    x_high <- sort(data_tmp$ratio_rank)
    x_low  <- sort(1 - data_tmp$ratio_rank)

    j_high <- sum(x_high <= alpha)
    j_low  <- sum(x_low  <= alpha)

    # n_eff 至少为 min(2, n)，但不超过实际 sgRNA 数 n
    n_eff <- max(j_high, min(2, n))
    n_eff <- min(n_eff, n)

    p_high <- if (j_high > 0) {
      k <- seq_len(j_high)
      p_all <- stats::pbeta(x_high[seq_len(j_high)],
                            shape1 = k, shape2 = n_eff + 1 - k)
      min(p_all)
    } else 1

    p_low <- if (j_low > 0) {
      k <- seq_len(j_low)
      p_all <- stats::pbeta(x_low[seq_len(j_low)],
                            shape1 = k, shape2 = n_eff + 1 - k)
      min(p_all)
    } else 1

    data.frame(
      gene = data_tmp$gene[1],
      sgRNA_count = n,
      sgRNA_high_in_alpha = j_high,
      sgRNA_low_in_alpha = j_low,
      n_eff = n_eff,
      p_high = p_high,
      p_low = p_low,
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, result_list)
  result_df$fdr_high <- stats::p.adjust(result_df$p_high, method = fdr_method)
  result_df$fdr_low  <- stats::p.adjust(result_df$p_low,  method = fdr_method)

  result_df
}


#' 加载 sgRNA library（mageck library.txt 风格）
#'
#' 按列位置读取三列：sgRNA id、sequence、gene_name。不依赖列名。
#' 仅做最基本的校验（三列齐全、无 NA id、无 NA gene_name）。
#'
#' @param path 文件路径
#' @param col_idx 长度 3 整数向量，默认 \code{c(1, 2, 3)}，
#'   分别对应 sgRNA id、sequence、gene_name 的列位置
#' @param sep 分隔符，默认 \code{","}
#' @param header 是否有表头，默认 \code{TRUE}
#'
#' @return data.frame，列名固定为 \code{sgrna_id, sequence, gene_name}
#' @export
#'
#' @examples
#' \dontrun{
#' lib <- load_sgrna_library("library.txt", sep = "\t")
#' }
load_sgrna_library <- function(path,
                               col_idx = c(1, 2, 3),
                               sep = ",",
                               header = TRUE) {
  stopifnot("'col_idx' 长度必须为 3" = length(col_idx) == 3)
  stopifnot("'col_idx' 必须是正整数" =
              all(col_idx == round(col_idx)) && all(col_idx >= 1))

  df <- utils::read.table(path, sep = sep, header = header,
                          stringsAsFactors = FALSE)
  if (ncol(df) < 3) {
    stop(sprintf("文件仅有 %d 列，library.txt 至少需要 3 列", ncol(df)))
  }

  out <- df[, col_idx]
  colnames(out) <- c("sgrna_id", "sequence", "gene_name")

  if (any(is.na(out$sgrna_id) | out$sgrna_id == "")) {
    stop("library 的 sgRNA id 列存在 NA 或空字符串")
  }
  if (any(is.na(out$gene_name) | out$gene_name == "")) {
    stop("library 的 gene_name 列存在 NA 或空字符串")
  }

  out
}


#' 对 shape 系数表跑 RRA 流程（join library 在此发生）
#'
#' 这是 stage 3 的主函数：把 shape 分解输出的 sgRNA-level 系数表
#' 与 library join，得到 gene 分组，然后对每个 score 列单独跑
#' \code{\link{alpha_rra}} 或 \code{\link{alpha_rra_v2}}。
#'
#' 输出包含 score/rank_high/rank_low/p_two/alpha_high/alpha_low/
#' genename_high/genename_low 等可视化辅助列，可直接喂给
#' \code{\link{plot_rra_volcano}}。
#'
#' @param sgrna_summary data.frame（来自 \code{\link{summarize_sgrna_shapes}}），
#'   必须含列 \code{sgrna}（sgRNA id）
#' @param score_cols 字符向量，要跑 RRA 的 shape 系数列名，
#'   默认 \code{c("linear", "quadratic")}
#' @param library_df data.frame（来自 \code{\link{load_sgrna_library}}），
#'   含 \code{sgrna_id, sequence, gene_name}
#' @param rra_alpha RRA alpha 阈值，默认 0.1
#' @param rra_variant \code{"v1"} 或 \code{"v2"}
#' @param fdr_method p 值校正方法，默认 \code{"BH"}
#' @param sig_threshold 显著性阈值，用于 alpha/genename 标注，默认 0.01
#'
#' @return list，长度 = length(score_cols)，每个元素是带 score/rank/
#'   p_two/alpha/genename 的 data.frame
#' @export
#'
#' @examples
#' \dontrun{
#' sgrna_summary <- summarize_sgrna_shapes(pheno, calls)
#' lib <- load_sgrna_library("library.txt", sep = "\t")
#' rra_out <- run_rra_pipeline(sgrna_summary,
#'                              score_cols = c("linear", "quadratic"),
#'                              library_df = lib)
#' }
run_rra_pipeline <- function(sgrna_summary,
                             score_cols = c("linear", "quadratic"),
                             library_df,
                             rra_alpha = 0.1,
                             rra_variant = c("v1", "v2"),
                             fdr_method = "BH",
                             sig_threshold = 0.01) {
  stopifnot("'sgrna_summary' 必须是含 sgrna 列的 data.frame" =
              is.data.frame(sgrna_summary) && "sgrna" %in% colnames(sgrna_summary))
  stopifnot("'library_df' 必须是含 sgrna_id, gene_name 列的 data.frame" =
              is.data.frame(library_df) &&
              all(c("sgrna_id", "gene_name") %in% colnames(library_df)))
  stopifnot("'score_cols' 必须存在于 sgrna_summary" =
              all(score_cols %in% colnames(sgrna_summary)))

  rra_variant <- match.arg(rra_variant)

  # join sgrna_summary 与 library（内连接）
  # 只取 library 的 sgrna_id 和 gene_name 列，并重命名为 sgrna/gene_name_lib
  # 避免 sgrna_summary 已含 gene_name 列时产生列名冲突
  lib_subset <- library_df[, c("sgrna_id", "gene_name")]
  colnames(lib_subset) <- c("sgrna", "gene_name")
  merged <- merge(sgrna_summary, lib_subset,
                  by = "sgrna", all.x = FALSE)
  if (nrow(merged) == 0) {
    stop("library 与 sgrna_summary 在 sgrna id 上没有交集，请检查 ID 一致性")
  }

  # gene 分组（用于算 score）
  gene_groups <- split(merged, merged$gene_name)

  results <- lapply(score_cols, function(col) {
    rra_df <- data.frame(
      gene = merged$gene_name,
      sgrna = merged$sgrna,
      rank = rank(merged[[col]], ties.method = "average")
    )

    rra_out <- if (rra_variant == "v1") {
      alpha_rra(rra_df, N_total = nrow(rra_df), alpha = rra_alpha,
                fdr_method = fdr_method)
    } else {
      alpha_rra_v2(rra_df, N_total = nrow(rra_df), alpha = rra_alpha,
                   fdr_method = fdr_method)
    }

    # score = mean(原始 shape 系数 per gene)
    rra_out$score <- vapply(rra_out$gene, function(g) {
      sub <- gene_groups[[g]]
      if (is.null(sub) || nrow(sub) == 0) NA_real_ else mean(sub[[col]])
    }, numeric(1))

    # rank_high, rank_low
    rra_out <- rra_out[order(rra_out$p_high), ]
    rra_out$rank_high <- seq_len(nrow(rra_out))
    rra_out <- rra_out[order(rra_out$p_low), ]
    rra_out$rank_low <- seq_len(nrow(rra_out))

    # p_two（双侧）
    rra_out$p_two <- ifelse(rra_out$p_high >= rra_out$p_low,
                            rra_out$p_low * 2,
                            rra_out$p_high * 2)

    # alpha/genename 标注（火山图用）
    rra_out$alpha_high <- ifelse(rra_out$p_high <= sig_threshold, 1, 0.3)
    rra_out$genename_high <- ifelse(rra_out$p_high <= sig_threshold,
                                    rra_out$gene, "")
    rra_out$alpha_low <- ifelse(rra_out$p_low <= sig_threshold, 1, 0.3)
    rra_out$genename_low <- ifelse(rra_out$p_low <= sig_threshold,
                                   rra_out$gene, "")
    rownames(rra_out) <- NULL
    rra_out
  })

  names(results) <- score_cols
  results
}
