#' 对每个 sgRNA 单独拟合 Poisson-Gaussian 混合模型，判定 positive cell
#'
#' 实现参考 Replogle et al. 2020 (Nature Biotechnology) 的 direct capture
#' perturb-seq guide calling 思路：对每个 guide，将 cell 分为 "background"
#' (低 UMI、ambient) 与 "expressing" (高 UMI、true transduction) 两群，
#' 求出阈值。一个 cell-sgRNA 对被判为 positive 当且仅当
#' \code{count >= min_umi} 且其后验概率属于 expressing 群。
#'
#' 实现细节：对每个 guide 的 cell UMI count 拟合两分量混合模型——Poisson 分量
#' 建模 background（ambient 低计数），Gaussian 分量建模 expressing（真实转导
#' 高计数），用 EM 算法估计参数并约束 Poisson 为低端分量，最后按后验概率
#' \code{P(expressing) > 0.5} 判定 positive（参考 Replogle et al. 2020 的
#' Poisson-Gaussian mixture guide calling）。当某 guide 拟合失败（如全零、
#' 单分量退化或低端分量不是 Poisson）时退回 \code{count >= min_umi} 的简单阈值。
#'
#' @param seu Seurat 对象（含 guide assay）或直接的 sgRNA × cell count 矩阵。
#'   传入 Seurat 时内部提取 guide assay 的 counts 层矩阵
#' @param guide_assay guide assay 名（仅在 \code{seu} 为 Seurat 对象时使用），
#'   默认 \code{"CRISPR Guide"}
#' @param method calling 方法：\code{"poisson_gaussian"}（默认，Poisson-Gaussian
#'   混合 EM）、\code{"per_cell"}（每个细胞内部保留 UMI 最高的前
#'   \code{max_guides_per_cell} 个 guide）、\code{"simple"}（简单阈值，
#'   \code{count >= min_umi} 即 positive）或 \code{"none"}（已 calling，直接
#'   \code{>0} 转 boolean）
#' @param min_umi positive cell 的最小 guide UMI 数，默认 3L（cellranger 推荐）
#' @param max_guides_per_cell 一个 cell 最多允许携带的 sgRNA 数。
#'   超过的 cell 视为 doublet 处理：\code{1L} 严格丢弃所有 multi-sgRNA cell；
#'   \code{Inf} 全保留；中间值保留前 N 个 UMI 最高的 sgRNA。
#' @param seed 随机种子，默认 123
#'
#' @return 同形状的 boolean 矩阵（cell × sgRNA），colnames = cell barcode，
#'   rownames = sgRNA id
#' @export
#'
#' @examples
#' \dontrun{
#' bool_mat <- call_sgrna(seu, guide_assay = "CRISPR Guide")
#' }
call_sgrna <- function(seu,
                       guide_assay = "CRISPR Guide",
                       method = c("poisson_gaussian", "per_cell", "simple", "none"),
                       min_umi = 3L,
                       max_guides_per_cell = 1L,
                       seed = 123) {
  method <- match.arg(method)

  ## 接受 Seurat 对象（内部提取 guide count）或直接的 count 矩阵
  if (inherits(seu, "Seurat")) {
    count_mat <- .get_guide_matrix(seu, guide_assay)
  } else {
    count_mat <- seu
  }
  stopifnot("count_mat 必须是 matrix 或 sparseMatrix" =
              is.matrix(count_mat) ||
              inherits(count_mat, "CsparseMatrix") ||
              inherits(count_mat, "TsparseMatrix"))
  stopifnot("min_umi 必须是非负整数" =
              (is.numeric(min_umi) && length(min_umi) == 1 &&
               min_umi >= 0 && min_umi == round(min_umi)))
  stopifnot("max_guides_per_cell 必须是 >=1 的整数或 Inf" =
              (is.numeric(max_guides_per_cell) && length(max_guides_per_cell) == 1 &&
               (max_guides_per_cell == Inf ||
                (max_guides_per_cell >= 1 &&
                 max_guides_per_cell == round(max_guides_per_cell)))))

  set.seed(seed)

  n_guides <- nrow(count_mat)
  n_cells <- ncol(count_mat)
  guide_ids <- rownames(count_mat)
  cell_ids <- colnames(count_mat)

  if (method == "none") {
    ## sparse 输入时 count_mat > 0 得到 lgCMatrix，storage.mode 无法处理 S4，
    ## 统一转 dense logical 与 poisson_gaussian 分支保持一致
    bool_mat <- as.matrix(count_mat > 0)
    storage.mode(bool_mat) <- "logical"
    dimnames(bool_mat) <- list(guide_ids, cell_ids)
  } else if (method == "simple") {
    ## 简单阈值：count >= min_umi 即 positive，不做任何混合拟合
    bool_mat <- as.matrix(count_mat >= min_umi)
    storage.mode(bool_mat) <- "logical"
    dimnames(bool_mat) <- list(guide_ids, cell_ids)
  } else if (method == "poisson_gaussian") {
    bool_mat <- matrix(FALSE, nrow = n_guides, ncol = n_cells,
                       dimnames = list(guide_ids, cell_ids))

    # 逐 guide 拟合
    for (g in seq_len(n_guides)) {
      counts <- as.numeric(count_mat[g, ])
      if (max(counts) == 0) next  # 全 0 跳过
      if (sum(counts > 0) < 5) next  # positive cell 太少跳过
      bool_mat[g, ] <- .call_one_guide(counts, min_umi = min_umi)
    }
  } else if (method == "per_cell") {
    bool_mat <- matrix(FALSE, nrow = n_guides, ncol = n_cells,
                       dimnames = list(guide_ids, cell_ids))
    # 每个细胞内部：保留 UMI 最高的前 max_guides_per_cell 个 guide（且 >= min_umi）
    for (j in seq_len(n_cells)) {
      cc <- as.numeric(count_mat[, j])
      if (max(cc) < min_umi) next
      ord <- order(cc, decreasing = TRUE)
      keep <- ord[seq_len(min(max_guides_per_cell, n_guides))]
      keep <- keep[cc[keep] >= min_umi]
      bool_mat[keep, j] <- TRUE
    }
  }

  # 处理 multi-guide cell
  if (!is.infinite(max_guides_per_cell)) {
    counts_per_cell <- colSums(bool_mat)
    multi_cells <- names(counts_per_cell)[counts_per_cell > max_guides_per_cell]
    if (length(multi_cells) > 0) {
      if (max_guides_per_cell == 1L) {
        # 严格丢弃：所有 multi-guide cell 全部置 FALSE
        bool_mat[, multi_cells] <- FALSE
      } else {
        # 保留前 N 个 UMI 最高的 guide
        for (cell in multi_cells) {
          cell_counts <- as.numeric(count_mat[, cell])
          keep_top <- order(cell_counts, decreasing = TRUE)[seq_len(max_guides_per_cell)]
          keep_mask <- rep(FALSE, n_guides)
          keep_mask[keep_top] <- TRUE
          bool_mat[, cell] <- bool_mat[, cell] & keep_mask
        }
      }
    }
  }

  bool_mat
}


#' 对单个 guide 拟合 Poisson-Gaussian 混合并判定 positive cell（内部函数）
#'
#' 在 log2 尺度上对每个 guide 的 cell UMI count 拟合两分量混合：Poisson 分量
#' （连续化，用 gamma 函数推广）建模 background（ambient 低计数），Gaussian
#' 分量建模 expressing（真实转导高计数）。用 EM 估计参数，约束 Poisson 为低端
#' 分量（\code{lambda < mu}），最后按后验概率 \code{P(expressing) > 0.5} 且
#' \code{count >= min_umi} 判定 positive。
#'
#' 注：guide UMI 高度偏斜（大量 ambient UMI=1-2，少量表达 UMI 达几十~几千），
#' 原始 count 尺度下 Gaussian 的 sigma 异常大导致分量不可分，故在 log2 尺度
#' 拟合（参考 Replogle et al. 2020）。
#'
#' @param counts 数值向量，每个 cell 对该 guide 的 UMI count（整数）
#' @param min_umi 最小 UMI 阈值（原始 count 尺度）
#' @return 逻辑向量，长度等于 counts
#' @keywords internal
.dpois_cont <- function(x, lambda) exp(x * log(lambda) - lambda - lgamma(x + 1))

.call_one_guide <- function(counts, min_umi) {
  x <- round(as.numeric(counts))              # 原始整数 UMI count
  n <- length(x)

  pos <- x[x > 0]
  if (length(pos) < 2) return(x >= min_umi)

  ## 只对非零 count 拟合（count=0 直接判 negative）；log2 尺度压缩长尾
  xn <- log2(pos)

  ## 初始化（分位数稳健估计两群；Poisson 为低端、Gaussian 为高端）
  lambda <- max(stats::median(xn), 0.01)        # Poisson 背景率（log2 尺度）
  mu <- stats::quantile(xn, 0.9, names = FALSE) # Gaussian 表达均值（log2 尺度）
  if (!is.finite(mu) || mu <= lambda) mu <- max(mean(xn), lambda + 0.5)
  sigma <- max(stats::sd(xn), 0.1)
  pi1 <- 0.5                                    # expressing 混合权重
  pi0 <- 1 - pi1

  ## EM 迭代（只在非零 cell 上进行）
  for (iter in seq_len(200L)) {
    dens0 <- pmax(.dpois_cont(xn, lambda), 1e-300)     # background
    dens1 <- pmax(stats::dnorm(xn, mu, sigma), 1e-300) # expressing

    gamma1 <- (pi1 * dens1) / (pi0 * dens0 + pi1 * dens1)
    gamma0 <- 1 - gamma1

    s0 <- sum(gamma0)
    s1 <- sum(gamma1)
    if (s0 < 1e-9 || s1 < 1e-9) break          # 单分量退化

    pi1_new <- s1 / length(xn)
    lambda_new <- max(sum(gamma0 * xn) / s0, 0.001)  # 下限防止 lambda=0 导致 log(0)
    mu_new <- sum(gamma1 * xn) / s1
    sigma_new <- sqrt(sum(gamma1 * (xn - mu_new)^2) / s1)
    sigma_new <- max(sigma_new, 0.1)

    if (abs(pi1_new - pi1) < 1e-8 &&
        abs(lambda_new - lambda) < 1e-8 &&
        abs(mu_new - mu) < 1e-8) {
      pi1 <- pi1_new; pi0 <- 1 - pi1_new
      lambda <- lambda_new; mu <- mu_new; sigma <- sigma_new
      break
    }
    pi1 <- pi1_new; pi0 <- 1 - pi1_new
    lambda <- lambda_new; mu <- mu_new; sigma <- sigma_new
  }

  ## 约束：Poisson 必须是低端分量（lambda < mu），否则退化退回简单阈值
  if (!is.finite(lambda) || !is.finite(mu) || lambda >= mu) {
    return(x >= min_umi)
  }

  ## 后验判定：log2 尺度算后验，原始 UMI 判 min_umi；count=0 判 negative
  dens0 <- pmax(.dpois_cont(xn, lambda), 1e-300)
  dens1 <- pmax(stats::dnorm(xn, mu, sigma), 1e-300)
  gamma1 <- (pi1 * dens1) / (pi0 * dens0 + pi1 * dens1)

  out <- logical(n)
  out[x > 0] <- gamma1 > 0.5 & (pos >= min_umi)
  out
}
