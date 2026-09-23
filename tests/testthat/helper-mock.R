# 测试辅助：mock 数据生成器（不导出）

# 生成 mock sgRNA library（mageck 三列风格）
mock_library <- function(n_guides = 20, n_genes = 5) {
  gene_names <- paste0("gene", seq_len(n_genes))
  data.frame(
    sgrna_id = paste0("g", seq_len(n_guides)),
    sequence = replicate(n_guides,
      paste0(sample(c("A", "C", "G", "T"), 20, replace = TRUE), collapse = "")),
    gene_name = sample(gene_names, n_guides, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# 生成 mock 表型向量（命名数值）
mock_phenotype <- function(n_cells = 200, shape = c("unimodal", "bimodal")) {
  shape <- match.arg(shape)
  vals <- if (shape == "unimodal") {
    rnorm(n_cells, mean = 10, sd = 1)
  } else {
    c(rnorm(n_cells / 2, mean = 9, sd = 0.5),
      rnorm(n_cells / 2, mean = 12, sd = 0.5))
  }
  setNames(vals, paste0("cell", seq_len(n_cells)))
}

# 生成 mock guide count 矩阵（cell × sgRNA）
mock_count_matrix <- function(n_cells = 200, n_guides = 20,
                              min_umi = 3) {
  cells <- paste0("cell", seq_len(n_cells))
  guides <- paste0("g", seq_len(n_guides))

  mat <- matrix(0L, nrow = n_guides, ncol = n_cells,
                dimnames = list(guides, cells))

  # 每个 cell 随机选 1 个 sgRNA 携带，count 来自 Poisson(8)
  for (i in seq_len(n_cells)) {
    g <- sample(n_guides, 1)
    mat[g, i] <- rpois(1, lambda = 8)
  }
  # 加少量背景 UMI
  bg <- matrix(rpois(n_guides * n_cells, lambda = 0.1),
               nrow = n_guides, ncol = n_cells)
  mat <- mat + bg
  dimnames(mat) <- list(guides, cells)
  mat
}

# 生成 mock sgrna_calls（boolean 矩阵）
mock_calls <- function(n_cells = 200, n_guides = 20) {
  cells <- paste0("cell", seq_len(n_cells))
  guides <- paste0("g", seq_len(n_guides))
  mat <- matrix(FALSE, nrow = n_guides, ncol = n_cells,
                dimnames = list(guides, cells))
  for (i in seq_len(n_cells)) {
    g <- sample(n_guides, 1)
    mat[g, i] <- TRUE
  }
  mat
}

# 生成 mock sgrna_summary（带 shape 系数，不包含 gene_name 列）
# 与真实 summarize_sgrna_shapes 输出列对齐：sgrna, linear, quadratic, freq
mock_sgrna_summary <- function(n_guides = 20, n_genes = 5) {
  lib <- mock_library(n_guides, n_genes)
  data.frame(
    sgrna     = lib$sgrna_id,
    linear    = rnorm(n_guides),
    quadratic = rnorm(n_guides),
    freq      = sample(20:100, n_guides, replace = TRUE),
    stringsAsFactors = FALSE
  )
}
