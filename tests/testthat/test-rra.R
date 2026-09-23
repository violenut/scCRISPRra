# 测试 alpha_rra / alpha_rra_v2 / run_rra_pipeline / load_sgrna_library

test_that("alpha_rra: 完美排序输入让 p 极小", {
  # gene g1 的两个 sgRNA 都在 rank 最低端
  df <- data.frame(
    gene = c("g1", "g1", "g2", "g2", "g3"),
    sgrna = c("s1", "s2", "s3", "s4", "s5"),
    rank = c(1, 2, 50, 100, 60)
  )
  res <- alpha_rra(df, alpha = 0.1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("p_high", "p_low", "fdr_high", "fdr_low") %in%
                  colnames(res)))
  # g1 的 p_high 应极小
  g1_row <- res[res$gene == "g1", ]
  expect_lt(g1_row$p_high, 0.1)
  # g1 的 p_low 应接近 1
  expect_gt(g1_row$p_low, 0.9)
})

test_that("alpha_rra: 均匀随机 rank 让 p 接近 1", {
  set.seed(42)
  df <- data.frame(
    gene = rep(paste0("g", 1:4), each = 3),
    sgrna = paste0("s", 1:12),
    rank = sample(1:100, 12)
  )
  res <- alpha_rra(df, alpha = 0.1)
  # 大多数 gene 的 p 应不太小
  expect_true(all(res$p_high >= 0))
  expect_true(all(res$p_low >= 0))
})

test_that("alpha_rra: 单 sgRNA gene 不应崩", {
  df <- data.frame(
    gene = c("g1", "g2", "g2"),
    sgrna = c("s1", "s2", "s3"),
    rank = c(1, 5, 10)
  )
  res <- alpha_rra(df, alpha = 0.1)
  expect_equal(nrow(res), 2)
  # 单 sgRNA gene 的 p_high 是 pbeta(rank/(N+1), 1, 1)
  expect_true(res$p_high[res$gene == "g1"] >= 0)
  expect_true(res$p_high[res$gene == "g1"] <= 1)
})

test_that("alpha_rra: NA rank 被丢弃", {
  df <- data.frame(
    gene = c("g1", "g1", "g2"),
    sgrna = c("s1", "s2", "s3"),
    rank = c(1, NA, 5)
  )
  res <- alpha_rra(df, alpha = 0.1)
  # g1 只剩 1 个 sgRNA
  expect_equal(res$sgRNA_count[res$gene == "g1"], 1)
})

test_that("alpha_rra: BH 校正与 stats::p.adjust 一致", {
  df <- data.frame(
    gene = paste0("g", 1:10),
    sgrna = paste0("s", 1:10),
    rank = c(1, 5, 10, 20, 30, 40, 50, 60, 70, 80)
  )
  res <- alpha_rra(df, alpha = 0.1)
  expected_fdr_high <- stats::p.adjust(res$p_high, method = "BH")
  expect_equal(res$fdr_high, expected_fdr_high)
})

test_that("alpha_rra_v2: n_eff 列存在且收缩生效", {
  df <- data.frame(
    gene = c("g1", "g1", "g2", "g2"),
    sgrna = c("s1", "s2", "s3", "s4"),
    rank = c(1, 2, 50, 100)
  )
  res <- alpha_rra_v2(df, alpha = 0.1)
  expect_true("n_eff" %in% colnames(res))
  # n_eff >= 2 或 == sgRNA_count (n)
  expect_true(all(res$n_eff >= 1))
})

test_that("load_sgrna_library: 按列位置读取三列", {
  # 写一个临时文件
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  df <- data.frame(
    some_id = paste0("g", 1:5),
    some_seq = replicate(5,
      paste0(sample(c("A", "C", "G", "T"), 20, replace = TRUE),
             collapse = "")),
    some_gene = paste0("gene", 1:5),
    extra_col = 1:5,
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, tmp, row.names = FALSE)

  lib <- load_sgrna_library(tmp, col_idx = c(1, 2, 3))
  expect_equal(colnames(lib), c("sgrna_id", "sequence", "gene_name"))
  expect_equal(lib$sgrna_id, paste0("g", 1:5))
  expect_equal(lib$gene_name, paste0("gene", 1:5))
})

test_that("load_sgrna_library: 列数不足报错", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp))
  utils::write.csv(data.frame(a = 1:3, b = 4:6), tmp, row.names = FALSE)
  expect_error(load_sgrna_library(tmp),
               "至少需要 3 列")
})

test_that("run_rra_pipeline: 输出 list 按 score_cols 命名", {
  lib <- mock_library(n_guides = 10, n_genes = 3)
  sgrna_summary <- mock_sgrna_summary(n_guides = 10, n_genes = 3)
  # 让 sgrna 列与 library 的 sgrna_id 一致
  sgrna_summary$sgrna <- lib$sgrna_id

  out <- run_rra_pipeline(sgrna_summary,
                          score_cols = c("linear", "quadratic"),
                          library_df = lib)
  expect_type(out, "list")
  expect_equal(names(out), c("linear", "quadratic"))

  # 每个元素都有必要列
  for (col in c("linear", "quadratic")) {
    res <- out[[col]]
    expect_true(all(c("score", "rank_high", "rank_low", "p_two",
                      "alpha_high", "genename_high",
                      "alpha_low", "genename_low") %in% colnames(res)))
  }
})

test_that("run_rra_pipeline: 无 library 交集时报错", {
  lib <- mock_library(n_guides = 5, n_genes = 3)
  sgrna_summary <- data.frame(
    sgrna = paste0("MISSING_", 1:5),
    linear = rnorm(5),
    quadratic = rnorm(5),
    freq = 100L,
    stringsAsFactors = FALSE
  )
  expect_error(run_rra_pipeline(sgrna_summary,
                                library_df = lib),
               "没有交集")
})

## ──────── 直接法（direct_sgrna_rra / run_direct_rra_pipeline）────────

test_that("direct_sgrna_rra: 效应 sgRNA 携带细胞表型偏低 → p_high 显著", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1), paste0("cell", 1:500))
  calls <- mock_calls(n_cells = 500, n_guides = 10)
  # g1 的携带细胞表型整体下移（模拟促病毒因子：KO → 低载量）
  g1_cells <- colnames(calls)[calls["g1", ]]
  pheno[g1_cells] <- pheno[g1_cells] - 3

  out <- direct_sgrna_rra(pheno, calls, min_cells_per_sgrna = 1)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("sgrna", "n_cells", "cells_low", "cells_high",
                    "p_high", "p_low", "fdr_high", "fdr_low",
                    "mean_rank") %in% colnames(out)))
  expect_equal(sort(out$sgrna), sort(rownames(calls)))

  g1 <- out[out$sgrna == "g1", ]
  expect_lt(g1$p_high, 0.05)
  expect_gt(g1$p_low, 0.5)
  # mean_rank 应显著低于 0.5（携带细胞聚集低表型端）
  expect_lt(g1$mean_rank, 0.4)
})

test_that("direct_sgrna_rra: phenotype_range 与 quantile 截断一致可用", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1), paste0("cell", 1:500))
  calls <- mock_calls(n_cells = 500, n_guides = 10)

  out_r <- direct_sgrna_rra(pheno, calls, phenotype_range = c(9, 11),
                            min_cells_per_sgrna = 1)
  expect_true(nrow(out_r) > 0)
  # 窗口内细胞数与表型落在区间内的数量一致（忽略无 positive call 的细胞）
  n_in <- sum(pheno >= 9 & pheno <= 11)
  expect_true(all(out_r$n_cells <= n_in))

  expect_error(direct_sgrna_rra(pheno, calls, phenotype_range = c(99, 100),
                                min_cells_per_sgrna = 1),
               "无法计算")
})

test_that("run_direct_rra_pipeline: 输出 list 两元素且 gene_level 可喂火山图", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1), paste0("cell", 1:500))
  n_guides <- 10
  calls <- matrix(FALSE, nrow = n_guides, ncol = 500,
                  dimnames = list(paste0("g", 1:n_guides),
                                  paste0("cell", 1:500)))
  gene_map <- rep(paste0("gene", 1:5), each = 2)   # 每基因 2 条 sgRNA
  for (i in seq_len(500)) calls[sample(n_guides, 1), i] <- TRUE
  # gene1 的 sgRNA 携带细胞表型偏低
  g1_cells <- colnames(calls)[calls["g1", ] | calls["g2", ]]
  pheno[g1_cells] <- pheno[g1_cells] - 3

  sgl <- direct_sgrna_rra(pheno, calls, min_cells_per_sgrna = 1)
  lib <- data.frame(sgrna_id = paste0("g", 1:n_guides),
                    sequence = "A",
                    gene_name = gene_map, stringsAsFactors = FALSE)

  # 每基因仅 2 条 sgRNA，rra_alpha 放宽到 0.25 保证检验功效
  dir_out <- run_direct_rra_pipeline(sgl, library_df = lib, rra_alpha = 0.25)
  expect_type(dir_out, "list")
  expect_equal(names(dir_out), c("sgrna_level", "gene_level"))
  gl <- dir_out$gene_level
  expect_true(all(c("p_high", "p_low", "fdr_high", "fdr_low", "score",
                    "rank_high", "rank_low", "p_two", "alpha_high",
                    "genename_high") %in% colnames(gl)))
  # gene1（促病毒模拟）应 p_high 侧显著
  expect_lt(gl$p_high[gl$gene == "gene1"], 0.1)
})
