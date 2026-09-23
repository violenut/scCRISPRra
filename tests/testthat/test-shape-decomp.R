# 测试 shape_decomp.R: compute_breaks / compute_global_hist / decompose_shape / summarize_sgrna_shapes

test_that("compute_breaks: 单峰数据生成正确区间", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1),
                    paste0("cell", 1:500))
  brks <- compute_breaks(pheno,
                         phenotype_quantile = c(0.01, 0.99),
                         break_step = 0.5)
  expect_type(brks, "double")
  expect_true(length(brks) >= 2)
  expect_true(all(diff(brks) > 0))
  expect_true(min(brks) >= quantile(pheno, 0.01) - 0.01)
  expect_true(max(brks) <= quantile(pheno, 0.99) + 0.01)
})

test_that("compute_breaks: 0.01 步长产生密度合理的 bin 数", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1),
                    paste0("cell", 1:500))
  brks <- compute_breaks(pheno,
                         break_step = 0.01)
  # [q01, q99] 范围约 2.5，0.01 步长约 250 个 bin
  expect_true(length(brks) > 100)
  expect_true(length(brks) < 500)
})

test_that("compute_global_hist: 归一化后 sum(freq) == 1", {
  set.seed(123)
  pheno <- setNames(rnorm(200, mean = 10, sd = 1),
                    paste0("cell", 1:200))
  brks <- compute_breaks(pheno, break_step = 0.5)
  gh <- compute_global_hist(pheno, brks)
  expect_equal(sum(gh$freq), 1, tolerance = 1e-10)
  expect_equal(nrow(gh), length(brks) - 1)
})

test_that("decompose_shape: 右偏 sgRNA linear 系数 > 0", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1),
                    paste0("cell", 1:500))
  brks <- compute_breaks(pheno, break_step = 0.5)
  gh <- compute_global_hist(pheno, brks)

  # 构造右偏 sgRNA 的 cell 表型（shift +2）
  sub_pheno <- setNames(pheno[1:50] + 2, names(pheno)[1:50])
  shape <- decompose_shape(sub_pheno, gh, degree = 2)
  expect_length(shape, 2)
  expect_named(shape, c("linear", "quadratic"))
  expect_gt(shape["linear"], 0)
})

test_that("decompose_shape: 左偏 sgRNA linear 系数 < 0", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1),
                    paste0("cell", 1:500))
  brks <- compute_breaks(pheno, break_step = 0.5)
  gh <- compute_global_hist(pheno, brks)

  sub_pheno <- setNames(pheno[1:50] - 2, names(pheno)[1:50])
  shape <- decompose_shape(sub_pheno, gh, degree = 2)
  expect_lt(shape["linear"], 0)
})

test_that("decompose_shape: 空 sgRNA 不崩（freq=0 处理）", {
  set.seed(123)
  pheno <- setNames(rnorm(100, mean = 10, sd = 1),
                    paste0("cell", 1:100))
  brks <- compute_breaks(pheno, break_step = 0.5)
  gh <- compute_global_hist(pheno, brks)

  # 全部 cell 表型都落在 breaks 之外（极端值）
  sub_pheno <- setNames(c(-100, -99), c("a", "b"))
  shape <- decompose_shape(sub_pheno, gh, degree = 2)
  expect_length(shape, 2)
  # 全 0 的 shape
  expect_equal(as.numeric(shape), c(0, 0))
})

test_that("summarize_sgrna_shapes: 默认 min_cells_per_sgrna=20", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1),
                    paste0("cell", 1:500))
  calls <- mock_calls(n_cells = 500, n_guides = 20)

  sgrna_summary <- summarize_sgrna_shapes(pheno, calls,
                                          break_step = 0.5,
                                          poly_degree = 2)
  expect_s3_class(sgrna_summary, "data.frame")
  expect_true(all(c("sgrna", "linear", "quadratic", "freq") %in%
                  colnames(sgrna_summary)))
  # 默认 min_cells_per_sgrna = 20，所有保留的 sgRNA 至少有 20 cell
  expect_true(all(sgrna_summary$freq >= 20))
})

test_that("summarize_sgrna_shapes: poly_degree=3 输出 cubic 列", {
  set.seed(123)
  pheno <- setNames(rnorm(500, mean = 10, sd = 1),
                    paste0("cell", 1:500))
  calls <- mock_calls(n_cells = 500, n_guides = 20)

  sgrna_summary <- summarize_sgrna_shapes(pheno, calls,
                                          poly_degree = 3,
                                          break_step = 0.5)
  expect_true("cubic" %in% colnames(sgrna_summary))
})
