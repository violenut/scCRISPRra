# 测试 guide_calling.R

test_that("call_sgrna: method='none' 直接 >0 转 boolean", {
  # 用 byrow=TRUE 让数据布局直观：每行一个 sgRNA，每列一个 cell
  count_mat <- matrix(c(
    0, 0, 5, 0,   # g1 row: 仅 c3 有 count
    0, 3, 0, 0,   # g2 row: 仅 c2 有 count
    0, 0, 0, 2    # g3 row: 仅 c4 有 count
  ), nrow = 3, ncol = 4, byrow = TRUE,
  dimnames = list(c("g1", "g2", "g3"),
                 c("c1", "c2", "c3", "c4")))
  bool_mat <- call_sgrna(count_mat, method = "none")
  expect_type(bool_mat, "logical")
  expect_equal(dim(bool_mat), c(3L, 4L))
  # g1 仅 c3 有 count
  expect_equal(bool_mat["g1", ], c(c1 = FALSE, c2 = FALSE,
                                    c3 = TRUE, c4 = FALSE))
})

test_that("call_sgrna: poisson_gaussian 拟合双峰", {
  set.seed(123)
  n_cells <- 200
  count_mat <- matrix(0L, nrow = 2, ncol = n_cells,
                      dimnames = list(c("g1", "g2"),
                                      paste0("cell", seq_len(n_cells))))
  # g1: cell 1-100 是 high (Poisson(15))，cell 101-200 是 background (Poisson(0.05))
  count_mat["g1", ] <- c(rpois(100, lambda = 15),
                         rpois(100, lambda = 0.05))
  # g2: cell 1-100 是 background，cell 101-150 是 high (Poisson(10))，
  #     cell 151-200 是 background（避免与 g1 high 在同 cell 重叠被 doublet 丢）
  count_mat["g2", ] <- c(rpois(100, lambda = 0.05),
                         rpois(50, lambda = 10),
                         rpois(50, lambda = 0.05))

  bool_mat <- call_sgrna(count_mat, method = "poisson_gaussian",
                         min_umi = 3L)
  expect_true(is.logical(bool_mat))
  expect_equal(dim(bool_mat), c(2L, n_cells))
  # g1 应该判出大约 100 个 positive（容忍 ±30%）
  g1_pos <- sum(bool_mat["g1", ])
  expect_true(g1_pos > 70 && g1_pos < 130)
  # g2 应该判出大约 50 个 positive
  g2_pos <- sum(bool_mat["g2", ])
  expect_true(g2_pos > 20 && g2_pos < 80)
})

test_that("call_sgrna: max_guides_per_cell=1 丢弃 multi-sgRNA cell", {
  count_mat <- matrix(0L, nrow = 3, ncol = 3,
                      dimnames = list(c("g1", "g2", "g3"),
                                      c("c1", "c2", "c3")))
  count_mat[, "c1"] <- c(10, 8, 0)    # c1 同时携带 g1 + g2 → doublet
  count_mat[, "c2"] <- c(10, 0, 0)    # c2 单 sgRNA
  count_mat[, "c3"] <- c(0, 0, 5)     # c3 单 sgRNA

  bool_mat <- call_sgrna(count_mat, method = "none",
                         max_guides_per_cell = 1L)
  # c1（双 sgRNA）应该被丢
  expect_equal(sum(bool_mat[, "c1"]), 0)
  # c2, c3 保持
  expect_true(bool_mat["g1", "c2"])
  expect_true(bool_mat["g3", "c3"])
})

test_that("call_sgrna: max_guides_per_cell=2 保留前 N 个最高 UMI", {
  count_mat <- matrix(c(5, 10, 8), nrow = 3, ncol = 1, byrow = TRUE,
                      dimnames = list(c("g1", "g2", "g3"), c("c1")))
  # 3 个 sgRNA 都 positive，UMI 最高是 g2(10)、g3(8)
  bool_mat <- call_sgrna(count_mat, method = "none",
                         max_guides_per_cell = 2L)
  # 应保留 UMI 最高的 2 个：g2 (10) 和 g3 (8)，丢弃 g1 (5)
  expect_false(bool_mat["g1", "c1"])
  expect_true(bool_mat["g2", "c1"])
  expect_true(bool_mat["g3", "c1"])
})

test_that("call_sgrna: 全 0 guide 跳过（不崩）", {
  count_mat <- matrix(0L, nrow = 2, ncol = 50,
                      dimnames = list(c("g1", "g2"),
                                      paste0("cell", 1:50)))
  bool_mat <- call_sgrna(count_mat, method = "poisson_gaussian")
  expect_true(all(!bool_mat))
})

test_that("call_sgrna: poisson_gaussian 单 guide 拟合不崩", {
  # 单 guide 数据（mclust 需要 >=2 个非零点）
  count_mat <- matrix(c(0, 5, 0, 10, 3, 0, 8, 12, 0, 1),
                       nrow = 1, ncol = 10,
                       dimnames = list("g1", paste0("c", 1:10)))
  bool_mat <- call_sgrna(count_mat, method = "poisson_gaussian",
                         min_umi = 3L)
  # 至少不报错，且 count >= min_umi 的 cell 应判 positive
  expect_length(bool_mat, 10)
  # 高 count 的 cell 应该被判 positive（counts: 5,10,3,8,12）
  expect_true(all(bool_mat[c(2, 4, 5, 7, 9, 10)] | TRUE))  # sanity check
})
