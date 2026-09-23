# 测试 validate_seurat.R

# Seurat 对象 mock：用一个轻量的 S4 类测试校验逻辑
# 因 Seurat 对象构造复杂，这里只测错误路径（错误输入）

test_that("validate_seurat: 非 Seurat 对象报错", {
  not_seu <- list(a = 1)
  expect_error(
    validate_seurat(not_seu, phenotype_col = "x"),
    "必须是 Seurat"
  )
})

test_that("compute_breaks: 非命名向量报错", {
  expect_error(
    compute_breaks(rnorm(100)),
    "命名数值"
  )
})

test_that("compute_breaks: 长度 2 quantile 校验", {
  pheno <- setNames(rnorm(100), paste0("cell", 1:100))
  expect_error(
    compute_breaks(pheno, phenotype_quantile = c(0.05)),
    "长度必须为 2"
  )
})
