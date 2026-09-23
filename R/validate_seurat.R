#' 校验 Seurat 对象是否符合 scCRISPRra 输入要求
#'
#' 检查 Seurat 对象是否同时满足以下条件：
#' \itemize{
#'   \item 是合法的 \link[Seurat]{Seurat} 对象
#'   \item 含有指定的 CRISPR Guide assay（默认名 \code{"CRISPR Guide"}）
#'   \item 该 assay 的 rownames 是 sgRNA id（用户在 cellranger feature reference
#'         中应将 sgRNA id 填入 \code{id} 字段，使 assay rownames 直接为 id）
#'   \item \code{meta.data} 中存在指定的表型评分列（连续数值）
#'   \item 可选：存在指定的 sample 列
#' }
#'
#' @param seu Seurat 对象
#' @param guide_assay CRISPR Guide assay 名，默认 \code{"CRISPR Guide"}
#' @param phenotype_col metadata 中连续表型评分的列名（必传）
#' @param sample_col 可选，metadata 中 sample 列名
#'
#' @return invisible(\code{TRUE})，校验通过；否则抛出错误
#' @export
#'
#' @examples
#' \dontrun{
#' validate_seurat(seu, phenotype_col = "virus_load")
#' }
validate_seurat <- function(seu,
                            phenotype_col,
                            guide_assay = "CRISPR Guide",
                            sample_col = NULL) {
  stopifnot("'seu' 必须是 Seurat 对象" = inherits(seu, "Seurat"))

  assays_avail <- Seurat::Assays(seu)
  if (!guide_assay %in% assays_avail) {
    stop(sprintf(
      "Seurat 对象中找不到 assay '%s'。可用 assay: %s",
      guide_assay, paste(assays_avail, collapse = ", ")
    ))
  }

  feat_names <- rownames(seu[[guide_assay]])
  if (is.null(feat_names) || length(feat_names) == 0) {
    stop(sprintf("assay '%s' 的 rownames 为空，无法用作 sgRNA id", guide_assay))
  }

  meta_cols <- colnames(seu@meta.data)
  if (missing(phenotype_col) || is.null(phenotype_col) ||
      !phenotype_col %in% meta_cols) {
    stop(sprintf(
      "metadata 中找不到表型评分列 '%s'。可用列: %s",
      phenotype_col, paste(meta_cols, collapse = ", ")
    ))
  }

  phenotype_vals <- seu@meta.data[[phenotype_col]]
  if (!is.numeric(phenotype_vals)) {
    stop(sprintf("表型评分列 '%s' 必须是 numeric，实际为 %s",
                 phenotype_col, class(phenotype_vals)[1]))
  }

  if (!is.null(sample_col)) {
    if (!sample_col %in% meta_cols) {
      stop(sprintf("metadata 中找不到 sample 列 '%s'", sample_col))
    }
  }

  invisible(TRUE)
}


#' 内部：从 Seurat 对象提取 CRISPR Guide count 矩阵
#'
#' @param seu Seurat 对象
#' @param guide_assay CRISPR Guide assay 名
#' @return 矩阵（cell × sgRNA），colnames = cell barcode，rownames = sgRNA id
#' @keywords internal
.get_guide_matrix <- function(seu, guide_assay = "CRISPR Guide") {
  ## 仅校验 Seurat + guide assay，不校验 phenotype（get_guide_matrix 与表型无关）
  stopifnot("'seu' 必须是 Seurat 对象" = inherits(seu, "Seurat"))
  assays_avail <- Seurat::Assays(seu)
  if (!guide_assay %in% assays_avail) {
    stop(sprintf(
      "Seurat 对象中找不到 assay '%s'。可用 assay: %s",
      guide_assay, paste(assays_avail, collapse = ", ")
    ))
  }

  ## Seurat 中 assay[[]] 返回 meta.features 注释表而非表达矩阵；
  ## 取 counts 必须用 GetAssayData（v5 等价 LayerData(layer="counts")）
  mat <- Seurat::GetAssayData(seu, assay = guide_assay, layer = "counts")
  if (!is.matrix(mat) && !inherits(mat, "dgCMatrix") &&
      !inherits(mat, "dgTMatrix")) {
    mat <- as.matrix(mat)
  }

  if (is.null(rownames(mat)) || length(rownames(mat)) == 0) {
    stop("guide assay 的 rownames 为空")
  }
  if (ncol(mat) == 0) {
    stop(sprintf("guide assay '%s' 的 counts 矩阵列数为 0（无 cell）", guide_assay))
  }
  mat
}


#' 内部：从 Seurat 对象提取表型评分向量
#'
#' @param seu Seurat 对象
#' @param phenotype_col metadata 中表型评分列名
#' @return 命名数值向量，names = cell barcode
#' @keywords internal
.get_phenotype_vector <- function(seu, phenotype_col) {
  ## 取表型向量与 guide assay 无关，这里只校验 Seurat + phenotype 列
  stopifnot("'seu' 必须是 Seurat 对象" = inherits(seu, "Seurat"))
  meta_cols <- colnames(seu@meta.data)
  if (missing(phenotype_col) || is.null(phenotype_col) ||
      !phenotype_col %in% meta_cols) {
    stop(sprintf(
      "metadata 中找不到表型评分列 '%s'。可用列: %s",
      phenotype_col, paste(meta_cols, collapse = ", ")
    ))
  }
  vals <- seu@meta.data[[phenotype_col]]
  if (!is.numeric(vals)) {
    stop(sprintf("表型评分列 '%s' 必须是 numeric，实际为 %s",
                 phenotype_col, class(vals)[1]))
  }
  names(vals) <- rownames(seu@meta.data)
  vals
}
