# Input readers. Every reader returns list(counts = dgCMatrix genes x spots, xy = matrix[spot, c("x","y")]).

.read_positions <- function(spatial_dir) {
  cand <- file.path(spatial_dir, c("tissue_positions.csv", "tissue_positions_list.csv", "tissue_positions.parquet"))
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) stop("no tissue_positions file in ", spatial_dir)
  if (grepl("parquet$", f)) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("tissue_positions.parquet needs the 'arrow' package")
    pos <- data.table::as.data.table(arrow::read_parquet(f))
  } else {
    pos <- data.table::fread(f, header = grepl("^barcode", readLines(f, n = 1)))
  }
  data.table::setnames(pos, 1:6, c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"))
  pos
}

.read_coords_file <- function(f) {
  d <- data.table::fread(f)
  need <- c("barcode", "x", "y")
  if (!all(need %in% names(d))) stop("coords file must have columns barcode, x, y: ", f)
  d
}

read_sample <- function(row, cfg) {
  fmt <- row$format; path <- row$path
  if (fmt %in% c("spaceranger", "10x_h5")) {
    h5 <- if (fmt == "spaceranger") file.path(path, "filtered_feature_bc_matrix.h5") else path
    counts <- Seurat::Read10X_h5(h5)
    if (is.list(counts)) counts <- counts[[1]]
    pos <- .read_positions(file.path(dirname(h5), "spatial"))
    pos <- pos[match(colnames(counts), pos$barcode)]
    xy <- cbind(x = pos$pxl_col, y = pos$pxl_row)
  } else if (fmt %in% c("seurat_qs", "seurat_rds")) {
    obj <- if (fmt == "seurat_qs") {
      if (!requireNamespace("qs", quietly = TRUE)) stop("format seurat_qs needs the 'qs' package")
      qs::qread(path)
    } else readRDS(path)
    assay <- if (is.null(cfg$input$assay)) Seurat::DefaultAssay(obj) else cfg$input$assay
    counts <- SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
    md <- obj[[]]
    xc <- cfg$input$x_col; yc <- cfg$input$y_col
    if (!is.null(xc) && !is.null(yc) && all(c(xc, yc) %in% names(md))) {
      xy <- cbind(x = md[[xc]], y = md[[yc]])
    } else {
      tc <- Seurat::GetTissueCoordinates(obj)
      tc <- tc[colnames(counts), ]
      num <- names(tc)[vapply(tc, is.numeric, logical(1))]
      xy <- cbind(x = tc[[num[1]]], y = tc[[num[2]]])
    }
  } else if (fmt == "rds_list") {
    x <- readRDS(path)
    counts <- x$counts
    xy <- as.matrix(x$coords[colnames(counts), c("x", "y")])
  } else stop("unknown format '", fmt, "' (use spaceranger, 10x_h5, seurat_qs, seurat_rds or rds_list)")
  if (!is.null(row$coords) && !is.na(row$coords) && nzchar(row$coords)) {
    d <- .read_coords_file(row$coords)
    d <- d[match(colnames(counts), d$barcode)]
    xy <- cbind(x = d$x, y = d$y)
  }
  counts <- methods::as(counts, "CsparseMatrix")
  rownames(xy) <- colnames(counts)
  keep <- stats::complete.cases(xy)
  if (!all(keep)) message(sprintf("  %s: %d spots without coordinates dropped", row$sample_id, sum(!keep)))
  list(counts = counts[, keep, drop = FALSE], xy = xy[keep, , drop = FALSE])
}

read_samplesheet <- function(f) {
  s <- data.table::fread(f, colClasses = "character")
  need <- c("sample_id", "path", "format")
  if (!all(need %in% names(s))) stop("sample sheet needs columns: ", paste(need, collapse = ", "))
  if (anyDuplicated(s$sample_id)) stop("duplicated sample_id in the sample sheet")
  if (any(grepl("|", s$sample_id, fixed = TRUE))) stop("sample_id must not contain '|'")
  miss <- s$path[!file.exists(s$path)]
  if (length(miss)) stop("input not found: ", paste(miss, collapse = ", "))
  s
}
