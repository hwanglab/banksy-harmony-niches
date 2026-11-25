# Cell-type deconvolution of spots with RCTD (spacexr), given a single-cell RNA-seq reference.

.clean_names <- function(x) gsub("[^A-Za-z0-9_.+-]", "_", as.character(x))

# Reference input -> list(counts = genes x cells integer matrix, types = factor named by cell)
read_reference <- function(cfg) {
  r <- cfg$reference; fmt <- r$format
  if (fmt %in% c("seurat_rds", "seurat_qs")) {
    obj <- if (fmt == "seurat_qs") qs::qread(r$path) else readRDS(r$path)
    assay <- if (is.null(r$assay)) Seurat::DefaultAssay(obj) else r$assay
    counts <- SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
    if (is.null(r$celltype_col) || !r$celltype_col %in% names(obj[[]])) stop("reference.celltype_col not found in the Seurat metadata")
    types <- stats::setNames(obj[[r$celltype_col]][, 1], colnames(obj))
  } else if (fmt == "rds_list") {
    x <- readRDS(r$path); counts <- x$counts
    types <- if (is.data.frame(x$cell_types)) stats::setNames(x$cell_types[[2]], x$cell_types[[1]]) else x$cell_types
    if (is.null(names(types))) names(types) <- colnames(counts)
  } else stop("reference.format must be seurat_rds, seurat_qs or rds_list")
  types <- types[colnames(counts)]
  keep <- !is.na(types)
  counts <- methods::as(counts[, keep, drop = FALSE], "CsparseMatrix"); types <- .clean_names(types[keep])
  if (any(abs(counts@x - round(counts@x)) > 1e-6)) stop("the reference must hold raw integer counts (found non-integers)")
  n <- table(types); small <- names(n)[n < r$min_cells_per_type]
  if (length(small)) message(sprintf("reference: %d cell type(s) with < %d cells dropped: %s", length(small), r$min_cells_per_type, paste(small, collapse = ", ")))
  k <- !types %in% small
  list(counts = counts[, k, drop = FALSE], types = factor(stats::setNames(types[k], colnames(counts)[k])))
}

build_reference <- function(ref, cfg) {
  suppressWarnings(spacexr::Reference(ref$counts, ref$types, stats::setNames(Matrix::colSums(ref$counts), colnames(ref$counts)), n_max_cells = cfg$reference$max_cells_per_type))
}

# One sample -> spots x cell types matrix of weights that sum to 1.
run_rctd_sample <- function(x, reference, cfg) {
  p <- cfg$rctd
  counts <- x$counts; counts@x <- round(counts@x)
  coords <- data.frame(x = x$xy[, "x"], y = x$xy[, "y"], row.names = colnames(counts))
  puck <- spacexr::SpatialRNA(coords, counts, Matrix::colSums(counts))
  obj <- spacexr::create.RCTD(puck, reference, max_cores = p$cores, gene_cutoff = p$gene_cutoff, fc_cutoff = p$fc_cutoff, gene_cutoff_reg = p$gene_cutoff_reg, fc_cutoff_reg = p$fc_cutoff_reg,
                              UMI_min = p$umi_min, UMI_max = p$umi_max, counts_MIN = p$counts_min, UMI_min_sigma = p$umi_min_sigma, CELL_MIN_INSTANCE = cfg$reference$min_cells_per_type,
                              CONFIDENCE_THRESHOLD = p$confidence_threshold, DOUBLET_THRESHOLD = p$doublet_threshold)
  obj <- spacexr::run.RCTD(obj, doublet_mode = p$mode)
  types <- obj@cell_type_info$info[[2]]; res <- obj@results
  if (p$mode == "multi") {
    bc <- colnames(obj@spatialRNA@counts)[seq_along(res)]
    W <- matrix(0, length(res), length(types), dimnames = list(bc, types))
    for (i in seq_along(res)) { w <- res[[i]]$all_weights; W[i, names(w)] <- as.numeric(w) }
  } else {
    W <- as.matrix(res$weights)[, types, drop = FALSE]
  }
  W[W < 0] <- 0
  W <- W / pmax(rowSums(W), 1e-12)
  extra <- if (p$mode == "doublet") data.table::data.table(barcode = rownames(res$results_df), res$results_df) else NULL
  list(weights = W, doublet_table = extra, genes_used = length(obj@internal_vars$gene_list_reg))
}

# Optional grouping of fine cell types into broader ones (config: groups: {Tcell: [CD4T, CD8T], ...}).
group_weights <- function(W, groups) {
  g <- vapply(names(groups), function(k) rowSums(W[, intersect(.clean_names(groups[[k]]), colnames(W)), drop = FALSE]), numeric(nrow(W)))
  if (is.null(dim(g))) g <- matrix(g, nrow = nrow(W), dimnames = list(rownames(W), names(groups)))
  rest <- setdiff(colnames(W), .clean_names(unlist(groups)))
  if (length(rest)) g <- cbind(g, other = rowSums(W[, rest, drop = FALSE]))
  g
}
