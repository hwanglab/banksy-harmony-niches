# Pipeline steps: QC -> depth equalisation -> HVG -> BANKSY -> subsample -> BANKSY PCA -> Harmony -> clustering.

qc_sample <- function(x, cfg) {
  n_feat <- Matrix::colSums(x$counts > 0)
  mt <- grepl(cfg$qc$mt_pattern, rownames(x$counts))
  tot <- Matrix::colSums(x$counts)
  pct_mt <- if (any(mt)) 100 * Matrix::colSums(x$counts[mt, , drop = FALSE]) / pmax(tot, 1) else rep(0, ncol(x$counts))
  keep <- n_feat > cfg$qc$min_features & pct_mt < cfg$qc$max_pct_mt
  list(counts = x$counts[, keep, drop = FALSE], xy = x$xy[keep, , drop = FALSE], n_in = length(keep), n_pass = sum(keep))
}

# Binomial thinning so that every sample has the median depth of the shallowest sample.
equalise_depth <- function(objs, seed = 1) {
  med <- vapply(objs, function(x) stats::median(Matrix::colSums(x$counts)), numeric(1))
  target <- min(med)
  set.seed(seed)
  for (s in names(objs)) {
    p <- min(1, target / med[[s]])
    if (p < 1) {
      m <- objs[[s]]$counts
      m@x <- as.numeric(stats::rbinom(length(m@x), round(m@x), p))
      objs[[s]]$counts <- Matrix::drop0(m)
    }
  }
  attr(objs, "depth_before") <- med
  attr(objs, "depth_target") <- target
  objs
}

pick_spots <- function(objs, n, seed) {
  set.seed(seed)
  lapply(objs, function(x) if (ncol(x$counts) > n) sort(sample(ncol(x$counts), n)) else seq_len(ncol(x$counts)))
}

pooled_counts <- function(objs, genes, idx) {
  do.call(cbind, lapply(names(objs), function(s) {
    m <- objs[[s]]$counts[genes, idx[[s]], drop = FALSE]
    colnames(m) <- paste(s, colnames(m), sep = "|")
    m
  }))
}

select_hvg <- function(objs, genes, idx, n) {
  so <- Seurat::CreateSeuratObject(pooled_counts(objs, genes, idx))
  so <- Seurat::NormalizeData(so, verbose = FALSE)
  Seurat::VariableFeatures(Seurat::FindVariableFeatures(so, nfeatures = n, verbose = FALSE))
}

# BANKSY neighbourhood features, computed on the WHOLE sample so that a subsampled spot still sees its true neighbours.
banksy_features <- function(x, sample_id, hvg, cfg) {
  ls <- Matrix::colSums(x$counts)
  nc <- as.matrix(Matrix::t(Matrix::t(x$counts[hvg, , drop = FALSE]) / ls) * stats::median(ls))
  se <- SpatialExperiment::SpatialExperiment(
    assays = list(normcounts = nc), spatialCoords = x$xy,
    colData = S4Vectors::DataFrame(sample_id = rep(sample_id, ncol(nc))))
  Banksy::computeBanksy(se, assay_name = "normcounts", compute_agf = cfg$banksy$use_agf, k_geom = cfg$banksy$k_geom, verbose = FALSE)
}

build_embedding <- function(objs, bk, genes, meta, cfg, seed) {
  idx <- pick_spots(objs, cfg$subsample$n_per_sample, seed)
  se <- do.call(cbind, lapply(names(objs), function(s) {
    z <- bk[[s]][, idx[[s]]]
    colnames(z) <- paste(s, colnames(z), sep = "|")
    z
  }))
  se <- Banksy::runBanksyPCA(se, use_agf = cfg$banksy$use_agf, lambda = cfg$banksy$lambda, group = "sample_id", npcs = cfg$pca$npcs, seed = 1)
  pcs <- SingleCellExperiment::reducedDim(se, grep("^PCA", SingleCellExperiment::reducedDimNames(se), value = TRUE)[1])
  rownames(pcs) <- colnames(se)
  so <- Seurat::CreateSeuratObject(pooled_counts(objs, genes, idx))
  so$sample_id <- sub("\\|.*", "", colnames(so))
  for (v in setdiff(names(meta), c("sample_id", "path", "format", "coords"))) so[[v]] <- meta[[v]][match(so$sample_id, meta$sample_id)]
  so <- Seurat::NormalizeData(so, verbose = FALSE)
  emb <- pcs[colnames(so), , drop = FALSE]
  vars <- cfg$harmony$vars
  if (length(vars) && !identical(vars, "none")) {
    emb <- harmony::RunHarmony(emb, so[[]], vars, verbose = FALSE)
    rownames(emb) <- colnames(so)
  }
  colnames(emb) <- paste0("integrated_", seq_len(ncol(emb)))
  so[["integrated"]] <- Seurat::CreateDimReducObject(emb, key = "integrated_", assay = Seurat::DefaultAssay(so))
  dims <- seq_len(min(cfg$cluster$dims, ncol(emb)))
  so <- Seurat::FindNeighbors(so, reduction = "integrated", dims = dims, verbose = FALSE)
  for (r in cfg$cluster$resolutions) so <- Seurat::FindClusters(so, resolution = r, random.seed = 1, verbose = FALSE)
  so
}

cluster_col <- function(so, res) {
  cand <- grep("_snn_res\\.", names(so[[]]), value = TRUE)
  hit <- cand[abs(as.numeric(sub(".*_snn_res\\.", "", cand)) - res) < 1e-8]
  if (!length(hit)) stop("resolution ", res, " was not computed")
  hit[1]
}

# ---- carry the labels to every spot ---------------------------------------------------------------------------
# Feature per spot: sqrt(1 - lambda) * z(own expression) and sqrt(lambda) * z(mean of k neighbours), z within sample.
spot_features <- function(x, hvg, lambda, k) {
  tot <- Matrix::colSums(x$counts)
  L <- log1p(Matrix::t(Matrix::t(x$counts[hvg, , drop = FALSE]) / tot) * 1e4)
  nn <- FNN::get.knn(x$xy, k = min(k, nrow(x$xy) - 1))$nn.index
  own <- as.matrix(L)
  nb <- Reduce(`+`, lapply(seq_len(ncol(nn)), function(j) own[, nn[, j], drop = FALSE])) / ncol(nn)
  z <- function(m) { s <- apply(m, 1, stats::sd); (m - rowMeans(m)) / ifelse(s > 0, s, 1) }
  f <- rbind(sqrt(1 - lambda) * z(own), sqrt(lambda) * z(nb))
  colnames(f) <- colnames(x$counts)
  f
}

assign_all_spots <- function(objs, hvg, labels, cfg) {
  feats <- lapply(objs, spot_features, hvg = hvg, lambda = cfg$banksy$lambda, k = cfg$banksy$k_geom)
  keep <- names(which(table(labels$cluster) >= cfg$assign$min_spots_per_cluster))
  centroid <- function(samples) vapply(keep, function(k) {
    rowMeans(do.call(cbind, lapply(samples, function(s) feats[[s]][, labels[sample_id == s & cluster == k]$barcode, drop = FALSE])))
  }, numeric(nrow(feats[[1]])))
  nearest <- function(f, cen) {
    r <- stats::cor(f, cen)
    list(cluster = colnames(cen)[max.col(r, ties.method = "first")], cor = apply(r, 1, max))
  }
  samples <- names(objs)
  loso <- data.table::rbindlist(lapply(samples, function(s) {
    b <- labels[sample_id == s & cluster %in% keep]
    if (!nrow(b) || length(samples) < 2) return(NULL)
    data.table::data.table(sample_id = s, truth = b$cluster, predicted = nearest(feats[[s]][, b$barcode, drop = FALSE], centroid(setdiff(samples, s)))$cluster)
  }))
  cen <- centroid(samples)
  all <- data.table::rbindlist(lapply(samples, function(s) {
    a <- nearest(feats[[s]], cen)
    data.table::data.table(sample_id = s, barcode = colnames(feats[[s]]), x = objs[[s]]$xy[, "x"], y = objs[[s]]$xy[, "y"], cluster = a$cluster, cor_to_centroid = round(a$cor, 3))
  }))
  list(all = all, loso = loso, kept = keep)
}
