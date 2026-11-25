#!/usr/bin/env Rscript
# BANKSY -> Harmony -> unsupervised clustering of spatial transcriptomics spots across many samples.
# Usage: Rscript run_pipeline.R --config config/example_config.yaml [--samples samples.tsv] [--outdir results]
suppressPackageStartupMessages({library(data.table); library(Matrix)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
for (f in c("io.R", "steps.R", "report.R")) source(file.path(here, "R", f))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
if ("--help" %in% a || is.null(opt("--config"))) { cat("Usage: Rscript run_pipeline.R --config <yaml> [--samples <tsv>] [--outdir <dir>]\n"); quit(status = 0) }

defaults <- list(
  input = list(assay = NULL, x_col = NULL, y_col = NULL),
  qc = list(min_features = 200, max_pct_mt = 20, mt_pattern = "^MT-|^mt-", min_spots = 50),
  depth = list(equalise = FALSE, min_umi_after = 100),
  hvg = list(n = 2000),
  banksy = list(lambda = 0.2, k_geom = 6, use_agf = TRUE),
  subsample = list(n_per_sample = 500, seed = 1, stability_seed = 2),
  pca = list(npcs = 30),
  harmony = list(vars = "sample_id"),
  cluster = list(dims = 20, resolutions = c(0.3, 0.6, 1.0), main_resolution = 0.6),
  assign = list(all_spots = TRUE, min_spots_per_cluster = 100),
  annotation = list(weights_pattern = NULL),
  report = list(marker_min_auc = 0.65, n_markers = 10, min_spots_per_sample = 10, colour_by = "sample_id", max_spatial_plots = 50))
cfg <- utils::modifyList(defaults, yaml::read_yaml(opt("--config")))
if (!is.null(opt("--samples"))) cfg$samples <- opt("--samples")
if (!is.null(opt("--outdir"))) cfg$outdir <- opt("--outdir")
stopifnot(!is.null(cfg$samples), !is.null(cfg$outdir), cfg$cluster$main_resolution %in% cfg$cluster$resolutions)
out <- cfg$outdir; dir.create(out, showWarnings = FALSE, recursive = TRUE)
yaml::write_yaml(cfg, file.path(out, "config_used.yaml"))
logf <- file(file.path(out, "log.txt"), "w"); say <- function(...) { m <- sprintf(...); message(m); writeLines(m, logf); flush(logf) }

suppressPackageStartupMessages({library(Seurat); library(SpatialExperiment); library(Banksy); library(harmony)})
meta <- read_samplesheet(cfg$samples)
say("samples: %d | lambda %.2f | k_geom %d | Harmony on: %s", nrow(meta), cfg$banksy$lambda, cfg$banksy$k_geom, paste(cfg$harmony$vars, collapse = ", "))

# 1. read + QC ----------------------------------------------------------------------------------------------------
objs <- list(); qc <- list()
for (i in seq_len(nrow(meta))) {
  s <- meta$sample_id[i]
  x <- qc_sample(read_sample(meta[i], cfg), cfg)
  qc[[s]] <- data.table(sample_id = s, spots_in = x$n_in, spots_pass = x$n_pass, median_umi = if (x$n_pass) median(Matrix::colSums(x$counts)) else NA_real_, used = x$n_pass >= cfg$qc$min_spots)
  if (x$n_pass >= cfg$qc$min_spots) objs[[s]] <- x[c("counts", "xy")] else say("  %s dropped: %d spots pass QC (< %d)", s, x$n_pass, cfg$qc$min_spots)
}
qc <- rbindlist(qc); fwrite(qc, file.path(out, "sample_qc.tsv"), sep = "\t")
stopifnot(length(objs) >= 1); meta <- meta[sample_id %in% names(objs)]
say("samples used: %d | spots per sample (median): %.0f", length(objs), median(qc[used == TRUE]$spots_pass))

# 2. depth equalisation (optional) ------------------------------------------------------------------------------------
if (isTRUE(cfg$depth$equalise)) {
  objs <- equalise_depth(objs, seed = cfg$subsample$seed)
  say("depth equalised: median UMI per spot %.0f to %.0f before, thinned to %.0f", min(attr(objs, "depth_before")), max(attr(objs, "depth_before")), attr(objs, "depth_target"))
  for (s in names(objs)) { k <- Matrix::colSums(objs[[s]]$counts) >= cfg$depth$min_umi_after; objs[[s]]$counts <- objs[[s]]$counts[, k, drop = FALSE]; objs[[s]]$xy <- objs[[s]]$xy[k, , drop = FALSE] }
}

# 3. variable genes, 4. BANKSY features ---------------------------------------------------------------------------------
genes <- Reduce(intersect, lapply(objs, function(x) rownames(x$counts)))
hvg <- select_hvg(objs, genes, pick_spots(objs, cfg$subsample$n_per_sample, cfg$subsample$seed), cfg$hvg$n)
say("genes shared by all samples: %d | variable genes: %d", length(genes), length(hvg))
bk <- lapply(setNames(names(objs), names(objs)), function(s) banksy_features(objs[[s]], s, hvg, cfg))
say("BANKSY features computed")

# 5-8. subsample, BANKSY PCA, Harmony, clustering; second draw for stability ---------------------------------------------
so <- build_embedding(objs, bk, genes, meta, cfg, cfg$subsample$seed)
so <- Seurat::RunUMAP(so, reduction = "integrated", dims = seq_len(min(cfg$cluster$dims, ncol(Embeddings(so, "integrated")))), verbose = FALSE)
cols <- vapply(cfg$cluster$resolutions, function(r) cluster_col(so, r), character(1))
say("spots in the build: %d | clusters per resolution: %s", ncol(so), paste(sprintf("%s = %d", cfg$cluster$resolutions, vapply(cols, function(cc) nlevels(so[[cc]][, 1]), integer(1))), collapse = ", "))
stab <- list()
if (!is.null(cfg$subsample$stability_seed) && requireNamespace("mclust", quietly = TRUE)) {
  so2 <- build_embedding(objs, bk, genes, meta, cfg, cfg$subsample$stability_seed)
  sh <- intersect(colnames(so), colnames(so2))
  for (j in seq_along(cols)) stab[[j]] <- data.table(resolution = cfg$cluster$resolutions[j], shared_spots = length(sh),
    ARI_between_draws = if (length(sh) >= 50) round(mclust::adjustedRandIndex(so[[cols[j]]][sh, 1], so2[[cols[j]]][sh, 1]), 3) else NA_real_)
  stab <- rbindlist(stab); fwrite(stab, file.path(out, "stability.tsv"), sep = "\t")
  say("stability (ARI between two random draws, on shared spots): %s", paste(sprintf("res %s = %s", stab$resolution, stab$ARI_between_draws), collapse = ", "))
  rm(so2)
}

# 9. tables and figures -----------------------------------------------------------------------------------------------------
main_col <- cluster_col(so, cfg$cluster$main_resolution)
res <- write_cluster_tables(so, main_col, out, cfg); so <- res$so
cl <- data.table(spot = colnames(so), sample_id = so$sample_id, barcode = sub("^[^|]*\\|", "", colnames(so)), cluster = so$cluster)
for (j in seq_along(cols)) cl[[paste0("res_", cfg$cluster$resolutions[j])]] <- as.character(so[[cols[j]]][, 1])
um <- Embeddings(so, "umap"); cl$umap_1 <- um[, 1]; cl$umap_2 <- um[, 2]
xy <- do.call(rbind, lapply(names(objs), function(s) { m <- objs[[s]]$xy; rownames(m) <- paste(s, rownames(m), sep = "|"); m })); cl$x <- xy[cl$spot, "x"]; cl$y <- xy[cl$spot, "y"]
fwrite(cl, file.path(out, "clusters.tsv"), sep = "\t")
plot_umap(so, out, cfg$report$colour_by)
say("main resolution %s: %d clusters; %d present (>= %d spots) in at least half of the samples", cfg$cluster$main_resolution, nrow(res$summary), sum(res$summary$samples_with_min_spots >= length(objs) / 2), cfg$report$min_spots_per_sample)

# 10. carry the labels to every spot ---------------------------------------------------------------------------------------
if (isTRUE(cfg$assign$all_spots)) {
  a <- assign_all_spots(objs, hvg, cl[, .(sample_id, barcode, cluster)], cfg)
  fwrite(a$all, file.path(out, "all_spots_clusters.tsv.gz"), sep = "\t")
  if (nrow(a$loso)) {
    acc <- a$loso[, .(spots = .N, agreement = round(mean(truth == predicted), 3)), by = sample_id]; fwrite(acc, file.path(out, "assignment_check.tsv"), sep = "\t")
    say("labels carried to all %d spots (%d clusters with >= %d labelled spots) | leave-one-sample-out agreement with the clustering: %.2f (chance %.2f)", nrow(a$all), length(a$kept), cfg$assign$min_spots_per_cluster, mean(a$loso$truth == a$loso$predicted), 1 / length(a$kept))
  }
  sh <- a$all[, .N, by = .(sample_id, cluster)]; sh[, share := N / sum(N), by = sample_id]; fwrite(dcast(sh, sample_id ~ cluster, value.var = "share", fill = 0), file.path(out, "sample_share_all_spots.tsv"), sep = "\t")
  plot_spatial(a$all, out, cfg$report$max_spatial_plots)
} else plot_spatial(cl, out, cfg$report$max_spatial_plots)
# 11. optional: describe the niches by cell-type weights from run_rctd.R ---------------------------------------------------
if (!is.null(cfg$annotation$weights_pattern)) {
  lab <- if (isTRUE(cfg$assign$all_spots)) a$all[, .(sample_id, barcode, cluster)] else cl[, .(sample_id, barcode, cluster)]
  W <- rbindlist(lapply(names(objs), function(s) { f <- gsub("{sample_id}", s, cfg$annotation$weights_pattern, fixed = TRUE); if (!file.exists(f)) { say("  no cell-type weights for %s (%s)", s, f); return(NULL) }; cbind(sample_id = s, fread(f)) }), fill = TRUE)
  if (nrow(W)) { m <- merge(lab, W, by = c("sample_id", "barcode")); ct <- setdiff(names(W), c("sample_id", "barcode"))
    cw <- m[, c(list(spots = .N), lapply(.SD, function(v) round(mean(v, na.rm = TRUE), 3))), by = cluster, .SDcols = ct]; cw[, top_cell_type := ct[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = ct]
    fwrite(cw[order(-spots)], file.path(out, "cluster_celltype_weights.tsv"), sep = "\t"); say("niches described by cell-type weights: %d spots matched, %d cell types", nrow(m), length(ct)) }
}
saveRDS(so, file.path(out, "seurat_build.rds"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
say("done: %s", normalizePath(out)); close(logf)
