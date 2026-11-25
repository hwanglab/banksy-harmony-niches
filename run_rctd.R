#!/usr/bin/env Rscript
# Cell-type fractions per spot with RCTD, from a single-cell RNA-seq reference. Intended for Visium-type spots.
# Usage: Rscript run_rctd.R --config config/example_rctd_config.yaml [--samples samples.tsv] [--outdir dir] [--sample <sample_id>]
#   --sample runs one sample only (use it for Slurm array jobs); the summary is rebuilt from every sample already finished.
suppressPackageStartupMessages({library(data.table); library(Matrix)})
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])))
for (f in c("io.R", "steps.R", "rctd.R")) source(file.path(here, "R", f))
a <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) { i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
if ("--help" %in% a || is.null(opt("--config"))) { cat("Usage: Rscript run_rctd.R --config <yaml> [--samples <tsv>] [--outdir <dir>] [--sample <sample_id>]\n"); quit(status = 0) }
defaults <- list(
  input = list(assay = NULL, x_col = NULL, y_col = NULL),
  qc = list(min_features = 200, max_pct_mt = 20, mt_pattern = "^MT-|^mt-", min_spots = 50),
  reference = list(format = "seurat_rds", assay = NULL, celltype_col = NULL, min_cells_per_type = 25, max_cells_per_type = 10000),
  rctd = list(mode = "full", cores = 4, gene_cutoff = 0.000125, fc_cutoff = 0.5, gene_cutoff_reg = 0.0002, fc_cutoff_reg = 0.75, umi_min = 100, umi_max = 2e7,
              counts_min = 10, umi_min_sigma = 300, confidence_threshold = 5, doublet_threshold = 20),
  groups = NULL)
cfg <- utils::modifyList(defaults, yaml::read_yaml(opt("--config")))
if (!is.null(opt("--samples"))) cfg$samples <- opt("--samples")
if (!is.null(opt("--outdir"))) cfg$outdir <- opt("--outdir")
stopifnot(!is.null(cfg$samples), !is.null(cfg$outdir), !is.null(cfg$reference$path), cfg$rctd$mode %in% c("full", "multi", "doublet"))
out <- cfg$outdir; dir.create(out, showWarnings = FALSE, recursive = TRUE); yaml::write_yaml(cfg, file.path(out, "rctd_config_used.yaml"))
say <- function(...) { m <- sprintf(...); message(m); cat(m, "\n", file = file.path(out, "rctd_log.txt"), append = TRUE) }
suppressPackageStartupMessages({library(Seurat); library(spacexr)})
meta <- read_samplesheet(cfg$samples); todo <- if (is.null(opt("--sample"))) meta$sample_id else opt("--sample"); stopifnot(all(todo %in% meta$sample_id))

ref <- read_reference(cfg); reference <- build_reference(ref, cfg)
say("reference: %d cells, %d cell types (%s)", ncol(ref$counts), nlevels(ref$types), paste(sprintf("%s %d", names(table(ref$types)), as.integer(table(ref$types))), collapse = ", "))
if (cfg$rctd$mode != "full") say("note: mode '%s' assumes few cell types per spot; 'full' is the usual choice for 55 um Visium spots", cfg$rctd$mode)

for (s in todo) {
  x <- qc_sample(read_sample(meta[sample_id == s], cfg), cfg)
  if (x$n_pass < cfg$qc$min_spots) { say("%s skipped: %d spots pass QC (< %d)", s, x$n_pass, cfg$qc$min_spots); next }
  t0 <- Sys.time(); r <- run_rctd_sample(x, reference, cfg); W <- r$weights
  fwrite(data.table(barcode = rownames(W), round(W, 5)), file.path(out, paste0(s, "_weights.tsv.gz")), sep = "\t")
  if (!is.null(cfg$groups)) fwrite(data.table(barcode = rownames(W), round(group_weights(W, cfg$groups), 5)), file.path(out, paste0(s, "_group_weights.tsv.gz")), sep = "\t")
  if (!is.null(r$doublet_table)) fwrite(r$doublet_table, file.path(out, paste0(s, "_doublet_calls.tsv.gz")), sep = "\t")
  say("%s: %d of %d QC-passing spots deconvolved (%d genes) in %.1f min | mean largest weight %.2f", s, nrow(W), x$n_pass, r$genes_used, as.numeric(difftime(Sys.time(), t0, units = "mins")), mean(apply(W, 1, max)))
}

# summary over every finished sample: mean weight per cell type, and the share of spots each type dominates
fs <- list.files(out, pattern = "_weights\\.tsv\\.gz$", full.names = TRUE); fs <- fs[!grepl("_group_weights", fs)]
S <- rbindlist(lapply(fs, function(f) { w <- as.matrix(fread(f), rownames = 1); data.table(sample_id = sub("_weights\\.tsv\\.gz$", "", basename(f)), spots = nrow(w), cell_type = colnames(w), mean_weight = round(colMeans(w), 4),
  dominant_share = round(as.numeric(table(factor(colnames(w)[max.col(w, ties.method = "first")], levels = colnames(w)))) / nrow(w), 4)) }))
fwrite(dcast(S, sample_id + spots ~ cell_type, value.var = "mean_weight"), file.path(out, "mean_weight_per_sample.tsv"), sep = "\t")
fwrite(dcast(S, sample_id + spots ~ cell_type, value.var = "dominant_share"), file.path(out, "dominant_type_share_per_sample.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(out, "rctd_sessionInfo.txt")); say("summary written for %d sample(s): %s", length(fs), normalizePath(out))
