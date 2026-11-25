#!/usr/bin/env Rscript
# Simulate a single-cell reference (4 cell types) and 2 samples of spots that are known mixtures of those types.
# Usage: Rscript scripts/make_toy_rctd_data.R [outdir]     (default: example_rctd)
out <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(out)) out <- "example_rctd"
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE); set.seed(7)
n_genes <- 800; genes <- sprintf("GENE%03d", seq_len(n_genes)); types <- c("TypeA", "TypeB", "TypeC", "TypeD")
base <- rgamma(n_genes, shape = 0.5, rate = 0.5) + 0.02
prof <- sapply(seq_along(types), function(k) { p <- base; i <- ((k - 1) * 80 + 1):(k * 80); p[i] <- p[i] * 8 + 1; p / sum(p) }); colnames(prof) <- types
cells <- rep(types, each = 150); umi <- round(rlnorm(length(cells), log(3000), 0.3))
ref <- sapply(seq_along(cells), function(i) rmultinom(1, umi[i], prof[, cells[i]])); dimnames(ref) <- list(genes, sprintf("cell%04d", seq_along(cells)))
saveRDS(list(counts = Matrix::Matrix(ref, sparse = TRUE), cell_types = setNames(cells, colnames(ref))), file.path(out, "data", "reference.rds"))
rows <- list()
for (s in c("mix1", "mix2")) {
  side <- 16; g <- expand.grid(x = seq_len(side), y = seq_len(side)); n <- nrow(g)
  a <- matrix(rgamma(n * 4, shape = 0.5), n, 4); a[, 1] <- a[, 1] * (1 + 3 * (g$x <= side / 2)); truth <- a / rowSums(a); colnames(truth) <- types
  depth <- round(rlnorm(n, log(4000), 0.3)); eff <- exp(rnorm(n_genes, 0, 0.3))                       # eff = platform effect, the thing RCTD has to absorb
  cnt <- sapply(seq_len(n), function(i) { p <- as.numeric(prof %*% truth[i, ]) * eff; rmultinom(1, depth[i], p / sum(p)) }); dimnames(cnt) <- list(genes, sprintf("spot%04d", seq_len(n)))
  f <- file.path(out, "data", paste0(s, ".rds")); saveRDS(list(counts = Matrix::Matrix(cnt, sparse = TRUE), coords = data.frame(x = g$x * 100, y = g$y * 100, row.names = colnames(cnt))), f)
  write.table(data.frame(barcode = colnames(cnt), round(truth, 4)), file.path(out, "data", paste0(s, "_truth.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  rows[[s]] <- data.frame(sample_id = s, path = f, format = "rds_list")
}
write.table(do.call(rbind, rows), file.path(out, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); cat("toy RCTD data written to", out, "\n")
