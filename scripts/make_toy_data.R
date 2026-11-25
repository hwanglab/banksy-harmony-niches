#!/usr/bin/env Rscript
# Simulate a small data set (4 samples, 3 spatial domains each, with a sample-specific offset) to test the pipeline.
# Usage: Rscript scripts/make_toy_data.R [outdir]     (default: example)
out <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(out)) out <- "example"
dir.create(file.path(out, "data"), recursive = TRUE, showWarnings = FALSE); set.seed(42)
n_genes <- 600; genes <- sprintf("GENE%03d", seq_len(n_genes)); side <- 30
prog <- list(1:60, 61:120, 121:180)                           # marker genes of the three domains
rows <- list()
for (s in sprintf("toy%d", 1:4)) {
  g <- expand.grid(x = seq_len(side), y = seq_len(side))
  cut1 <- sample(9:11, 1); cut2 <- sample(19:21, 1)            # domain borders differ between samples
  dom <- ifelse(g$x <= cut1, 1L, ifelse(g$x <= cut2, 2L, 3L))
  base <- rgamma(n_genes, shape = 0.6, rate = 0.3) * exp(rnorm(n_genes, 0, 0.35))   # second factor = sample (batch) offset
  mu <- sapply(dom, function(d) { m <- base; m[prog[[d]]] <- m[prog[[d]]] * 6 + 3; m })
  depth <- rlnorm(nrow(g), 0, 0.3)
  counts <- matrix(rpois(length(mu), t(t(mu) * depth)), nrow = n_genes, dimnames = list(genes, sprintf("spot%04d", seq_len(nrow(g)))))
  coords <- data.frame(x = g$x * 100, y = g$y * 100, row.names = colnames(counts))
  f <- file.path(out, "data", paste0(s, ".rds"))
  saveRDS(list(counts = Matrix::Matrix(counts, sparse = TRUE), coords = coords), f)
  write.table(data.frame(barcode = colnames(counts), domain = dom), file.path(out, "data", paste0(s, "_truth.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  rows[[s]] <- data.frame(sample_id = s, path = f, format = "rds_list", group = ifelse(s %in% c("toy1", "toy2"), "A", "B"))
}
write.table(do.call(rbind, rows), file.path(out, "samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("toy data written to", out, "\n")
