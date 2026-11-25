#!/usr/bin/env bash
# End-to-end test of the RCTD step on simulated mixtures (a few minutes). Run from the repository root.
set -euo pipefail
Rscript scripts/make_toy_rctd_data.R example_rctd
cat > example_rctd/config.yaml <<YAML
samples: example_rctd/samples.tsv
outdir: example_rctd/results
reference: {path: example_rctd/data/reference.rds, format: rds_list, min_cells_per_type: 25}
qc: {min_features: 20, max_pct_mt: 100, min_spots: 50}
rctd: {mode: full, cores: ${RCTD_CORES:-4}}
groups: {AB: [TypeA, TypeB]}
YAML
Rscript run_rctd.R --config example_rctd/config.yaml
Rscript - <<'RS'
r <- sapply(c("mix1", "mix2"), function(s) { w <- as.matrix(data.table::fread(sprintf("example_rctd/results/%s_weights.tsv.gz", s)), rownames = 1)
  t <- as.matrix(data.table::fread(sprintf("example_rctd/data/%s_truth.tsv", s)), rownames = 1); k <- intersect(rownames(w), rownames(t)); sapply(colnames(t), function(ct) cor(w[k, ct], t[k, ct])) })
print(round(r, 2)); cat(sprintf("toy check: lowest correlation between estimated and true fraction %.2f\n", min(r)))
if (min(r) < 0.9) { cat("FAILED: expected >= 0.9\n"); quit(status = 1) }
cat("PASSED\n")
RS
