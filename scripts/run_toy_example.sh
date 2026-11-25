#!/usr/bin/env bash
# End-to-end test on simulated data (about one minute). Run from the repository root.
set -euo pipefail
Rscript scripts/make_toy_data.R example
cat > example/config.yaml <<YAML
samples: example/samples.tsv
outdir: example/results
qc: {min_features: 20, max_pct_mt: 100, min_spots: 50}
hvg: {n: 300}
banksy: {lambda: 0.2, k_geom: 4, use_agf: false}
subsample: {n_per_sample: 400, seed: 1, stability_seed: 2}
cluster: {dims: 10, resolutions: [0.05, 0.3], main_resolution: 0.05}
assign: {all_spots: true, min_spots_per_cluster: 30}
report: {colour_by: [sample_id, group], marker_min_auc: 0.7}
YAML
Rscript run_pipeline.R --config example/config.yaml
Rscript - <<'RS'
a <- data.table::fread("example/results/all_spots_clusters.tsv.gz")
t <- data.table::rbindlist(lapply(unique(a$sample_id), function(s) cbind(sample_id = s, data.table::fread(file.path("example/data", paste0(s, "_truth.tsv"))))))
m <- merge(a, t, by = c("sample_id", "barcode"))
ari <- mclust::adjustedRandIndex(m$cluster, m$domain)
cat(sprintf("toy check: %d clusters for 3 simulated domains, adjusted Rand index %.2f\n", data.table::uniqueN(m$cluster), ari))
if (ari < 0.8) { cat("FAILED: expected ARI >= 0.8\n"); quit(status = 1) }
cat("PASSED\n")
RS
