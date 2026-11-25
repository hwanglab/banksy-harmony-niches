# banksy-harmony-niches

Two steps for multi-sample spatial transcriptomics (Visium-type spots):

1. **`run_pipeline.R`** - find tissue niches that are shared across many samples (BANKSY -> Harmony -> clustering).
2. **`run_rctd.R`** - estimate cell-type fractions per spot from a single-cell RNA-seq reference (RCTD), and use them to describe the niches.

The two steps are independent; either can be run alone.

## 1. Shared niches

Each spot is described by its own expression **and** that of its neighbours (BANKSY), samples are integrated with
Harmony, and the spots are clustered without supervision. The result is one set of niche labels that means the same
thing in every sample, so that niche shares can be compared between patients or conditions.

```
counts + coordinates per sample
   -> QC                         (genes per spot, mitochondrial fraction)
   -> depth equalisation         (optional; binomial thinning to the shallowest sample)
   -> variable genes             (2,000 by default)
   -> BANKSY features            (own expression + neighbour mean + neighbourhood gradient), per sample
   -> equal draw per sample      (so large samples do not dominate)
   -> BANKSY PCA                 (lambda sets the weight of the neighbourhood)
   -> Harmony                    (over sample, or any sample-sheet column)
   -> graph clustering           (several resolutions)
   -> labels for every spot      (nearest niche centroid, with a leave-one-sample-out check)
```

## Requirements

R >= 4.3 with: `Seurat` (v5), `Banksy` (>= 1.0, Bioconductor), `harmony`, `SpatialExperiment`, `presto`, `FNN`,
`mclust`, `data.table`, `Matrix`, `yaml`, `ggplot2`, `patchwork`; for step 2 also `spacexr`. Optional: `hdf5r` (10x `.h5` input), `qs`
(`.qs` Seurat objects), `arrow` (`tissue_positions.parquet`).

```r
install.packages(c("Seurat", "harmony", "FNN", "mclust", "data.table", "yaml", "ggplot2", "patchwork", "hdf5r", "remotes"))
BiocManager::install(c("Banksy", "SpatialExperiment"))
remotes::install_github("immunogenomics/presto")
remotes::install_github("dmcable/spacexr")            # step 2 only
```

Tested with R 4.5.2, Seurat 5.6.0, Banksy 1.6.0, harmony 2.1.0, spacexr 2.2.1.

## Quick start

```bash
git clone https://github.com/hwanglab/banksy-harmony-niches.git && cd banksy-harmony-niches
bash scripts/run_toy_example.sh          # niches: simulates 4 samples with 3 domains, runs everything, checks the answer
bash scripts/run_toy_rctd.sh             # cell types: simulates a reference and known mixtures, checks the fractions
```

Each should end with `PASSED`; `example/results/` and `example_rctd/results/` then hold every output described below.

## Running on your data

1. Write a **sample sheet** (tab-separated):

   | column | meaning |
   |---|---|
   | `sample_id` | unique name (no `|`) |
   | `path` | input file or folder |
   | `format` | `spaceranger`, `10x_h5`, `seurat_qs`, `seurat_rds` or `rds_list` |
   | `coords` | optional: a file with columns `barcode`, `x`, `y` that overrides the coordinates |
   | anything else | kept as metadata (e.g. `batch`, `cohort`); usable for Harmony and for colouring the UMAP |

   | format | `path` points to | coordinates come from |
   |---|---|---|
   | `spaceranger` | Space Ranger `outs/` folder | `outs/spatial/tissue_positions*.csv` (or `.parquet`) |
   | `10x_h5` | a `filtered_feature_bc_matrix.h5` | `spatial/` next to the file |
   | `seurat_qs`, `seurat_rds` | a Seurat object with raw counts | metadata columns named in `input.x_col` / `input.y_col`, else `GetTissueCoordinates()` |
   | `rds_list` | `list(counts = <genes x spots sparse matrix>, coords = <data.frame x, y; row names = spots>)` | the list |

2. Copy `config/example_config.yaml`, set `samples` and `outdir`, and change what you need. Every other key has a default.

3. Run:

   ```bash
   Rscript run_pipeline.R --config my_config.yaml
   # or on a cluster: edit slurm/run_pipeline.sbatch, then
   sbatch slurm/run_pipeline.sbatch my_config.yaml
   ```

   `--samples` and `--outdir` on the command line override the config.

## Outputs

| file | content |
|---|---|
| `clusters.tsv` | the drawn spots: sample, barcode, coordinates, UMAP, cluster at every resolution (`cluster` = main resolution) |
| `all_spots_clusters.tsv.gz` | every spot of every sample with its niche and the correlation to that niche's centroid |
| `cluster_summary.tsv` | size, number of samples in which the niche is present, share of its largest sample, top markers |
| `markers.tsv` | marker genes per niche (`presto`, AUC) |
| `sample_share.tsv`, `sample_share_all_spots.tsv` | niche share per sample (drawn spots / all spots) |
| `stability.tsv` | adjusted Rand index between two independent random draws, per resolution |
| `assignment_check.tsv` | leave-one-sample-out agreement of the nearest-centroid labels with the clustering |
| `sample_qc.tsv` | spots in, spots passing QC, median depth, used or not |
| `umap.png`, `spatial/<sample>.png` | figures |
| `seurat_build.rds`, `config_used.yaml`, `sessionInfo.txt`, `log.txt` | object and provenance |

## 2. Cell-type fractions with RCTD

Given a single-cell RNA-seq reference with cell-type labels, `run_rctd.R` estimates, for every spot, the fraction of each cell type
(RCTD, `spacexr`). It reads the same sample sheet as step 1.

```bash
Rscript run_rctd.R --config my_rctd_config.yaml                    # all samples, one after the other
Rscript run_rctd.R --config my_rctd_config.yaml --sample S01       # one sample (Slurm arrays: slurm/run_rctd_array.sbatch)
```

**Reference** (`reference:` in `config/example_rctd_config.yaml`):

| format | what it is |
|---|---|
| `seurat_rds`, `seurat_qs` | a Seurat object with raw counts; `celltype_col` names the metadata column with the labels |
| `rds_list` | `list(counts = <genes x cells raw counts>, cell_types = <named vector, or data.frame(cell, type)>)` |

Raw integer counts are required. Cell types with fewer than `min_cells_per_type` cells (25) are dropped; each type is capped at `max_cells_per_type` cells.
Characters other than letters, digits, `_ . + -` in type names are replaced by `_`.

**Mode** (`rctd.mode`): `full` allows any number of cell types per spot and is the usual choice for 55 um Visium spots, which hold several cells.
`multi` (up to 4 types) and `doublet` (1-2 types) are meant for finer platforms.

**Outputs** (in `outdir`):

| file | content |
|---|---|
| `<sample>_weights.tsv.gz` | spots x cell types; each row sums to 1 |
| `<sample>_group_weights.tsv.gz` | the same, summed into the broader classes given under `groups:` (optional) |
| `<sample>_doublet_calls.tsv.gz` | `doublet` mode only: RCTD's per-spot classification |
| `mean_weight_per_sample.tsv` | average fraction of each cell type per sample |
| `dominant_type_share_per_sample.tsv` | share of spots in which each type has the largest weight |
| `rctd_log.txt`, `rctd_config_used.yaml`, `rctd_sessionInfo.txt` | provenance |

**Using the fractions to describe the niches.** Set `annotation.weights_pattern` in the step-1 config to
`"<rctd outdir>/{sample_id}_weights.tsv.gz"`; step 1 then also writes `cluster_celltype_weights.tsv` (mean cell-type weight per niche and the leading type).

**Cautions**

- The fractions are only as good as the reference: a cell type missing from the reference is assigned to its nearest neighbour, and a reference from untreated
  tissue applied to treated tissue cannot be checked from the data alone.
- Weights are fractions of RNA attributed to a type, not cell counts.
- Spots below `rctd.umi_min` are not deconvolved and are absent from the output.
- Do not define groups on these labels and then test the same labels; use them to describe, and confirm with a direct measurement.

## Choices that change the niche answer

- **`banksy.lambda`** - 0.2 keeps clusters close to cell-type composition; 0.8 gives broader tissue domains.
- **`banksy.k_geom`** - the first ring of neighbours: 6 on a Visium hexagonal grid, 8 on a square grid (binned Visium HD, pseudo-spots). With `use_agf: true`
  the borders between domains can come out as clusters of their own.
- **`harmony.vars`** - integrate over `sample_id` when the tumour compartment is patient-specific; integrating only over a coarser variable (cohort, batch)
  tends to leave one cluster per patient. Check `largest_sample_share` in `cluster_summary.tsv`: a value near 1 means a single-sample cluster.
- **`subsample.n_per_sample`** - an equal draw per sample; niche shares from the draw are unbiased, and labels are then carried to all spots.
- **`depth.equalise`** - turn on when median depth differs several-fold between samples; otherwise shallow samples cluster by depth.
- **`cluster.main_resolution`** - read `stability.tsv` first. A resolution whose clusters are not reproduced by a second draw should not be interpreted.

## Reading the result with care

- Niche labels are derived from the data they are later compared on. Use them to describe tissue; treat a between-group difference in niche share as a lead to be
  confirmed with a direct measurement.
- The unit for any statistical comparison is the patient (or sample), not the spot.
- A spot holds several cells; a niche is a neighbourhood, not a cell type.
- `all_spots_clusters.tsv.gz` is a nearest-centroid transfer; `assignment_check.tsv` tells you how well it reproduces the clustering in a held-out sample.

## Citation

If you use this pipeline, please cite the methods it is built on:

- Singhal V, et al. BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. *Nat Genet* 2024.
- Korsunsky I, et al. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 2019.
- Cable DM, et al. Robust decomposition of cell type mixtures in spatial transcriptomics. *Nat Biotechnol* 2022 (RCTD).
- Hao Y, et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. *Nat Biotechnol* 2024 (Seurat v5).

## License

MIT - see `LICENSE`.
