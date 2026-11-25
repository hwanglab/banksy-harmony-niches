# Tables and figures.

write_cluster_tables <- function(so, main_col, out, cfg) {
  so$cluster <- paste0("C", so[[main_col]][, 1])
  mk <- data.table::as.data.table(presto::wilcoxauc(so, group_by = "cluster"))
  data.table::fwrite(mk[auc >= cfg$report$marker_min_auc][order(group, -auc)], file.path(out, "markers.tsv"), sep = "\t")
  cls <- names(sort(table(so$cluster), decreasing = TRUE))
  summ <- data.table::rbindlist(lapply(cls, function(k) {
    i <- so$cluster == k
    data.table::data.table(cluster = k, spots = sum(i), samples_with_min_spots = sum(table(so$sample_id[i]) >= cfg$report$min_spots_per_sample),
                           largest_sample_share = round(max(table(so$sample_id[i])) / sum(i), 2),
                           top_markers = paste(mk[group == k][order(-auc)][seq_len(min(cfg$report$n_markers, .N))]$feature, collapse = ", "))
  }))
  data.table::fwrite(summ, file.path(out, "cluster_summary.tsv"), sep = "\t")
  tb <- data.table::as.data.table(table(sample_id = so$sample_id, cluster = so$cluster))
  tb[, share := N / sum(N), by = sample_id]
  data.table::fwrite(data.table::dcast(tb, sample_id ~ cluster, value.var = "share"), file.path(out, "sample_share.tsv"), sep = "\t")
  list(so = so, summary = summ)
}

.theme <- function() ggplot2::theme_classic(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(size = 16))

plot_umap <- function(so, out, colour_by) {
  d <- data.table::data.table(Seurat::Embeddings(so, "umap"), cluster = so$cluster, sample_id = so$sample_id)
  for (v in colour_by) if (v %in% names(so[[]])) d[[v]] <- as.character(so[[v]][, 1])
  names(d)[1:2] <- c("umap_1", "umap_2")
  panels <- lapply(unique(c("cluster", colour_by)), function(v) {
    if (!v %in% names(d)) return(NULL)
    ggplot2::ggplot(d, ggplot2::aes(umap_1, umap_2, colour = .data[[v]])) + ggplot2::geom_point(size = 0.2) + .theme() +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4), title = v)) + ggplot2::ggtitle(paste("coloured by", v))
  })
  panels <- Filter(Negate(is.null), panels)
  grDevices::png(file.path(out, "umap.png"), width = 750 * length(panels), height = 680, res = 110)
  print(patchwork::wrap_plots(panels, nrow = 1))
  grDevices::dev.off()
}

plot_spatial <- function(tab, out, max_samples) {
  dir.create(file.path(out, "spatial"), showWarnings = FALSE)
  lv <- names(sort(table(tab$cluster), decreasing = TRUE))
  pal <- stats::setNames(grDevices::hcl.colors(length(lv), "Dynamic"), lv)
  for (s in utils::head(unique(tab$sample_id), max_samples)) {
    d <- tab[sample_id == s]
    p <- ggplot2::ggplot(d, ggplot2::aes(x, -y, colour = factor(cluster, levels = lv))) + ggplot2::geom_point(size = 0.4) + ggplot2::coord_equal() +
      ggplot2::scale_colour_manual(values = pal, drop = FALSE, name = "cluster") + .theme() + ggplot2::theme(axis.title = ggplot2::element_blank()) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4))) + ggplot2::ggtitle(s)
    ggplot2::ggsave(file.path(out, "spatial", paste0(s, ".png")), p, width = 8, height = 7, dpi = 130)
  }
}
