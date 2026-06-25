library(SummarizedExperiment)

set.seed(42)

n_rows        <- 100L
n_cols        <- 20L
txps_per_gene <- 10L
n_genes       <- n_rows / txps_per_gene

row_names <- paste0("txp",    seq_len(n_rows))
col_names <- paste0("sample", seq_len(n_cols))

counts    <- matrix(rpois(n_rows * n_cols, lambda = 100), nrow = n_rows, ncol = n_cols,
                    dimnames = list(row_names, col_names))
abundance <- matrix(exp(rnorm(n_rows * n_cols)), nrow = n_rows, ncol = n_cols,
                    dimnames = list(row_names, col_names))
length    <- matrix(round(runif(n_rows * n_cols, min = 1e3, max = 1e4)), nrow = n_rows, ncol = n_cols,
                    dimnames = list(row_names, col_names))

col_data <- DataFrame(
  condition = rep(c("ctrl", "treat"), each = n_cols / 2),
  batch     = rep(c("A", "B"), times = n_cols / 2),
  size_factor = runif(n_cols, 0.5, 2.0),
  row.names = col_names
)

row_data <- DataFrame(
  gene       = rep(paste0("gene", seq_len(n_genes)), each = txps_per_gene),
  gene_type  = sample(c("protein_coding", "lncRNA", "pseudogene"), n_rows, replace = TRUE),
  chrom      = rep(sample(paste0("chr", c(1:22, "X", "Y")), n_genes, replace = TRUE), each = txps_per_gene),
  gc_content = runif(n_rows, 0.3, 0.7),
  row.names  = row_names
)

se <- SummarizedExperiment(
  assays   = list(counts = counts, abundance = abundance, length = length),
  colData  = col_data,
  rowData  = row_data
)

se

colData(se)

rowData(se)

random_hclust <- function(n) {
  merge_mat <- matrix(0L, n - 1L, 2L)
  available <- -(seq_len(n))
  for (i in seq_len(n - 1L)) {
    pick <- sample(length(available), 2L)
    merge_mat[i, ] <- sort(available[pick])
    available <- c(available[-pick], i)
  }
  leaf_order <- local({
    get_leaves <- function(node) {
      if (node < 0L) return(-node)
      c(get_leaves(merge_mat[node, 1L]), get_leaves(merge_mat[node, 2L]))
    }
    get_leaves(n - 1L)
  })
  structure(
    list(
      merge  = merge_mat,
      height = as.double(seq_len(n - 1L)),
      order  = leaf_order,
      labels = NULL,
      method = "random",
      dist.method = NULL
    ),
    class = "hclust"
  )
}

# a little plyxp to make a new assay

library(plyxp)
xp <- se |> new_plyxp()
xp <- xp |>
  mutate(rel_abundance = 1e6 * t(t(abundance) / colSums(abundance)))
colSums(assay(xp, "rel_abundance"))

library(treeClustDemo)
tree_xp <- xp |>
  group_by(rows(gene)) |>
  mutate(rows(
    nodes = random_hclust(n()) |> tree_vctr()
  )) |>
  ungroup()

# Drop assays not needed downstream
tree_xp <- tree_xp |>
  select(
    counts, rel_abundance,
    rows(everything()), 
    cols(everything())
)

# Aggregate assays across inner nodes; carry gene and chrom from first child
inner_xp <- group_by_nodes(tree_xp, rows(inner_nodes(nodes))) |>
  summarize(
    across(everything(), colSums),
    rows(gene = gene[1], chrom = chrom[1])
  )

# rbind leaves and inner nodes after aligning to shared rowData columns
common_rdcols <- intersect(
  names(rowData(se(tree_xp))),
  names(rowData(se(inner_xp)))
)

full_xp <- rbind(
  new_plyxp(se(tree_xp)) |>
    select(counts, rel_abundance, rows(all_of(common_rdcols)), cols(everything())) |>
    ungroup() |>
    se(),
  se(inner_xp)
) |> new_plyxp() |> ungroup()

# ---- Benchmark ----------------------------------------------------------------

make_se_with_tree <- function(n_genes, txps_per_gene = 10L, n_cols = 20L) {
  n_rows    <- n_genes * txps_per_gene
  row_names <- paste0("txp",    seq_len(n_rows))
  col_names <- paste0("sample", seq_len(n_cols))

  counts    <- matrix(rpois(n_rows * n_cols, lambda = 100), nrow = n_rows, ncol = n_cols,
                      dimnames = list(row_names, col_names))
  abundance <- matrix(exp(rnorm(n_rows * n_cols)), nrow = n_rows, ncol = n_cols,
                      dimnames = list(row_names, col_names))

  se <- SummarizedExperiment(
    assays  = list(counts = counts, abundance = abundance),
    colData = DataFrame(
      condition   = rep(c("ctrl", "treat"), each = n_cols / 2L),
      batch       = rep(c("A", "B"),        times = n_cols / 2L),
      size_factor = runif(n_cols, 0.5, 2.0),
      row.names   = col_names
    ),
    rowData = DataFrame(
      gene  = rep(paste0("gene", seq_len(n_genes)), each = txps_per_gene),
      chrom = rep(sample(paste0("chr", c(1:22, "X", "Y")), n_genes, replace = TRUE),
                  each = txps_per_gene),
      row.names = row_names
    )
  )

  new_plyxp(se) |>
    mutate(rel_abundance = 1e6 * t(t(abundance) / colSums(abundance))) |>
    select(counts, rel_abundance, rows(everything()), cols(everything())) |>
    group_by(rows(gene)) |>
    mutate(rows(nodes = random_hclust(n()) |> tree_vctr())) |>
    ungroup()
}

bench_params <- expand.grid(
  n_genes       = as.integer(500,1000,1500,2000),
  txps_per_gene = c(10L, 20L)
)

# Build all inputs outside the timed region
bench_inputs <- Map(make_se_with_tree, bench_params$n_genes, bench_params$txps_per_gene)

library(microbenchmark)

bench_results <- lapply(seq_len(nrow(bench_params)), function(i) {
  xp <- bench_inputs[[i]]
  mb <- microbenchmark(
    group_by_nodes(xp, rows(inner_nodes(nodes))) |>
      summarize(
        across(everything(), colSums),
        rows(gene = gene[1], chrom = chrom[1])
      ),
    times = 5L,
    unit  = "s"
  )
  cbind(bench_params[i, ], summary(mb)[, c("min", "lq", "median", "uq", "max")])
})

bench_df <- do.call(rbind, bench_results)
bench_df$txps_per_gene <- factor(bench_df$txps_per_gene)

library(ggplot2)
ggplot(bench_df, aes(x = n_genes, y = median, colour = txps_per_gene)) +
  geom_line() +
  geom_point() +
  labs(
    x      = "Number of genes",
    y      = "Time (s)",
    colour = "Txps per gene",
    title  = "group_by_nodes |> summarize timing"
  ) +
  theme_minimal()

