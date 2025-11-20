## code to prepare `glinos` dataset goes here

google_url <- function(id) {
  sprintf("https://drive.google.com/uc?export=download&id=%s", id)
}
url <- google_url("1ATHgwFlIMmI651C5TYtiBxm2atTvbBYn")
filename <- "quantification_flair_filter.counts.txt.gz"

path_counts <- tempfile(fileext = filename)
download.file(url, path_counts)

# GTF
url <- google_url("1S0TRXoXsfnXwJ0Dh5bLqVmyvFAAALFMU")
filename <- "flair_filter_transcripts.gtf.gz"

path_gtf <- tempfile(fileext = filename)
download.file(url, path_gtf)

# metadata
url <- google_url("1LpYPXXhMtUV3vMG7SMLV69-JxZfc59o2")
filename <- "glinos_metadata.csv"

#' @describeIn local_files The file location of the metadata. see
#' [download][.__module__.]
#' @export
path_metadata <- tempfile(fileext = filename)
download.file(url, path_metadata)
download <- list(
  path_gtf = path_gtf,
  path_metadata = path_metadata,
  path_counts = path_counts
)
txdb <- txdbmaker::makeTxDbFromGFF(download$path_gtf)
exons <- ensembldb::exonsBy(txdb, by = "tx")
txps <- AnnotationDbi::select(
  txdb, AnnotationDbi::keys(txdb, "TXID"),
  c("TXNAME", "GENEID"), "TXID"
) |>
  tibble::as_tibble() |>
  dplyr::mutate(TXID = as.character(TXID))

# better names for exons:
stopifnot(all.equal(names(exons), txps$TXID))

names(exons) <- dplyr::case_when(
  stringr::str_detect(txps$TXNAME, "ENST") ~ txps$TXNAME,
  TRUE ~ paste0(txps$GENEID, "-", stringr::str_sub(txps$TXNAME, 1, 8))
)

counts <- readr::read_delim(download$path_counts)
meta <- readr::read_delim(download$path_metadata)
colnames(meta)[10] <- "total_reads"
meta <- meta |>
  dplyr::mutate(condition = stringr::str_extract(sample_id, "(ctrl|exp)")) |>
  dplyr::select(
    sample_id, sample_name, condition,
    dplyr::contains("read") | dplyr::contains("3_prime")
  ) |>
  dplyr::arrange(condition, sample_name)

library(plyxp)

rnames <- counts$transcript
counts <- as.matrix(counts[, meta$sample_id]) # subset to the knockdown
rownames(counts) <- rnames
glinos <- SummarizedExperiment(
  list(counts = counts),
  colData = meta
) |>
  new_plyxp() |>
  filter(
    rows(
      rowSums(.assays_asis$counts >= 10) >= 6
    )
  ) |>
  mutate(
    cols(sums = colSums(.assays_asis$counts)),
    abundance = counts / .cols$sums * 1e6,
    length = rep(1000, n()),
    cols(
      condition = dplyr::case_when(
        condition == "ctrl" ~ "WT",
        condition == "exp" ~ "KD"
      ),
      # note that protein knockdown will be used as reference
      # with this choice of factor levels
      condition = factor(condition, c("KD", "WT"))
    ),
    rows(
      .features = stringr::str_replace(.features, "_ENSG\\d+\\.\\d+$", "")
    )
  ) |>
  se()

slice_indx <- match(rownames(glinos), txps$TXNAME, nomatch = 0L)
txps <- txps[slice_indx, ]
# txps <- txdb$txps |>
# filter(tx_nms)

glinos <- glinos[txps$TXNAME, ]
stopifnot(all.equal(rownames(glinos), txps$TXNAME))
exons <- exons[slice_indx]
rownames(glinos) <- names(exons)
tx2gene <- S4Vectors::DataFrame(
  isoform_id = rownames(glinos),
  gene_id = txps$GENEID,
  txps_full_name = txps$TXNAME
)
rowData(glinos) <- tx2gene
rowData(glinos)$exons <- exons
usethis::use_data(glinos, overwrite = TRUE)
