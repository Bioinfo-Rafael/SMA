#!/usr/bin/env Rscript
# rds_to_h5ad_bridge.R <input.rds> <out_dir>
#
# Reads an RDS object, identifies its class, and writes a CellRanger-style
# MTX triplet (matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz) plus meta.csv
# that the Python side (src/io_rds_bridge.py) assembles into an AnnData.
#
# Supported:
#   * Seurat               -> counts slot + meta.data + feature names
#   * SingleCellExperiment -> counts assay + colData + rowData
# Anything else: write class_report.txt and stop() (non-zero exit).

suppressWarnings(suppressMessages({
  ok_matrix <- requireNamespace("Matrix", quietly = TRUE)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: rds_to_h5ad_bridge.R <input.rds> <out_dir>")
}
rds_path <- args[[1]]
out_dir  <- args[[2]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

report_path <- file.path(out_dir, "class_report.txt")
write_report <- function(lines) writeLines(lines, report_path)

if (!ok_matrix) {
  write_report(c("ERROR: R package 'Matrix' is required but not installed."))
  stop("Matrix package not available")
}

message("reading: ", rds_path)
obj <- readRDS(rds_path)
cls <- class(obj)
slot_names <- tryCatch(slotNames(obj), error = function(e) character(0))
nm <- tryCatch(names(obj), error = function(e) character(0))

write_report(c(
  paste0("rds: ", rds_path),
  paste0("class: ", paste(cls, collapse = ", ")),
  paste0("slotNames: ", paste(slot_names, collapse = ", ")),
  paste0("names: ", paste(utils::head(nm, 50), collapse = ", "))
))

write_triplet <- function(counts, features, barcodes, meta) {
  counts <- methods::as(counts, "CsparseMatrix")          # genes x cells
  mtx_path <- file.path(out_dir, "matrix.mtx")
  Matrix::writeMM(counts, mtx_path)
  R.utils_ok <- requireNamespace("R.utils", quietly = TRUE)
  # gzip if possible (Python side reads either .mtx or .mtx.gz)
  if (R.utils_ok) {
    R.utils::gzip(mtx_path, overwrite = TRUE, remove = TRUE)
  }
  writeLines(as.character(barcodes), file.path(out_dir, "barcodes.tsv"))
  feat_df <- data.frame(id = features, symbol = features, type = "Gene Expression")
  utils::write.table(feat_df, file.path(out_dir, "features.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  if (!is.null(meta)) {
    utils::write.csv(as.data.frame(meta), file.path(out_dir, "meta.csv"))
  }
  message("wrote MTX triplet: ", nrow(counts), " genes x ", ncol(counts), " cells")
}

is_seurat <- inherits(obj, "Seurat")
is_sce <- inherits(obj, "SingleCellExperiment")

if (is_seurat) {
  message("detected Seurat object")
  if (!requireNamespace("SeuratObject", quietly = TRUE) &&
      !requireNamespace("Seurat", quietly = TRUE)) {
    write_report(c("ERROR: Seurat/SeuratObject not installed to read counts."))
    stop("Seurat not available")
  }
  getdata <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
    SeuratObject::GetAssayData
  } else {
    get("GetAssayData", envir = asNamespace("Seurat"))
  }
  counts <- tryCatch(getdata(obj, slot = "counts"),
                     error = function(e) getdata(obj, layer = "counts"))
  if (is.null(counts) || nrow(counts) == 0) {
    write_report(c("ERROR: empty counts slot in Seurat object"))
    stop("empty counts")
  }
  meta <- obj@meta.data
  write_triplet(counts, rownames(counts), colnames(counts), meta)

} else if (is_sce) {
  message("detected SingleCellExperiment object")
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    write_report(c("ERROR: SingleCellExperiment not installed."))
    stop("SingleCellExperiment not available")
  }
  assays_avail <- SummarizedExperiment::assayNames(obj)
  assay_name <- if ("counts" %in% assays_avail) "counts" else assays_avail[[1]]
  counts <- SummarizedExperiment::assay(obj, assay_name)
  meta <- as.data.frame(SummarizedExperiment::colData(obj))
  write_triplet(counts, rownames(obj), colnames(obj), meta)

} else {
  write_report(c(
    paste0("UNSUPPORTED class: ", paste(cls, collapse = ", ")),
    paste0("slotNames: ", paste(slot_names, collapse = ", ")),
    paste0("names: ", paste(utils::head(nm, 50), collapse = ", "))
  ))
  stop(paste0("Unsupported RDS class: ", paste(cls, collapse = ", "),
              " (see ", report_path, ")"))
}

message("done")
