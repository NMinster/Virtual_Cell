#!/usr/bin/env Rscript
# Convert a Seurat or SingleCellExperiment .rds (the usual GEO supplementary
# object) into a 10x-style mtx folder that `ovary-fm qc` reads:
#
#     Rscript R/rds_to_mtx.R object.rds data/raw/<dataset>/ [assay] [layer]
#
# Writes matrix.mtx.gz (genes x cells, raw counts), features.tsv.gz,
# barcodes.tsv.gz and meta.csv (the object's cell metadata, one row per
# barcode) into the folder. No SeuratDisk, no h5ad writer: Matrix and
# whichever of Seurat / SingleCellExperiment the object needs.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript R/rds_to_mtx.R object.rds outdir [assay] [layer]")
rds <- args[1]; outdir <- args[2]
assay <- if (length(args) >= 3) args[3] else NULL
layer <- if (length(args) >= 4) args[4] else "counts"
suppressPackageStartupMessages(library(Matrix))
obj <- readRDS(rds)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (inherits(obj, "Seurat")) {
  suppressPackageStartupMessages(library(Seurat))
  if (is.null(assay)) assay <- DefaultAssay(obj)
  counts <- tryCatch(GetAssayData(obj, assay = assay, layer = layer),
                     error = function(e) GetAssayData(obj, assay = assay, slot = layer))
  meta <- obj@meta.data
} else if (inherits(obj, "SingleCellExperiment")) {
  suppressPackageStartupMessages(library(SingleCellExperiment))
  counts <- SummarizedExperiment::assay(obj, layer)
  meta <- as.data.frame(SummarizedExperiment::colData(obj))
} else if (inherits(obj, "dgCMatrix") || is.matrix(obj)) {
  counts <- as(obj, "dgCMatrix"); meta <- data.frame(row.names = colnames(counts))
} else {
  stop("not a Seurat, SingleCellExperiment or matrix object: ", class(obj)[1])
}
counts <- as(counts, "dgCMatrix")
if (any(counts@x != round(counts@x))) warning("the matrix is not integer valued: is this the counts layer?")
writeMM(counts, file.path(outdir, "matrix.mtx"))
system2("gzip", c("-f", file.path(outdir, "matrix.mtx")))
feat <- data.frame(id = rownames(counts), name = rownames(counts), type = "Gene Expression")
con <- gzfile(file.path(outdir, "features.tsv.gz"), "w"); write.table(feat, con, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE); close(con)
con <- gzfile(file.path(outdir, "barcodes.tsv.gz"), "w"); writeLines(colnames(counts), con); close(con)
write.csv(meta, file.path(outdir, "meta.csv"))
cat(sprintf("wrote %d genes x %d cells to %s (meta.csv: %d columns)\n", nrow(counts), ncol(counts), outdir, ncol(meta)))
