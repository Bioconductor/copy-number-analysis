## This script and helper functions implements an R-only version of
## the work flow from http://tiger.dbs.nus.edu.sg/cnv-seq/. It uses
## GenomicRanges utilities to perform read counts across bins, and the
## 'cnv' package available at the URL above for additional
## analysis. No intermediate files are generates.

## user arguments; provide full paths
root <- "~/benchmark/copynumber"
files <- file.path(root, c("tumorA.chr4.bam", "normalA.chr4.bam"))
names(files) <- c("test", "ref")

log2 <- .6
annotate <- TRUE

## script
source("../../R/RCNV_seq-helper.R")
suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(cnv)
})

window_size <- windowSize(files, log2=log2)
tiles <- tileGenomeOverlap(files, window_size)
hits <- summarizeOverlaps(tiles, files, binCounter)

chr4 <- subset(hits, seqnames %in% "chr4")
cnv <- cnv.cal(as.countsfile(chr4), log2=log2, annotate=annotate)
plot.cnv(cnv)
