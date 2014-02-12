library(GenomicAlignments)
library(RUnit)

countBinOverlaps <-
  function(features, reads, ...)
  {
    reads <- resize(granges(reads), width=1)
    countOverlaps(features, reads)
  }

summarizeBins <-
  function(file, seqnames, tilewidth=10000)
  {
    stopifnot(is(file, "character") && length(file) > 0)
    stopifnot(is(seqnames, "character") && length(seqnames) > 0)
    
    seqlengths <- seqlengths(BamFile(file[[1]]))
    tiles <- tileGenome(seqlengths[names(seqlengths) %in% seqnames],
                        tilewidth=tilewidth,
                        cut.last.tile.in.chrom=TRUE)
    
    summarizeOverlaps(tiles, fls, countBinOverlaps)
  }

fls <- dir("~/benchmark/copynumber/", pattern="bam$", full=TRUE)
counts <- summarizeBins(fls, "chr4")


