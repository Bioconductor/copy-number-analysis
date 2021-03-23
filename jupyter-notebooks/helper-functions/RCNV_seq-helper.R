## samtools view -F 4 tumorA.chr4.bam |\
##     perl -lane 'print "$F[2]\t$F[3]"' >tumor.hits
## samtools view -F 4 normalA.chr4.bam |\
##     perl -lane 'print "$F[2]\t$F[3]"' >normal.hits
## perl cnv-seq/cnv-seq.pl --test tumor.hits --ref normal.hits \
##     --genome human --Rexe "~/bin/R-devel/bin/R"

genomeSize <- function(files)
    sum(as.numeric(seqlengths(BamFile(files[[1]]))))

windowSize <-
    function(bam_files, pvalue=0.001, log2=0.6, bigger=1.5,
             genome_size, param)
{
    if (missing(genome_size))
        genome_size=genomeSize(bam_files)
    if (missing(param))
        param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE))

    total <- sapply(bam_files, function(...) {
        countBam(...)$records
    }, param=param)

    bt <- qnorm(1 - pvalue / 2)
    st <- qnorm(pvalue / 2)
    log2 <- abs(log2)
    brp <- 2^log2
    srp <- 1 / (2^log2)
    
    bw <- (total[["test"]] * brp^2 + total[["ref"]]) * genome_size * bt^2 /
        ((1-brp)^2 * total[["test"]] * total[["ref"]])
    sw <- (total[["test"]] * srp^2 + total[["ref"]]) * genome_size * st^2 /
        ((1-srp)^2 * total[["test"]] * total[["ref"]])

    window_size = floor(max(bw, sw) * bigger)
}

tileGenomeOverlap <- function(file, tilewidth) {
    ## overlapping tiles
    lengths <- seqlengths(BamFile(file[[1]]))
    tile0 <- tileGenome(lengths, tilewidth=tilewidth,
                        cut.last.tile.in.chrom=TRUE)
    tile1 <- tile0[width(tile0) >= tilewidth]
    tile1 <- shift(tile1[-cumsum(runLength(seqnames(tile1)))], tilewidth / 2)
    sort(c(tile0, tile1))
}

binCounter <- function(features, reads, ignore.strand, ...) {
    countOverlaps(features, resize(granges(reads), 1),
                  ignore.strand=ignore.strand) 
}

# as.countsfile <- function(hits, file=tempfile()) {
#     df <- with(rowData(hits), {
#         cbind(data.frame(chromosome=as.character(seqnames),
#                          start=start, end=end),
#               assay(hits))
#     })
#     write.table(df, file, quote=FALSE, row.names=FALSE, sep="\t")
#     file
# }


################### Modified as.countsfile function #################
as.countsfile <- function(hits, fileNameAndLocation) {
    df <- with(rowData(hits), {
        chrome = hits@rowRanges@seqnames
        start = hits@rowRanges@ranges
        end = hits@rowRanges@ranges@start + hits@rowRanges@ranges@width - 1
        cbind(data.frame(chromosome=chrome,
                         start=start,
                         end=end),
                         assay(hits))
    })

    df = subset(df, select = -c(end))
    names(df)[names(df) == "start.start"] <- "start"
    names(df)[names(df) == "start.end"] <- "end"
    names(df)[names(df) == "start.width"] <- "width"

    write.table(df, file = fileNameAndLocation, quote=FALSE, row.names=FALSE, sep="\t")
    return(fileNameAndLocation)
}

################### Added function toGRanges #################
toGRanges <- function(bed_dataframe){
    GRanges(seqnames = bed_dataframe$V1, 
            ranges = IRanges(start = bed_dataframe$V2,
                             end = bed_dataframe$V3))
}