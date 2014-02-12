library(GenomicAlignments)
library(RUnit)
    
cnv_countOverlaps <- 
    function(bamfilename, chrom, binSize=10000)
{
    bf <- BamFile(bamfilename)
    aln<- readGAlignments(bf)
    sub <- resize(granges(aln), width=1)
    tiles <- tileGenome(seqlengths(bf), tilewidth=binSize, 
                        cut.last.tile.in.chrom=TRUE)
    tiles$cnt <- countOverlaps(tiles, sub)
    tiles[seqnames(tiles)==chrom]
}


test_cnv_countOverlaps <- 
    function()
{
    binSize=10000  
    tumor.bins <- cnv_countOverlaps("tumorA.chr4.bam", 
                                    chrom="chr4", 
                                    binSize=binSize)
    
    checkEquals(min(tumor.bins$cnt), 0)
    checkEqualsNumeric(mean(tumor.bins$cnt), 56.66311, tolerance=10e-2)
    checkEqualsNumeric(max(tumor.bins$cnt),  8638,   tolerance=10e-2)
}