## change seqnameStyle from UCSC to NCBI - done in unix.
## sed s/chr//g paul_tumAlleleCounts.tsv > paul_tumAlleleCounts_ncbi.tsv
## sed s/chr//g chr4_tumor_reads.wig > chr4_tumor_reads_ncbi.wig
## sed s/chr//g chr4_normal_reads.wig > chr4_normal_reads_ncbi.wig

library(TitanCNA)

## load the files
id <- "test"
infile <-"paul_tumAlleleCounts_ncbi.tsv" 
tumWig <- "chr4_tumor_reads_ncbi.wig" 
normWig <- "chr4_normal_reads_ncbi.wig"
gcWig <- "GRCh37-lite.gc.ws_1000.wig"
mapWig <- "GRCh37-lite.map.ws_1000.wig"

## load the tumor allele read counts
data <- loadAlleleCountsFromFile(infile)

## correct gc and mappability bias
cnData <- correctReadDepth(tumWig,normWig,gcWig,mapWig)
