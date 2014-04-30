## This file shows an implementation of copy number using HMMcopy. 
## This is run on our chosen dataset
## tumor file : tumorA.chr4.bam
## normal file : normalA.chr4.bam

## outside R - generate readCounts file 
## bin/readCounter tumorA.chr4.bam > chr4_tum_reads.wig
## bin/readCounter normalA.chr4.bam > chr4_norm_reads.wig


tum_readfile <-"chr4_tum_reads.wig"
nor_readfile <-"chr4_norm_reads.wig"
mapfile <-"data/map_hg19.wig"
gcfile <-"data/gc_hg19.wig"

tum_uncorrected_reads <- wigsToRangedData(tum_readfile, gcfile, mapfile)
norm_uncorrected_reads <- wigsToRangedData(nor_readfile, gcfile, mapfile)

##correct read counts

tum_corrected_copy <- correctReadcount(tum_uncorrected_reads)
norm_corrected_copy <- correctReadcount(norm_uncorrected_reads)

# Normalizing Tumour by Normal
tum_corrected_copy$copy <- tum_corrected_copy$copy - norm_corrected_copy$copy

# Export to SEG format for CNAseq segmentation
rangedDataToSeg(tum_corrected_copy, file = "paul_tum_corrected_copy.seg")

# Segmenting
## use default segmentation 
seg_copy_def <- HMMsegment(tum_corrected_copy)

##get parametrs
realparam <- HMMsegment(tum_corrected_copy, getparam = TRUE) # retrieve converged parameters via EM

## Adjust parameters - case1
param1 <- realparam
param1$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
param1$m <- param1$mu
segmented_copy_case1 <- HMMsegment(tum_corrected_copy, param1) # perform segmentation via Viterbi

##adjust parameters - case2 ## to decrease no of segments/ 
param2 <- realparam
param2$strength <- 1e30
param2$e <- 0.99999999999999
segmented_copy_case2 <- HMMsegment(tum_corrected_copy, param2)

##adjust parameters - case2 ## to increase no of segments/ 
param3 <- realparam
param3$strength <- 0.1
param3$e <- 0.1
segmented_copy_case3 <- HMMsegment(tum_corrected_copy, param3)

# visualization

plotBias(tum_corrected_copy)
plotCorrection(tum_corrected_copy)
plotSegments(tum_corrected_copy, segmented_copy)

png("paul-plotBias.png");plotBias(tum_corrected_copy);dev.off()

png("paul-plotCorrection.png");plotCorrection(tum_corrected_copy); dev.off()

png("paul-plotCorrection.png");plotCorrection(tum_corrected_copy); dev.off()
