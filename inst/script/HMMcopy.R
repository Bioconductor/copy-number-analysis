## This file shows an implementation of copy number using HMMcopy. 
## This is run on our chosen dataset
## tumor file : tumorA.chr4.bam
## normal file : normalA.chr4.bam

## outside R - generate readCounts file 
## bin/readCounter tumorA.chr4.bam > chr4_tum_reads.wig
## bin/readCounter normalA.chr4.bam > chr4_norm_reads.wig

## currently the wig files for gc content and mappability have NCBI style
## of seqnames/ chromsomes. so convert tumor and normal wig files to same. 
## sed s/=chr/=/ chr4_tumor_reads.wig > chr4_tumor_reads_ncbi.wig
## sed s/=chr/=/ chr4_tumor_reads.wig > chr4_tumor_reads_ncbi.wig

tum_readfile <-"chr4_tumor_reads_ncbi.wig"
nor_readfile <-"chr4_normal_reads_ncbi.wig"

## Note - these files are distributed along with TitanCNA
## the files distributed along with HMMcopy had inconsistent seqname style
mapfile <-"GRCh37-lite.map.ws_1000.wig"
gcfile <-"GRCh37-lite.gc.ws_1000.wig"

## create a RangedData object.
tum_uncorrected_reads <- wigsToRangedData(tum_readfile, gcfile, mapfile)
norm_uncorrected_reads <- wigsToRangedData(nor_readfile, gcfile, mapfile)

## subset to have reads only from chr4
tum_uc_reads<- tum_uncorrected_reads["4"]
norm_uc_reads<- norm_uncorrected_reads["4"]

##correct read counts
tum_corrected_copy <- correctReadcount(tum_uc_reads)
norm_corrected_copy <- correctReadcount(norm_uc_reads)

## Normalizing Tumour by Normal
tum_corrected_copy$copy <- tum_corrected_copy$copy - norm_corrected_copy$copy

## Export to SEG format for CNAseq segmentation
rangedDataToSeg(tum_corrected_copy, file = "paul_tum_corrected_copy.seg")

## Segmenting
## use default segmentation 
seg_copy_def <- HMMsegment(tum_corrected_copy)

## get parametrs
realparam <- HMMsegment(tum_corrected_copy, getparam = TRUE) # retrieve converged parameters via EM

## Adjust parameters - case1
param1 <- realparam
param1$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
param1$m <- param1$mu
segmented_copy_case1 <- HMMsegment(tum_corrected_copy, param1) # perform segmentation via Viterbi

## adjust parameters - case2 ## to decrease no of segments/ 
param2 <- realparam
param2$strength <- 1e30
param2$e <- 0.99999999999999
segmented_copy_case2 <- HMMsegment(tum_corrected_copy, param2)

## adjust parameters - case2 ## to increase no of segments/ 
param3 <- realparam
param3$strength <- 0.1
param3$e <- 0.1
segmented_copy_case3 <- HMMsegment(tum_corrected_copy, param3)

## visualization

plotBias(tum_corrected_copy)
plotCorrection(tum_corrected_copy)
plotSegments(tum_corrected_copy, seg_copy_def)
plotSegments(tum_corrected_copy, segmented_copy_case1)
plotSegments(tum_corrected_copy, segmented_copy_case2)
plotSegments(tum_corrected_copy, segmented_copy_case3)



