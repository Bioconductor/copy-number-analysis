library(cn.mops)
library(RUnit)

tumor_gr <- getReadCountsFromBAM("tumorA.chr4.bam",
                                 refSeqName="chr4",WL=10000,mode="unpaired")
normal_gr <- getReadCountsFromBAM("normalA.chr4.bam",
                                  refSeqName="chr4",WL=10000,mode="unpaired")

ref_analysis_norm0 <- referencecn.mops(tumor_gr, normal_gr,norm=0) 

ref_analysis_norm0 <- referencecn.mops(tumor_gr, normal_gr,norm=0) 
resCNMOPS <- calcIntegerCopyNumbers(ref_analysis_norm0)
resCNMOPS <-cn.mops:::.replaceNames(resCNMOPS, "tumorA.chr4.bam","tumor")

png("cn.mops-segplot.png")
segplot(resCNMOPS)
dev.off()
##Here the x-axis represents the genomic position and the y-axis represents the 
##log ratio of read counts and copy number call of each segment(red)

cnvr(resCNMOPS) #look at CNV regions

cnvs(resCNMOPS) #look at individual CNV regions


#-------------------------------------------------------------------------------
test_cnv_cn.mops <-
  function()
{	
    BamFiles <- list.files(system.file("extdata", package="cn.mops") , 
                     pattern=".bam$", full.names=TRUE)
    bamDataRanges <- getReadCountsFromBAM(BamFiles, 
                         sampleNames=paste("Sample",1:2), mode="unpaired")
    checkEquals(856, bamDataRanges$Sample.1[1])
    
    data(cn.mops)
    got <- cn.mops(XRanges[,1:3])
    checkEquals(6,length(cnvs(got)))
    checkEquals(1775001 ,start(ranges(cnvr(got))[1]))
    checkEquals(1850000, end(ranges(cnvr(got))[1]))
    checkEquals(75000, width(ranges(cnvr(got))[1]))
}