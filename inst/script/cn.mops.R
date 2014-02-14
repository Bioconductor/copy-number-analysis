library(cn.mops)
library(RUnit)

tumor_gr <- getReadCountsFromBAM("tumorA.chr4.bam",
     refSeqName="chr4",WL=10000,mode="unpaired")
normal_gr <- getReadCountsFromBAM("normalA.chr4.bam",
     refSeqName="chr4",WL=10000,mode="unpaired")

# We need a special normalization because the tumor has made large CNVs
X <- tumor_gr
values(X) <- cbind(values(tumor_gr),values(normal_gr)) 
X <- normalizeGenome(X,normType="mode")
 
# Parameter settings for tumor: 
# - norm=0, because we already have normalized
# - integer copy numbers higher than 8 allowed
# - DNAcopy as segmentation algorithm.
ref_analysis_norm0 <- referencecn.mops(X[,1], X[,2],
     norm=0, 
     I = c(0.025, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 8,16,32,64), 
     classes = c("CN0", "CN1", "CN2", "CN3", "CN4", "CN5",
                    "CN6", "CN7","CN8","CN16","CN32","CN64","CN128"),
     segAlgorithm="DNAcopy")


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
