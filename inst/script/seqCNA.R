library(seqCNA)

tumor <- file.path(datafiles, "tumorA.chr4.bam")
normal <- file.path(datafiles, "normalA.chr4.bam")
##runSeqsumm is run only once.
#runSeqsumm(summ.win=10,
#           file=tumor,
#           output.file="tumorA_chr4_seqsumm_out.txt",
#           samtools.path="samtools")

# runSeqsumm(summ.win=10,
#            file=normal,
#            output.file="normalA_chr4_seqsumm_out.txt",
#            samtools.path="samtools")

##note: by default it has written the two files in the same folder
## as the bam files!!

tumordata<- read.table(file.path(datafiles,"tumorA_chr4_seqsumm_out.txt"),
                       sep="\t", header=TRUE)
normaldata<- read.table(file.path(datafiles,"normalA_chr4_seqsumm_out.txt"),
                       sep="\t", header=TRUE)
tumordata<- tumordata[which(tumordata$chrom=="chr4"),]
normaldata<- normaldata[which(normaldata$chrom=="chr4"),]
head(tumordata)

rco = readSeqsumm(build="hg19",
                    tumour.data=tumordata,
                    normal.data=normaldata)

##apply the trimming and mapping quality filters
rco =applyFilters(rco, trim.filter=1, mapq.filter=2)
rco= runSeqnorm(rco)
rco=runGLAD(rco)

plotCNProfile(rco)

rco = applyThresholds(rco, seq(-0.8,4,by=0.8), 1)

plotCNProfile(rco)

summary(rco)
