library(CNAnorm)
library(RUnit)

#outside R session - convert perl file
##perl bam2windows.pl "tumorA.chr4.bam" "normalA.chr4.bam" > perloutput.txt

test_counts <- function(chr4data)
{
    test_indices <- c(8000001, 8010001,10000001 , 10010001, 1000001)
    test_res <- c(67,62,67,74,47) #from samtools
    res <- sapply( test_indices, function(x) chr4data[which(chr4data[,2]==x),3])
    checkEquals(test_res,res)
}

data<- read.table("perloutput.txt",sep="\t",header=TRUE)

#subset to chromosome 4
chr4data <- data[which(data[,1]=="chr4"),]

#check if the raw counts are similiar to counts from samtools
test_counts() #FALSE!

#create an object of class CNAnorm
cn <- dataFrame2object(chr4data)

#smooth the signal to decrease noise without losing resolution.
cn <- addSmooth(cn,lambda=7)

#estimate peaks and ploidy
cn <- peakPloidy(cn) 

#produce a plot
png("cna-norm.png")
plotPeaks(cn)
dev.off()