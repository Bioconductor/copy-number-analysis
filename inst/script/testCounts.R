## Our motivation for this gist comes from the fact that we would like
## to check the counts calculated by each of the methods.

## Here we provide a simple unitTest which takes as input a GRanges
## object which contains the reads from the "tumorA.chr4.bam" file.

## we choose 5 regions and get their counts from samtools using a simple
## command :
## samtools view  tumorA.chr4.bam chr4:8000001-8010000   | wc -l      67

## we can then use this function to check if  the counts coming out of
## a given method are equal to the counts given by samtools

## we assume that the metadata column containing the counts in the GRanges
## object is called as "cnt"

testCounts <- function(grTumor)
{
    test_indices <- c(8000001, 8010001,10000001 , 10010001, 1000001)
    test_res <- c(67,62,67,74,47) #from samtools
    counts <- sapply(test_indices , 
        function (x) grTumor[which(start(ranges(grTumor))==x)]$cnt  )
    checkEquals(test_res,indices,tolerance=2)
}

