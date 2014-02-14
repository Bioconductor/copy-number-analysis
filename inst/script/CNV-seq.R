#Here we show an implementation of CNV-seq

##step 1 - includes generating best-hit location files for each mapped 
##sequence read. The authors provide a perl script for BLAT psl file 
## and SOLiD maching pipeline. For BAM files, they suggest to extract 
##locations using the following command

#sarora@rhino01:~/copynumber$ samtools view -F 4 tumorA.chr4.bam |perl -lane 
#'print "F[2]\t$F[3]"' >tumor.hits

#sarora@rhino01:~/copynumber$ samtools view -F 4 normalA.chr4.bam |perl 
#-lane 'print "F[2]\t$F[3]"' >normal.hits


##cnv-seq.pl is used to calculate sliding window size, to count number of 
##mapped hits in each window, and to call cnv R package to calculate log2 
## ratios and annotate CNV

# perl cnv-seq.pl --test tumor.hits --ref normal.hits --genome human


##two output files are produced. They can be found under "result files":
##tumor.hits-vs-normal.hits.window-10000.minw-4.cnv
##tumor.hits-vs-normal.hits.window-10000.minw-4.count

#One can visualize the cnv inside R using the following code snippet
## the plot can be found under "image" folder

library(cnv)
data <- read.delim("tumor.hits-vs-normal.hits.window-10000.minw-4.cnv")
cnv.summary(data)
png("cnv-seq-plot.png")
plot.cnv(data)
dev.off()