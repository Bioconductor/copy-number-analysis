copy number analysis (CNA)
===========================

Explore, compare, and evaluate Bioconductor packages related to genomic copy number analysis

Genomic amplifications and deletions are found in most (all?) tumor genomes.  A common practice today is to do low coverage DNA sequencing (0.5x, for instance) of a tumor genome, and a matched normal genome (from the same subject).  Judicious comparison of the the two sequence genomes illuminates structural changes in the tumor.

Copy number changes in tumors vary from broad (an entire chromosome arm) to focal (i.e., a 10kb amplification, loss of heterozygosity or gain).   Detection methods should be sensitive enough to detect these very different phenomena in noisy low-coverage data.

Our purpose here is to provide

* A tumor/normal single chromosome pair of bam files (with accompanying index files)
* A reference analysis, using the popular SeqSeg matlab program from the Broad Institute
* A tutorial on the exploratory data analysis of these files using "native" Bioconductor capabilities
* Demonstrate (and evaluate) the capabilities of many of the Bioconductor copy number analysis packages

