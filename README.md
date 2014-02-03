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

Sample Data
===========
* http://s3.amazonaws.com/copy-number-analysis/tumorA.chr4.bam
* http://s3.amazonaws.com/copy-number-analysis/tumorA.chr4.bam.bai
* http://s3.amazonaws.com/copy-number-analysis/normalA.chr4.bam
* http://s3.amazonaws.com/copy-number-analysis/normalA.chr4.bam.bai

Use, e.g.,
<pre><code> 
download.file(url="http://s3.amazonaws.com/copy-number-analysis/tumorA.chr4.bam.bai",
              destfile="tumorA.chr4.bam.bai")
</code></pre>

gists
=====
* use "native" Bioconductor to get read counts in bins of a specified length across a tumor chromosome: 
 <a  href="https://gist.github.com/pshannon-bioc/8677930">binnedCounts.R</a>

Exploratory Data Analysis
==========================
We have done some primary exploratory data analysis on the Normal and Tumor Sample Datasets using binnedCounts.R
* a plot of raw counts of Tumor and Normal datasets is shown below: ![] (https://raw.github.com/Bioconductor/copy-number-analysis/4ebd109ee23e5630a631104bb845a7daca59abad/image/count2.png)
* a plot of the normalized counts of Tumor and Normal datasets is shown below: ![] (https://raw.github.com/Bioconductor/copy-number-analysis/4ebd109ee23e5630a631104bb845a7daca59abad/image/normalised-count2.png)
