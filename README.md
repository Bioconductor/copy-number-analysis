Copy Number Analysis 
=====================

Explore, compare, and evaluate Bioconductor packages related to genomic copy number analysis

Genomic amplifications and deletions are found in most (all?) tumor genomes.  A common practice today is to do low coverage DNA sequencing (0.5x, for instance) of a tumor genome, and a matched normal genome (from the same subject).  Judicious comparison of the the two sequence genomes illuminates structural changes in the tumor.

Copy number changes in tumors vary from broad (an entire chromosome arm) to focal (i.e., a 10kb amplification, loss of heterozygosity or gain).   Detection methods should be sensitive enough to detect these very different phenomena in noisy low-coverage data.

Our purpose here is to provide

* A tumor/normal single chromosome pair of bam files (with accompanying index files)
* A reference analysis, using the popular SeqSeg matlab program from the Broad Institute
* A tutorial on the exploratory data analysis of these files using "native" Bioconductor capabilities
* Demonstrate (and evaluate) the capabilities of many of the Bioconductor copy number analysis packages

Literature Resources
=========================
* Alkan, C., et al. (2011). <a href="http://www.ncbi.nlm.nih.gov/pubmed/21358748">"Genome structural variation discovery and genotyping."</a> Nat Rev Genet 12(5): 363-376. 
* Duan J, Zhang J-G, Deng H-W, Wang Y-P (2013) <a href="http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0059128">Comparative Studies of Copy Number Variation Detection Methods for Next-Generation Sequencing Technologies.</a> PLoS ONE 8(3): e59128. doi:10.1371/journal.pone.0059128

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


Exploratory Data Analysis
==========================
We have done some primary <a href="https://github.com/Bioconductor/copy-number-analysis.wiki.git">Exploratory Data Analysis</a> on the Normal and Tumor Sample Datasets.

List of Tools used
===================
Bioconductor Packages
* <a href="https://github.com/Bioconductor/copy-number-analysis/wiki/CountOverlaps-method-from-IRanges-Package">countOverlaps</a>
* <a href="https://github.com/Bioconductor/copy-number-analysis/wiki/cn.mops">cn.mops</a>
* <a href="https://github.com/Bioconductor/copy-number-analysis/wiki/CNAnorm">CNAnorm</a>
* <a href="https://github.com/Bioconductor/copy-number-analysis/wiki/seqCNA">seqCNA</a>
<br>

Non Biocondcutor packages
* <a href="https://github.com/Bioconductor/copy-number-analysis/wiki/CNV-seq">CNV-seq</a>
* <a href="https://github.com/Bioconductor/copy-number-analysis/wiki/SegSeq">Seg-seq</a>
