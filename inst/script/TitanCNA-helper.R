###################################################################
##Utility R script to extract tumour allele read counts for input 
##   to TitanCNA. 
## Author: Gavin Ha (gavin.ha@gmail.com)
## Date: May 12, 2014
#####################################################################
## The script makes use of the R Bioconductor package 
##       Rsamtools (>= 1.17.11).  

## The inputs to this script are:
## 1. File path to tumour BAM file (tumbamFile)
## 2. File path to tumour BAM index file (tumIndexFile)
## 3. File path to heterozygous germline SNP positions (vcfFile)
##     This file can be generated using a variety of solutions. 
##     One such solution is posted on 
##     http://compbio.bccrc.ca/software/titan/titan-running-titancna/
##     (see Step 5: Input files).
## The output to this script:
##   A data.frame (countMat).
##   Alternatively, this can be written to a tab-delimited text file with
##   path "outFile". The format of this text file matches that required
##   by TitanCNA ("loadAlleleCountsFromFile()"). 
#####################################################################
#####################################################################

## THE FOLLOWING SCRIPT ONLY WORKS WITH Rsamtools (>= 1.17.11)
library(Rsamtools)

#####################################################################
################ USERS - PLEASE MODIFY THESE PATHS #################
#####################################################################
tumbamFile <- "tumorA.chr4_sorted.bam"
tumIndexFile <- paste(tumbamFile,".bai",sep = "")
vcfFile <-  "titanCNA_sorted_HetSNPs.vcf" ##"titanCNAHetSNPs.vcf" 
outFile <- "paul_tumAlleleCounts.tsv"

#####################################################################
#################### LOAD HET POSITIONS VCF FILE ####################
#####################################################################
## read in vcf file of het positions
vcf <- BcfFile(vcfFile)
vcfPosns <- scanBcf(vcf)

#####################################################################
####################### SETUP PILEUP PARAMETERS #####################
#####################################################################
## setup PileupParam using sequence read filters
pp <- PileupParam(min_base_quality = 10, min_mapq = 20, 
		min_nucleotide_depth = 10, max_depth = 20, 
		distinguish_strands = FALSE, 
		distinguish_nucleotides = TRUE)
## setup the positions of interest to generate the pileup for
which <- GRanges(as.character(vcfPosns$CHROM), 
		IRanges(vcfPosns$POS, width = 1))
## setup addition BAM filters, such as excluding duplicate reads
sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), which = which)

#####################################################################
########################## GENERATE PILEUP ##########################
#####################################################################
## generate pileup using function (Rsamtools >= 1.17.11)
## this step can take a while
tumbamObj <- BamFile(tumbamFile, index = tumIndexFile)
counts <- pileup(tumbamObj, scanBamParam = sbp,  pileupParam = pp)

## set of command to manipulate the "counts" data.frame output
##     by pileup() such that multiple nucleotides are in a single
##     row rather than in multiple rows.
countsMerge <- xtabs(count ~ which_label + nucleotide, counts)
label <- do.call(rbind, strsplit(rownames(countsMerge), ":"))
posn <- do.call(rbind, strsplit(label[, 2],"-"))

countsMerge <- cbind( position = posn[, 1], countsMerge)
mode(countsMerge) <- "numeric"
countsMerge <- data.frame(chr=label[,1],countsMerge, check.names=FALSE)


#####################################################################
############### GET REFERENCE AND NON-REF READ COUNTS ###############
#####################################################################
## this block of code is used to match up the reference and 
##   non-reference nucleotide when assigning read counts
##   final output data.frame is "countMat"
## setup output data.frame
countMat <- data.frame(chr = vcfPosns$CHROM, 
			position = as.numeric(vcfPosns$POS), 
			ref = vcfPosns$REF, refCount = 0, 
			Nref = vcfPosns$ALT, NrefCount = 0, 
			stringsAsFactors = FALSE)

## match rows with vcf positions of interest
countMat <- merge(countMat, countsMerge, by = c("chr","position"), 
		sort = FALSE)

## assign the flattened table of nucleotide counts to ref, Nref
## note that non-reference (Nref) allele is sum of other bases
##    that is not matching the ref.
NT <- c("A", "T", "C", "G")
for (n in 1:length(NT)){	
	indRef <- countMat$ref == NT[n]
	countMat[indRef, "refCount"] <- countMat[indRef, NT[n]]
	countMat[indRef, "NrefCount"] <- rowSums(countMat[indRef, NT[-n]])
}
		
countMat <- countMat[,1:6]

#####################################################################
####################### OUTPUT TO TEXT FILE #########################
#####################################################################
## output text file will have the same format required by TitanCNA
write.table(countMat, file = outFile, row.names = FALSE, 
	col.names = TRUE, quote = FALSE, sep = "\t")
			
