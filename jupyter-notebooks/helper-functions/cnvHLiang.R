cnv.cal <- function(file, log2.threshold=NA, chromosomal.normalization=TRUE, annotate=FALSE, minimum.window=4)
{
	if(is.na(log2.threshold))
	{
		log2.threshold <- as.numeric(sub('\\.pvalue-.+$','', sub('^.+\\.log2-', '', file, perl=TRUE), perl=TRUE))
	}
	data <- read.delim(file)
	data$position <- round((data$end+data$start)/2)
	data$log2 <- NA
	data$p.value <- NA
	test.lambda <- list()
	ref.lambda <- list()
	chrom <- as.character(unique(data$chromosome))
	if(chromosomal.normalization)
	{
		lambda.test.allchr <- mean(data$test)
		lambda.ref.allchr <- mean(data$ref)
		for(chr in chrom)
		{
			sub <- subset(data, chromosome==chr)
			lambda.test <- mean(sub$test)
			lambda.ref <- mean(sub$ref)

			if(lambda.test==0 | lambda.ref==0){
				lambda.test <- lambda.test.allchr
				lambda.ref <- lambda.ref.allchr
			}

			test.lambda[[chr]] <- lambda.test
			ref.lambda[[chr]] <- lambda.ref
			norm <- lambda.test/lambda.ref
			ratio <- sub$test/sub$ref
			sub$log2 <- log2(sub$test/sub$ref/norm)
			data[rownames(sub), 'log2'] <- sub$log2

			subp <- subset(sub, log2>=0)
			subp$ref[subp$ref==0]<- 0.1
			if(nrow(subp)>0) data[rownames(subp),'p.value'] <- pnorm(z2t(subp$test/subp$ref,lambda.test,lambda.ref), lower.tail=FALSE)
			subn <- subset(sub, log2<0)
			if(nrow(subn)>0) data[rownames(subn),'p.value'] <- pnorm(z2t(subn$test/subn$ref,lambda.test,lambda.ref), lower.tail=TRUE)
		}


	}else
	{
		lambda.test <- mean(data$test)
		lambda.ref <- mean(data$ref)
		for(chr in chrom)
		{
			test.lambda[[chr]]<-lambda.test
			ref.lambda[[chr]]<-lambda.ref
		}
		norm <- lambda.test/lambda.ref
		ratio <- data$test/data$ref
		data$log2 <- log2(ratio/norm)
		subp <- subset(data, log2>=0)
		data[rownames(subp),'p.value'] <- pnorm(z2t(subp$test/subp$ref,lambda.test,lambda.ref), lower.tail=FALSE)
		subn <- subset(data, log2<0)
		data[rownames(subn),'p.value'] <- pnorm(z2t(subn$test/subn$ref,lambda.test,lambda.ref), lower.tail=TRUE)
	}
	if(annotate)
	{
		data$cnv <- 0
		data$cnv.size <- NA
		data$cnv.log2 <- NA
		data$cnv.p.value <- NA

		window.size <- data[1,'end']-data[1,'start']+1;
		step <- window.size/2;
		for(chr in chrom)
		{
			print(paste('chromosome: ',chr))
			sub <- subset(data, chromosome==chr)
			cnv.p <- cnv.ANNO(subset(sub,log2>=abs(log2.threshold)), max(sub$cnv, data$cnv), step, minimum.window, log2.threshold)
			if(nrow(subset(cnv.p, cnv>0))) 
			{
				data[rownames(cnv.p),'cnv'] <- cnv.p$cnv
			}
			cnv.n <- cnv.ANNO(subset(sub,log2<=-1*abs(log2.threshold)), max(sub$cnv, data$cnv), step, minimum.window, log2.threshold)
			if(nrow(subset(cnv.n, cnv>0))) 
			{
				data[rownames(cnv.n),'cnv'] <- cnv.n$cnv
			}
		}

		if(max(data$cnv)>0)
		{
			for(id in seq(1,max(data$cnv)))
			{
				print(paste('cnv_id: ', id, ' of ', max(data$cnv)))
				sub <- subset(data, cnv==id)
				start <- ceiling(mean(c(min(sub$start), min(sub$position))))
				end <- floor(mean(c(max(sub$end), max(sub$position))))
				size <- end - start +1
				chr <- as.character(unique(sub$chromosome))
				ratio <- sum(sub$test)/sum(sub$ref)
				norm <- test.lambda[[chr]]/ref.lambda[[chr]]
				cnv.p.value <- 2*pnorm(z2t(ratio, test.lambda[[chr]]*nrow(sub), ref.lambda[[chr]]*nrow(sub)), lower.tail=ratio<norm)
				index <- which(data$cnv==id)
				data[index,'cnv.log2'] <- log2(ratio/norm)
				data[index,'cnv.p.value'] <- cnv.p.value
				data[index,'cnv.size'] <- size
			}
		}
	}
	data
}

z2t <- function(z, lambdax, lambday)
{
	(lambday*z-lambdax)/sqrt(lambday*z^2+lambdax)
}

cnv.ANNO <- function(cnv.sub, anno.last, step, minimum.window, threshold)
{
	if(nrow(cnv.sub)>1)
	{
		distance <- cnv.sub[2:nrow(cnv.sub),'start']-cnv.sub[1:(nrow(cnv.sub)-1),'start']
		cnv.end <- c(which(distance>step),nrow(cnv.sub))
		last <- 1
		for(i in cnv.end)
		{
			if(i-last+1>=minimum.window)
			{
				anno.last <- anno.last+1
				cnv.sub[last:i,'cnv'] <- anno.last
			}
			last <- i+1
		}
	}else
	{
		if(nrow(cnv.sub)==1)
		{
			if(minimum.window<=1)
			{
				anno.last <- anno.last+1
				cnv.sub$cnv <- anno.last
			}
		}
	}
	cnv.sub
}

cnv.print <- function(cnv, file="")
{
	cat('cnv', 'chromosome', 'start', 'end', 'size', 'log2', 'p.value', sep="\t", file=file,fill=TRUE, append=TRUE)
	for(i in seq(max(min(cnv$cnv),1), max(cnv$cnv)))
	{
		sub <- subset(cnv, cnv==i)
		start <- ceiling(mean(c(min(sub$start), min(sub$position))))
		end <- floor(mean(c(max(sub$end), max(sub$position))))
		cat(paste('CNVR_',i,sep=''), paste('chr', unique(sub$chromosome), sep=''), start, end, end-start+1, unique(sub$cnv.log2), unique(sub$cnv.p.value), sep="\t", file=file, fill=TRUE, append=TRUE)
	}
}

cnv.summary <- function(cnv)
{
	true <- subset(cnv, cnv>0)
	
	max.size <- max(true$cnv.size)
	min.size <- min(true$cnv.size)
	count <- length(unique(true$cnv))
	nt.size <- 0
	for(i in unique(true$cnv)) nt.size <- nt.size + max(true[which(true$cnv==i),]$cnv.size)
	mean.size <- round(nt.size/count, 0)
	median.size <- median(unique(true$cnv.size))
	percen <- round(100*nrow(true)/nrow(cnv),1)
	
	cat('CNV percentage in genome: ', percen, "%\n", sep='')
	cat('CNV nucleotide content:', nt.size, fill=TRUE)
	cat('CNV count:', count, fill=TRUE)
	cat('Mean size:', mean.size, fill=TRUE)
	cat('Median size:', median.size, fill=TRUE)
	cat('Max Size:', max.size, fill=TRUE)
	cat('Min Size:', min.size, fill=TRUE)
}

plot.cnv <- function(data, chromosome=NA, CNV=NA, ...)
{
	if(!is.na(CNV))
	{
		plot.cnv.cnv(data, CNV, ...)
	}else
	{
		if(length(unique(data$chromosome))==1 || !is.na(chromosome)) plot.cnv.chr(data, chromosome, ...)
		else plot.cnv.all(data, ...)
	}
} 

plot.cnv.all <- function(data, chrom.gap=2e7, colour=5, title=NA, ylim=c(-2,2), xlabel='Chromosome')
{
	if(nrow(subset(data, log2>max(ylim) | log2<min(ylim)))>0) warning('missed some data points due to small ylim range')
	level <- sort(as.character(unique(data$chromosome)))
	ok.to.convert <- !is.na(suppressWarnings(as.numeric(level)))
	level <- c(sort(as.numeric(level[ok.to.convert])), level[!ok.to.convert])
	labels<-breaks<-c()
	last.pos <- 0
	temp <- data.frame()
	for(chr in level)
	{
		sub <- subset(data, chromosome==chr)
		if(nrow(sub)==0) next
		sub$gpos <- last.pos+sub$position
		labels <- c(labels, chr)
		breaks <- c(breaks, (min(sub$gpos)+max(sub$gpos))/2)
		last.pos <- max(sub$gpos)+chrom.gap
		temp <- rbind(temp, sub)
	}
	data <- temp
	data$colour <- as.character(match(data$chromosome,level)%%colour)
	p <- ggplot()+geom_point(data=data, map=aes(x=gpos, y=log2, colour=colour))
	p <- p + scale_y_continuous(expression(paste(Log[2],' Ratio')), lim=ylim)
	p <- p + scale_x_continuous(xlabel, breaks=breaks, labels=labels)
#	cat(breaks,"\n", sep=',')
#	cat(labels,"\n",sep=',')
	if(colour==2) p <- p + scale_colour_discrete(l=30)
	else p <- p + scale_colour_brewer(palette='Set1')
	p$legend.position <- 'none'
	if(!is.na(title)) p$title <- title
	p
}

plot.cnv.chr <- function(data, chromosome=NA, from=NA, to=NA, title=NA, ylim=c(-4,4), glim=c(NA,NA), xlabel='Position (bp)')
{
	chrom.name <- chromosome
	if(!is.na(chrom.name)) data <- subset(data, chromosome==chrom.name)
	if(!is.na(from) & !is.na(to)) data <- subset(data, position>=from & position<=to)
	if(length(unique(data$chromosome))>1) warning('More than one chromosome! use plot.cnv.all() instead?')
	if(nrow(subset(data, log2>max(ylim) | log2<min(ylim)))>0) warning('missed some data points due to small ylim range')
	if(is.na(glim[1])) glim[1] <- median(subset(data, log2>=0.6 & log2<=1)$p.value)
	if(is.na(glim[1])) glim[1] <- quantile(data$p.value, 0.8, na.rm = T)
	if(glim[1] == 0) glim[1] <- min(data$p.value[data$p.value > 0], na.rm = T)
	if(is.na(glim[2])) glim[2] <- max(na.omit(data$p.value))

	p <- ggplot() + geom_point(data=subset(data, p.value >= glim[1]), map=aes(x=position, y=log2, colour=p.value))
	p <- p + scale_y_continuous(expression(paste(Log[2], ' Ratio')), lim=ylim)
	p <- p + scale_x_continuous(xlabel)
	p <- p + scale_colour_gradient(limits=glim, trans='log10', low="red", high="white")
	temp = subset(data, p.value < min(glim))
	if(nrow(temp)>0) p <- p + geom_point(data=temp, map=aes(x=position,y=log2), colour="red")
	if(!is.na(title)) p$title <- title
	p
}

plot.cnv.cnv <- function(data, CNV, upstream=NA, downstream=NA, ...)
{
	sub <- subset(data, cnv==CNV)
	if(nrow(sub)==0)
	{
		return('CNV number given is not avaliable from data')
	}
	if(is.na(upstream)) upstream=5*sub$cnv.size[1]
	if(is.na(downstream)) downstream=5*sub$cnv.size[1]
	chrom = sub$chromosome[1]
	a = min(sub$start)-upstream
	b = max(sub$end)+downstream
	to.plot <- subset(data, chromosome==chrom & position>=a & position<=b)
	plot.cnv.chr(to.plot, chromosome=chrom, ...)
}
