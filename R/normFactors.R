getFactors <- function(data, method = c("sizefactor", "DESeq", "RLE", "UQ", "TMM", "TU"), lib.size=NULL,
                      pre_ratio=0.5, lower_trim=0.05, upper_trim=0.65, min_ubq=100)
{	
	if(!is.matrix(data)) data <- as.matrix(data);
	if(any(is.na(data))) stop("NA counts not permitted")
	if(is.null(lib.size)) libsize <- colSums(data) else libsize <- lib.size
	if(any(is.na(libsize))) stop("NA libsizes not permitted")
	
	method <- match.arg(method);
	
	i <- apply(data <= 0, 1, all);
	if (any(i)) data <- data[!i, , drop = FALSE];

	f <- switch(method,
		sizefactor = 1e6/libsize,
		DESeq = 1/estimateSizeFactorsForMatrix(data, p = 0.5),
		RLE = calcFactorRLE(data, p = 0.5) / libsize,
		UQ = calcFactorUpperquartile(data, lib.size = libsize, p = 0.75),
		TMM = {
			fq <- calcFactorUpperquartile(data = data, lib.size = libsize, p = 0.75);
			refColumn <- which.min(abs(fq - mean(fq)));
			if (length(refColumn)==0 | refColumn < 1 | refColumn > ncol(data)) refColumn <- 1;
			f <- rep(NA,ncol(data));
			for(i in 1:ncol(data)){
				f[i] <- calcFactorWeighted(obs=data[,i], ref=data[,refColumn], libsize.obs=libsize[i], libsize.ref=libsize[refColumn], logratioTrim=0.3, sumTrim=0.05, doWeighting=TRUE, Acutoff=-1e10)
			}
			f
		},
		TU = {
			if(!is.data.frame(data)) data <- data.frame(data);
			ubq_genes <- identifyUbq(data, lower_trim=lower_trim, upper_trim=upper_trim, pre_ratio=pre_ratio, min_ubq=min_ubq);
			ubq_sums <- colSums(data[ubq_genes,]);
			mean(ubq_sums)/ubq_sums;
		},
		)
	if (method == "RLE" || method == "UQ" || method == "TMM"){
		f <- 1e6 / libsize / f;
	}
	norm.factors <- f / exp(mean(base::log(f)));
	round(norm.factors, digits = 5);
}

.log <- function(x) {
	x[] <- structure(vapply(x, function(x) if (x == 0) NA else base::log(x), numeric(1)), dim=dim(x))
	return(x)
}

estimateSizeFactorsForMatrix <- function (data, p = p) { 
	loggeomeans <- rowMeans(.log(data), na.rm = TRUE);
	apply(data, 2, function(cnts) exp(quantile(.log(cnts) - loggeomeans, na.rm = TRUE, p = p)));
}

calcFactorRLE <- function (data, p = p) {
	gm <- exp(rowMeans(.log(data), na.rm = TRUE));
	apply(data, 2, function(u) quantile((u / gm)[u!=0], na.rm = TRUE, p = p));
}

calcFactorUpperquartile <- function (data, lib.size, p = p) {
	y <- t(t(data) / lib.size);
	f <- apply(y, 2, function(x) quantile(x[x!=0], p = p));
}

calcFactorWeighted <- function(obs, ref, libsize.obs=NULL, libsize.ref=NULL, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e10)
{
	if (all(obs == ref)) return(1);
	obs <- as.numeric(obs);
	ref <- as.numeric(ref);
	if( is.null(libsize.obs) ) nO <- sum(obs) else nO <- libsize.obs
	if( is.null(libsize.ref) ) nR <- sum(ref) else nR <- libsize.ref
	logR <- log2((obs / nO) / (ref / nR));
	absE <- (log2(obs / nO) + log2(ref / nR)) / 2;
	v <- (nO - obs) / nO / obs + (nR - ref) / nR / ref;
	fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff); 
	logR <- logR[fin];
	absE <- absE[fin];
	v <- v[fin];
	
	n <- length(logR);
	loL <- floor(n * logratioTrim) + 1;
	hiL <- n + 1 - loL;
	loS <- floor(n * sumTrim) + 1;
	hiS <- n + 1 - loS;
	keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >= loS & rank(absE) <= hiS);
	if (doWeighting){
		2 ^ (sum(logR[keep] / v[keep], na.rm = TRUE) / sum(1 / v[keep], na.rm = TRUE))
	}else{
		2 ^ (mean(logR[keep], na.rm = TRUE))
	}
}

identifyUbq <- function(data, pre_ratio=0.5, lower_trim=0.05, upper_trim=0.65, min_ubq=100){
	qlower <- apply(data, 2, function(x) quantile(x[x!=0], p = lower_trim));
	qupper <- apply(data, 2, function(x) quantile(x[x!=0], p = upper_trim));
	ubq_genes <- NULL;
	for(i in 1:nrow(data)){
		genes_finded <- findGenes(data[i,], qlower=qlower, qupper=qupper, pre_ratio=pre_ratio);
		ubq_genes <- c(ubq_genes, genes_finded);
	}
	if (length(ubq_genes) < min_ubq){
		cat("Parameters range", lower_trim, "-", upper_trim, "...identified too few ubiquitous genes (",length(ubq_genes),"), trying range 5-95  instead","\n");
		ubq_genes <- identifyUbqRepeat(data, pre_ratioC=pre_ratio, lower_trimC=0.05, upper_trimC=0.95);
	}
	return(ubq_genes);
}

findGenes <- function(g, qlower=NULL, qupper=NULL, pre_ratio=NULL){
	gene_name <- rownames(g);
	g <- unlist(g);
	seen <- which(g>=qlower & g<=qupper);
	counts <- length(seen);
	if (counts>=pre_ratio*length(g)){
		gene_name;
	}
}

identifyUbqRepeat <- function(data, pre_ratioC=NULL, lower_trimC=NULL, upper_trimC=NULL){
	qlower <- apply(data, 2, function(x) quantile(x[x!=0], p = lower_trimC));
	qupper <- apply(data, 2, function(x) quantile(x[x!=0], p = upper_trimC));
	ubq_genes <- NULL;
	for(i in 1:nrow(data)){
		genes_finded <- findGenes(data[i,], qlower=qlower, qupper=qupper, pre_ratio=pre_ratioC);
		ubq_genes <- c(ubq_genes, genes_finded);
	}
	return(ubq_genes);
}

getNormMatrix <- function(data, norm.factors){
		data*matrix(rep(norm.factors,dim(data)[1]), 
					nrow = dim(data)[1], 
					ncol = length(norm.factors), 
					byrow=T);
}
