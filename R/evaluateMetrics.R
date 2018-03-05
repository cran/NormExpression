gridAUCVC <- function(data, dataType=c("bk","sc"), HG7=NULL, ERCC=NULL, TN=NULL, TC=NULL, CR=NULL, NR=NULL, DESeq=NULL, 
					UQ=NULL, TMM=NULL, TU=0, GAPDH=NULL, nonzeroRatios=c(0.7,0.8,0.9,1), cvNorm=TRUE, cvResolution=0.005)
{
	grid_result <- NULL;
	if (length(TU)==1 && TU==1){
		colnames_paraMatrix <- c("nonzeroRatio", "pre_ratio", "lower_trim", "upper_trim");
		write.table(t(as.matrix(colnames_paraMatrix)), file = "bestPara.txt", sep = "\t", row.names=FALSE, col.names=FALSE);
	}
	for (i in nonzeroRatios){
		if (dataType=="sc"){
			if ((ncol(data)*i) <= 100){
				cat("nonzeroRatio:", i, " is too small!\n");
				stop("We suggest that the minimal counts of nonzero samples should be greater than 100!");
			} 
		}
		result <- nonzeroRatio2AUCVC(data=data, dataType=dataType, HG7=HG7, ERCC=ERCC, TN=TN, TC=TC, CR=CR, NR=NR, DESeq=DESeq, 
			UQ=UQ, TMM=TMM, TU=TU, GAPDH=GAPDH, nonzeroRatio=i, cvNorm=cvNorm, cvResolution=cvResolution);
		#grid_record <- c(i, result);
		#names(grid_record)[1] <- "NonzeroRatio"; 
		#grid_result <- c(grid_result, names(grid_record), grid_record);
		nonzeroM <- matrix(i,1,1,TRUE);
		colnames(nonzeroM) <- "NonzeroRatio";
		grid_record <- cbind(nonzeroM, result);
		grid_result <- rbind(grid_result, grid_record);
	}
	#aucvcMatrix <- matrix(grid_result, ncol=length(result)+1, byrow=TRUE);
	#return(aucvcMatrix);
	return(grid_result);
}

gridAUCVC4Matrices <- function(..., nonzeroRatios=NULL, cvNorm=TRUE, cvResolution=0.005){ 
	if (is.null(nonzeroRatios)){
		stop("Please provide nonzeroRatios!");
	}
	matrices <- list(...);
	numMethod <- length(matrices);
	grid_result <- NULL;
	for (i in nonzeroRatios){
		result.sorted <- getAUCVCs(..., nonzeroRatio=i, cvNorm=cvNorm, cvResolution=cvResolution);
		grid_record <- c(i, result.sorted);
		names(grid_record)[1] <- "NonzeroRatio"; 
		grid_result <- c(grid_result, names(grid_record), grid_record);
	}
	grid_result2 <- matrix(grid_result, ncol=numMethod+1, byrow=TRUE);
	return(grid_result2);
}

getCorMedians <- function(data){
	if(!is.data.frame(data)) data <- data.frame(data);
	if(is.factor(data$Value)) data$Value <- as.numeric(as.character(data$Value)); # attention!!!
	sorted_result <- sort(tapply(data$Value,data$Methods,median), decreasing=FALSE);
	return(sorted_result);
}

gatherCors <- function(data, cor_method=c("spearman", "pearson", "kendall"), HG7=NULL, ERCC=NULL, TN=NULL, TC=NULL, CR=NULL, NR=NULL, 
					DESeq=NULL, UQ=NULL, TMM=NULL, TU=NULL, GAPDH=NULL, pre_ratio=0.5, lower_trim=0.05, upper_trim=0.65, rounds=1000000)
{
	methodsList <- list(HG7=HG7, ERCC=ERCC, TN=TN, TC=TC, CR=CR, NR=NR, DESeq=DESeq, UQ=UQ, TMM=TMM, TU=TU, GAPDH=GAPDH);
	specifiedMethods <- methodsList[!unlist(lapply(methodsList,is.null))];
	numMethod <- length(specifiedMethods);
	method_range <- seq(1, numMethod, 1);
	ubq_genes <- identifyUbq(data, pre_ratio=pre_ratio, lower_trim=lower_trim, upper_trim=upper_trim, min_ubq=100);
	cor_value_method <- NULL;
	for (j in method_range){
		norm.matrix <- getNormMatrix(data, specifiedMethods[[j]]);
		dataUse2Cor <- norm.matrix[ubq_genes,];
		cor.result <- getCor(dataUse2Cor, method=cor_method, rounds=rounds);
		cor_vm <- cbind(cor.result, rep(names(specifiedMethods)[j], times=round(rounds)));
		cor_value_method <- rbind(cor_value_method, cor_vm);
	}
	colnames(cor_value_method) <- c("Value", "Methods");
	return(cor_value_method);
}

gatherCors4Matrices <- function(..., raw_matrix, cor_method=c("spearman", "pearson", "kendall"), 
								pre_ratio=0.5, lower_trim=0.05, upper_trim=0.65, rounds=1000000)
{
	matrices <- list(...);
	numMethod <- length(matrices);
	method_range <- seq(1, numMethod, 1);
	ubq_genes <- identifyUbq(raw_matrix, pre_ratio=pre_ratio, lower_trim=lower_trim, upper_trim=upper_trim, min_ubq=100);
	cor_value_method <- NULL;
	for (j in method_range){
		dataUse2Cor <- matrices[[j]][ubq_genes,];
		cor.result <- getCor(dataUse2Cor, method=cor_method, rounds=rounds);
		cor_vm <- cbind(cor.result, rep(names(matrices)[j], times=round(rounds)));
		cor_value_method <- rbind(cor_value_method, cor_vm);
	}
	colnames(cor_value_method) <- c("Value", "Methods");
	return(cor_value_method);
}

gatherCVs <- function(data, nonzeroRatio=NULL, HG7=NULL, ERCC=NULL, TN=NULL, TC=NULL, CR=NULL, NR=NULL, 
					DESeq=NULL, UQ=NULL, TMM=NULL, TU=NULL, GAPDH=NULL, cvNorm=TRUE, cvResolution=0.005)
{
	if (is.null(nonzeroRatio)){
		stop("Please provide nonzeroRatio!");
	}
	methodsList <- list(HG7=HG7, ERCC=ERCC, TN=TN, TC=TC, CR=CR, NR=NR, DESeq=DESeq, UQ=UQ, TMM=TMM, TU=TU, GAPDH=GAPDH);
	specifiedMethods <- methodsList[!unlist(lapply(methodsList,is.null))];
	numMethod <- length(specifiedMethods);
	method_range_tmp <- seq(1, numMethod, 1);
	cv_range_tmp <- seq(0, 1, cvResolution);
	method_range_times <- length(cv_range_tmp);
	cv_range_times <- length(method_range_tmp);
	method_range <- rep(method_range_tmp, each=round(method_range_times));
	cv_range <- rep(cv_range_tmp, times=round(cv_range_times));
	nozeroIndex <- filteredZero(data, nonzeroRatio=nonzeroRatio);
	for (j in method_range_tmp){
		norm.matrix <- getNormMatrix(data, specifiedMethods[[j]]);
		dataUse2CV <- norm.matrix[nozeroIndex,];
		cv.result <- getCV(dataUse2CV, cvNorm=cvNorm);
		assign(paste(names(specifiedMethods)[j], ".cv", sep = ""), cv.result);
	}
	cv_uniform <- NULL;
	cv_uniform_all <- mapply(function(i,j){
		cv.result <- paste(names(specifiedMethods)[j], ".cv", sep = "");
		gene_number <- length(which(get(cv.result)<=i));
		cv_uniform_row <- c(i, gene_number, names(specifiedMethods)[j]);
		rbind(cv_uniform, cv_uniform_row);
	}, cv_range, method_range);
	cv_uniform_all <- t(cv_uniform_all);
	colnames(cv_uniform_all) <- c("Cutoff", "Counts", "Methods");
	return(cv_uniform_all);
}

gatherCVs4Matrices <- function(..., raw_matrix, nonzeroRatio=NULL, cvNorm=TRUE, cvResolution=0.005)
{
	if (is.null(nonzeroRatio)){
		stop("Please provide nonzeroRatio!");
	}
	matrices <- list(...);
	matrices_name <- names(matrices);
	numMethod <- length(matrices);
	method_range_tmp <- seq(1, numMethod, 1);
	cv_range_tmp <- seq(0, 1, cvResolution);
	method_range_times <- length(cv_range_tmp);
	cv_range_times <- length(method_range_tmp);
	method_range <- rep(method_range_tmp, each=round(method_range_times));
	cv_range <- rep(cv_range_tmp, times=round(cv_range_times));
	nozeroIndex <- filteredZero(raw_matrix, nonzeroRatio=nonzeroRatio);
	for (j in method_range_tmp){
		dataUse2CV <- matrices[[j]][nozeroIndex,];
		cv.result <- getCV(dataUse2CV, cvNorm=cvNorm);
		assign(paste(matrices_name[j], ".cv", sep = ""), cv.result);
	}
	cv_uniform <- NULL;
	cv_uniform_all <- mapply(function(i,j){
		cv.result <- paste(matrices_name[j], ".cv", sep = "");
		gene_number <- length(which(get(cv.result)<=i));
		cv_uniform_row <- c(i, gene_number, matrices_name[j]);
		rbind(cv_uniform, cv_uniform_row);
	}, cv_range, method_range);
	cv_uniform_all <- t(cv_uniform_all);
	colnames(cv_uniform_all) <- c("Cutoff", "Counts", "Methods");
	return(cv_uniform_all);
}

gatherFactors <- function(data, methods=c("HG7", "ERCC", "TN", "TC", "CR", "NR", "DESeq", "UQ", "TMM", "TU"), 
						HG7.size=NULL, ERCC.size=NULL, TN.size=NULL, TC.size=NULL, CR.size=NULL, NR.size=NULL,
						pre_ratio=0.5, lower_trim=0.05, upper_trim=0.65, min_ubq=100)
{
	method1 <- as.list(methods);
	numMethod <- length(method1);
	method_range <- seq(1, numMethod, 1);
	for (i in method_range){
		if (method1[[i]] == "HG7"||method1[[i]] == "ERCC"||method1[[i]] == "TN"||method1[[i]] == "TC"||
			method1[[i]] == "CR"||method1[[i]] == "NR")
		{
			size.name <- paste(method1[[i]], ".size", sep = "");
			out.name1 <- paste(method1[[i]], ".factors", sep = "");
			if (is.null(size.name)){
				stop("Please provide", size.name, "!");
			}else{
				assign(out.name1, getFactors(data, method="sizefactor", lib.size=get(size.name)));
			}
		}
		if (method1[[i]] == "DESeq" || method1[[i]] == "RLE" || method1[[i]] == "UQ" || method1[[i]] == "TMM"){
			out.name2 <- paste(method1[[i]], ".factors", sep = "");
			assign(out.name2, getFactors(data, method=method1[[i]]));
		}
		if (method1[[i]] == "TU"){
			TU.factors <- getFactors(data, method="TU", pre_ratio=pre_ratio, lower_trim=lower_trim, upper_trim=upper_trim, min_ubq=min_ubq);
		}
	}
	factors.list <- NULL;
	for (m in methods){
		m.factors <- paste(m, ".factors", sep = "");
		factors.list <- c(factors.list, m.factors);
	}
	factors.result <- NULL;
	for (i in method_range){
		factors.result <- cbind(factors.result, get(factors.list[i]));
	}
	colnames(factors.result) <- methods;
	return(factors.result);
}

optTU <- function(data, nonzeroRatio=NULL, pre_ratio_range=c(0.2,0.6), prResolution=0.1,
						lower_range=c(0.05,0.4), upper_range=c(0.6,0.95), qResolution=0.05, min_ubq=100, cvNorm=TRUE, cvResolution=0.005)
{
	if (is.null(nonzeroRatio)){
		stop("Please provide nonzeroRatios!");
	}
	pre_ratio_times <- (pre_ratio_range[2]-pre_ratio_range[1]+prResolution)*10;
	lower_times <- (upper_range[2]-upper_range[1]+qResolution)/qResolution;
	lower_range_tmp <- rep(seq(lower_range[1], lower_range[2], qResolution), each=round(lower_times));
	lower_range2 <- rep(lower_range_tmp, times=round(pre_ratio_times));
	upper_times <- (lower_range[2]-lower_range[1]+qResolution)/qResolution;
	upper_range_tmp <- rep(seq(upper_range[1] , upper_range[2], qResolution), times=round(upper_times));
	upper_range2 <- rep(upper_range_tmp, times=round(pre_ratio_times));
	lower_upper_tmp_len <- length(lower_range_tmp);
	pre_ratio_range2 <- rep(seq(pre_ratio_range[1], pre_ratio_range[2], 0.1), each=round(lower_upper_tmp_len));
	nozeroIndex <- filteredZero(data, nonzeroRatio=nonzeroRatio);
	all_aucvc <- mapply(function(lower_trim, upper_trim, pre_ratio){
								factors.TU <- getFactors(data, method="TU", lower_trim=lower_trim, 
														upper_trim=upper_trim, pre_ratio=pre_ratio, min_ubq=min_ubq);
								norm.TU <- getNormMatrix(data, factors.TU);
								dataUse2CV <- norm.TU[nozeroIndex,];
								cv.TU <- getCV(dataUse2CV, cvNorm=cvNorm);
								TU.AUCVC <- CV2AUCVC(cv.TU, cvResolution=cvResolution);
								return(c(TU.AUCVC=TU.AUCVC, lower=lower_trim, upper=upper_trim, ratio=pre_ratio));
								}, lower_range2, upper_range2, pre_ratio_range2);
	all_aucvc2 <- t(all_aucvc);
	max_index <- which(max(all_aucvc2[,"TU.AUCVC"])==all_aucvc2[,"TU.AUCVC"]);
	return(all_aucvc2[max_index,]);
}

nonzeroRatio2AUCVC <- function(data, dataType=c("bk","sc"), HG7=NULL, ERCC=NULL, TN=NULL, TC=NULL, CR=NULL, 
	NR=NULL, DESeq=NULL, UQ=NULL, TMM=NULL, TU=0, GAPDH=NULL, nonzeroRatio=NULL, cvNorm=TRUE, cvResolution=0.005)
{
	# get the filtered gene indexes
	nozeroIndex <- filteredZero(data, nonzeroRatio=nonzeroRatio);
	methodsList <- list(HG7=HG7, ERCC=ERCC, TN=TN, TC=TC, CR=CR, NR=NR, DESeq=DESeq, UQ=UQ, TMM=TMM, TU=TU, GAPDH=GAPDH);
	specifiedMethods <- methodsList[!unlist(lapply(methodsList,is.null))];
	if (length(TU)==1 && TU==0){
		specifiedMethods$TU <- NULL;
	}
	if (length(TU)==1 && TU==1){
		if (dataType == "bk"){
			optimalPara <- optTU(data, nonzeroRatio=nonzeroRatio,  pre_ratio_range=c(1,1), prResolution=0.1, 
			lower_range=c(0.05,0.4), upper_range=c(0.6,0.95), qResolution=0.05, min_ubq=1000, cvNorm=cvNorm, cvResolution=cvResolution);
		}else{
			optimalPara <- optTU(data, nonzeroRatio=nonzeroRatio, pre_ratio_range=c(0.2,0.6), prResolution=0.1, 
			lower_range=c(0.05,0.4), upper_range=c(0.6,0.95), qResolution=0.05, min_ubq=100, cvNorm=cvNorm, cvResolution=cvResolution);
		}
		optimalPara <- as.matrix(optimalPara);
		lower_trim <- optimalPara["lower",1]; 
		upper_trim <- optimalPara["upper",1];
		pre_ratio <- optimalPara["ratio",1];
		para <- c(nonzeroRatio, pre_ratio, lower_trim, upper_trim);
		names(para)[1] <- "nonzeroRatio";
		paraMatrix <- t(as.matrix(para));
		write.table(paraMatrix, file = "bestPara.txt", sep = "\t", row.names=FALSE, col.names=FALSE, append = TRUE);
		# get the optimal TU normalization factors
		TU.factors <- getFactors(data, method="TU", lower_trim=lower_trim, upper_trim=upper_trim, pre_ratio=pre_ratio, min_ubq=100);
		norm.matrix <- getNormMatrix(data, TU.factors);
		dataUse2CV <- norm.matrix[nozeroIndex, ];
		cv.result <- getCV(dataUse2CV, cvNorm=cvNorm);
		TU.AUCVC <- CV2AUCVC(cv.result, cvResolution=cvResolution)
		specifiedMethods$TU <- NULL;
	}
	numMethod <- length(specifiedMethods);
	if (numMethod>=1){
		method_range <- seq(1, numMethod, 1);
		for (i in method_range){
			norm.matrix <- getNormMatrix(data, specifiedMethods[[i]]);
			dataUse2CV <- norm.matrix[nozeroIndex, ];
			cv.result <- getCV(dataUse2CV, cvNorm=cvNorm);
			assign(names(specifiedMethods)[i] ,CV2AUCVC(cv.result, cvResolution=cvResolution));
		}
		AUCVC.result <- NULL;
		for (i in method_range){
			AUCVC.result <- cbind(AUCVC.result, get(names(specifiedMethods)[i]));
		}
		colnames(AUCVC.result) <- names(specifiedMethods);
		if (length(TU)==1 && TU==1){
			AUCVC.result <- cbind(AUCVC.result, TU.AUCVC); # add TU
			colnames(AUCVC.result) <- c(names(specifiedMethods), "TU");
		}
	}
	if (numMethod==0 && TU==0) stop("Please specify at least one method!");
	if (numMethod==0 && TU==1){ # only TU
		AUCVC.result <- as.matrix(TU.AUCVC); # add TU
		colnames(AUCVC.result) <- "TU";
	}
	#result.sorted <- sort(colSums(AUCVC.result),decreasing=TRUE);
	#return(result.sorted);
	return(AUCVC.result);
}

getAUCVCs <- function(..., nonzeroRatio=NULL, cvNorm=TRUE, cvResolution=0.005){ 
	matrices <- list(...);
	numMethod <- length(matrices);
	method_range <- seq(1, numMethod, 1);
	result <- NULL;
	for (i in method_range){
		AUCVC.result <- getAUCVC(matrices[[i]], nonzeroRatio=nonzeroRatio, cvNorm=cvNorm, cvResolution=cvResolution);
		result <- c(result, AUCVC.result);
		names(result)[i] <- names(matrices)[i];
	}
	sorted_AUCVCs <- sort(result, decreasing=TRUE);
	return(sorted_AUCVCs);
}

getAUCVC <- function(data, nonzeroRatio=NULL, cvNorm=TRUE, cvResolution=0.005){
	nozeroIndex <- filteredZero(data, nonzeroRatio=nonzeroRatio);
	dataUse2CV <- data[nozeroIndex, ];
	cv.result <- getCV(dataUse2CV, cvNorm=cvNorm);
	CV2AUCVC(cv.result, cvResolution=cvResolution)
}

filteredZero <- function(data, nonzeroRatio){
	nozeroCount <- apply(data, 1, function(x) length(which(x!=0)));
	geneIndex <- which(nozeroCount>=ncol(data)*nonzeroRatio);
	return(geneIndex);
}

getCor <- function(data, method=c("spearman", "pearson", "kendall"), rounds=1000000){
	sp_result <- NULL;
	method <- match.arg(method);
	for (i in 1:rounds){
		rg1 <- sample(1:nrow(data),size=1);
		rg2 <- sample(1:nrow(data),size=1);
		while (rg1 == rg2){
			rg2 <- sample(1:nrow(data),size=1);
		}
		gene1 <- unlist(data[rg1,]);
		gene2 <- unlist(data[rg2,]);
		sp_value <- cor(gene1, gene2, method=method);
		sp_result <- c(sp_result, sp_value);
	}
	return(sp_result);
}

getCV <- function(data, cvNorm=TRUE){
	if(!is.matrix(data)) data <- as.matrix(data);
	if (cvNorm){
		rawCV <- apply(data, 1, function(x){sd(log2(x[x!=0]))/mean(log2(x[x!=0]))});
		#new_cv <- rawCV[-which(is.na(rawCV))];
		#(new_cv-min(new_cv))/(max(new_cv)-min(new_cv));
		(rawCV-min(rawCV))/(max(rawCV)-min(rawCV));
	}else{
		apply(data, 1, function(x){sd(x)/mean(x)});
	}
}

CV2AUCVC <- function(data, cvResolution=0.005){
	cv_cutoff <- NULL;
	uniform_genes_counts <- NULL;
	for (i in seq(0, 1, cvResolution)){
		cv_cutoff <- c(cv_cutoff, i);
		gene_number <- length(which(data<=i));
		uniform_genes_counts <- c(uniform_genes_counts, gene_number);
	}
	getArea(cv_cutoff, uniform_genes_counts);
}

# This function refers to the trapz function from pracma package
getArea <- function(x, y) {
	x <- x/max(x);
	y <- y/max(y);
	if (!(is.numeric(x) || is.complex(x)) ||
            !(is.numeric(y) || is.complex(y)) ){
        stop("Arguments 'x' and 'y' must be real or complex vectors.");
	}
	if (length(x) != length(y)){
        stop("The length of two input vectors should be equal!");
	}
	m <- length(x);
	n <- 2*m;
	xp <- c(x, x[m:1]);
    yp <- c(numeric(m), y[m:1]);
    p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1];
    p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n];
    return(0.5*(p1-p2));
}
