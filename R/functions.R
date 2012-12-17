# computes the composite kernel, that aggregates proximity and prior clustering info
computeCompositeKernel <- function(data, alpha=3, method="neighbors", pgaussian=TRUE) {
	# controls
	if (alpha < 1) {
		stop("alpha should be greater than 1")
	}

	if (!any(method == c("neighbors", "simple"))) {
		stop("method should be either 'neighbors' or 'simple'")
	}
	
	nelts <- dim(data)[1]
	d <- dim(data)[2] - 1 # last column is reserved for the labels
	Ksim <- matrix(0,nrow=nelts, ncol=nelts)
	
	# first, standardize data
	means <- apply(data[,1:d], 2, mean)
	sds <- apply(data[,1:d], 2, sd)
	sds[sds==0] <- 1 # manage NULL case
	data[,1:d] <- sweep(data[,1:d], 2, means, "-")
	data[,1:d] <- sweep(data[,1:d], 2, sds, "/")
	
	# compute p and sigma parameters, that will be used later for the p-gaussian kernel	
	# (see Francois2005 for details)
	# NB : D_N = 5% quantile, and D_F = 95% quantile (necessary for p to remain positive)
	# from there, expressions for sigma seem to have been swapped (D_N with 0.95 and reciprocally)
	# => add as a footnote
	if (pgaussian) {
		dists <- dist(data[,1:d]) # contains all unique non-trivial pairwise distances (n.(n-1)/2)
		quants <- quantile(dists, probs=c(0.05, 0.95))
		p <- log(log(0.05)/log(0.95)) / log(quants[2] / quants[1])
		sigma <- quants[2] / (-log(0.05))^(1/p)
		# should also equal quants[1] / (-log(0.95))^(1/p)		
	} else {
		# alternate : sigma propto maximal distance between any two elts (see Haykin1999).
		p <- 2
		dists <- dist(data[,1:d])
		sigma <- max(dists)
	}

	setSimKernelChunk <- function(data) {	
		Ksim <- sweep(data[,1:d,drop=FALSE], 2, data[1,1:d], "-")
		Ksim <- apply(Ksim, 1, function(x) sqrt(sum(x^2)))
		Ksim <- exp((-Ksim^p) / sigma^p)
		return(Ksim)
	}
		
	getMaxLabeled <- function(sim, labeledexamples) {
		themin <- which.max(sim[labeledexamples])
		themin <- labeledexamples[themin]
		return(themin)
	}

	modKernel <- function(data, labs) {
		sim <- data[1:nelts]
		curlab <- data[nelts+1]
		match <- (!is.na(labs)) & (!is.na(curlab)) & (labs == curlab)
		nomatch <- (!is.na(labs)) & (!is.na(curlab)) & (labs != curlab)
		
		sim[match] <- sim[match]^(1 / alpha)
		sim[nomatch] <- sim[nomatch]^alpha
		return(sim)
	}
		
	
	# first compute standard Kernel similarity matrix
	# avoiding unnecessary calculations
	for (i in 1:nelts) {
		Ksim[i,i:nelts] <- Ksim[i:nelts,i] <- setSimKernelChunk(data[i:nelts,,drop=FALSE]) # handle 1-elt case
	}
	
	# for each element, compute the closest labeled example
	labelvec <- data[,d+1]
	labeledexamples <- which(!is.na(labelvec))
	Kmod <- Ksim
	# find closest indexes of all elements
	closestindex <- apply(Ksim, 1, function(x) getMaxLabeled(x, labeledexamples))
	closestlab <- labelvec[closestindex]
	if (length(labeledexamples) > 0) {
	
		Kmod <- switch(method, 
			# reestablish apply usage - 1 applies row wise, and 2 column wise, but result is always returned as columns.
			simple=t(apply(cbind(Kmod, labelvec), 1, function(x) modKernel(x, labelvec))),
			neighbors=t(apply(cbind(Kmod, closestlab), 1, function(x) modKernel(x, closestlab))))
	}
	return(Kmod)
}

# computes the kernel projection of given data
computeKernelProjection <- function(data, dims=2, alpha=3, method="neighbors", pgaussian=TRUE) {
	nelts <- dim(data)[1]
	if (dims > nelts) stop("dims cannot be greater than data rank")
	kernel <- computeCompositeKernel(data, alpha=alpha, method=method, pgaussian=pgaussian)
	ones <- matrix(1/nelts, nrow=nelts, ncol=nelts)
	kernel <- kernel - ones%*%kernel - kernel%*%ones + ones %*%kernel%*%ones # normalize according to Bishop (12.85)
	eigendecomp <- eigen(kernel, symmetric=TRUE)
	# restric processing to dimensions actually used (ie dims)
	eigendecomp$values <- eigendecomp$values[1:dims]
	eigendecomp$vectors <- eigendecomp$vectors[,1:dims,drop=FALSE]
	# normalize eigenvectors according to (12.81)
	eigendecomp$vectors <- sweep(eigendecomp$vectors, 2, sqrt(eigendecomp$values*nelts), "/") 
	# then compute projections according to (12.82), in matrix form
	projection <- kernel %*% eigendecomp$vectors # transposed version of Bishop's description
}


# perform the kernel PCA projection step from an existing kernel matrix
# see computeKernelProjection for details.
computeProjectionFromKernel <- function(kernel, dims=2) {
	nelts <- dim(kernel)[1]
	if (dims > nelts) stop("dims cannot be greater than data rank")
	ones <- matrix(1/nelts, nrow=nelts, ncol=nelts)
	kernel <- kernel - ones%*%kernel - kernel%*%ones + ones %*%kernel%*%ones
	eigendecomp <- eigen(kernel, symmetric=TRUE)
	eigendecomp$values <- eigendecomp$values[1:dims]
	eigendecomp$vectors <- eigendecomp$vectors[,1:dims,drop=FALSE]
	eigendecomp$vectors <- sweep(eigendecomp$vectors, 2, sqrt(eigendecomp$values*nelts), "/") 
	projection <- kernel %*% eigendecomp$vectors
}


computeStandardKernel <- function(data, pgaussian=TRUE) {
	nelts <- dim(data)[1]
	d <- dim(data)[2]
	Ksim <- matrix(0,nrow=nelts, ncol=nelts)
	# first, standardize data
	means <- apply(data[,1:d], 2, mean)
	sds <- apply(data[,1:d], 2, sd)
	sds[sds==0] <- 1 # manage NULL case
	data[,1:d] <- sweep(data[,1:d], 2, means, "-")
	data[,1:d] <- sweep(data[,1:d], 2, sds, "/")

	# see computeCompositeKernel for explanations of this part
	if (pgaussian) {
		dists <- dist(data[,1:d]) # contains all unique non-trivial pairwise distances (n.(n-1)/2)
		quants <- quantile(dists, probs=c(0.05, 0.95))
		p <- log(log(0.05)/log(0.95)) / log(quants[2] / quants[1])
		sigma <- quants[2] / (-log(0.05))^(1/p)	
	} else {
		p <- 2
		dists <- dist(data[,1:d])
		sigma <- max(dists)
	}

	setSimKernelChunk <- function(data) {	
		Ksim <- sweep(data[,1:d,drop=FALSE], 2, data[1,1:d], "-")
		Ksim <- apply(Ksim, 1, function(x) sqrt(sum(x^2)))
		Ksim <- exp((-Ksim^p) / sigma^p)
		return(Ksim)
	}

	for (i in 1:nelts) {
		Ksim[i,i:nelts] <- Ksim[i:nelts,i] <- setSimKernelChunk(data[i:nelts,,drop=FALSE]) # handle 1-elt case
	}
	return(Ksim)
}


distortions <- function(origdata, projdata) {
	# standardize both data sets
	origmeans <- apply(origdata, 2, mean)
	projmeans <- apply(projdata, 2, mean)
	origsds <- apply(origdata, 2, sd)
	projsds <- apply(projdata, 2, sd)
	origsds[origsds==0] <- 1 # manage NULL case
	projsds[projsds==0] <- 1
	origdata <- sweep(origdata, 2, origmeans, "-")
	projdata <- sweep(projdata, 2, projmeans, "-")
	origdata <- sweep(origdata, 2, origsds, "/")
	projdata <- sweep(projdata, 2, projsds, "/")
	
	# compute distance matrices
	origmat <- as.matrix(dist(origdata))
	projmat <- as.matrix(dist(projdata))
	# normalize in [0,1]
	origmat <- origmat / max(origmat)
	projmat <- projmat / max(projmat)
	
	# compute D^+ and D^- matrices
	posdiffs <- negdiffs <- diffs <- origmat - projmat
	posdiffs[posdiffs<0] <- 0
	negdiffs[negdiffs>0] <- 0
	negdiffs <- -negdiffs
	
	# mu aggregates, and results
	mucompress <- apply(posdiffs, 2, sum)
	mustretch <- apply(negdiffs, 2, sum)
	
	rescompress <- apply(as.data.frame(mucompress), 1, function(x) (x - min(mucompress)) / (max(mucompress) - min(mucompress)))
	resstretch <- apply(as.data.frame(mustretch), 1, function(x) (x - min(mustretch)) / (max(mustretch) - min(mustretch)))
	return(list(compress=rescompress, stretch=resstretch))
}






