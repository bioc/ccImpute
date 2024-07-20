#' Impute Dropout Values in Single-Cell RNA Sequencing Data 
#' 
#' Performs imputation of dropout values in single-cell RNA sequencing
#' (scRNA-seq) data using a consensus clustering-based algorithm (ccImpute). 
#' This implementation includes performance enhancements over the original
#' ccImpute method described in the paper "ccImpute: an accurate and scalable
#' consensus clustering based algorithm to impute dropout events in the
#' single-cell RNA-seq data" (DOI: https://doi.org/10.1186/s12859-022-04814-8).
#'
#' @param object A \code{SingleCellExperiment} class object containing the
#' scRNA-seq data. The \code{logcounts} assay should contain matrix with
#' log-normalized expression values. This code supports both dense and sparse
#' (dgCMatrix) matrix format storage.
#' @param dist (Optional) A distance matrix used for cell similarity.
#' calculations. If not provided, a weighted Spearman correlation matrix is
#' calculated. 
#' @param nCeil (Optional) The maximum number of cells used to compute the
#' proportion of singular vectors (default: \code{2000}).
#' @param svdMaxRatio (Optional) The maximum proportion of singular vectors
#' used for generating subsets (default: \code{0.08}).
#' @param maxSets (Optional) The maximum number of sub-datasets used for
#' consensus clustering (default: \code{8}).
#' @param k (Optional) The number of clusters (cell groups) in the data. If not
#' provided, it is estimated using the Tracy-Widom Bound.
#' @param consMin (Optional) The low-pass filter threshold for processing the
#' consensus matrix (default: \code{0.75}).
#' @param kmNStart nstart parameter passed to \code{\link[stats]{kmeans}}.
#' function. Can be set manually. By default it is \code{1000} for up to
#' \code{2000} cells and \code{50} for more than \code{2000} cells.
#' @param kmMax iter.max parameter passed to \code{\link[stats]{kmeans}}.
#' @param fastSolver (Optional) Whether to use mean of
#' non-zero values for calculating dropout values or a linear equation solver
#' (much slower and did show empirical difference in imputation performance)
#' (default: \code{TRUE}).
#' @param BPPARAM (Optional) A \code{BiocParallelParam} object for parallel
#' processing (default: \code{bpparam()}).
#' @param verbose (Optional) Whether to print progress messages
#' (default: \code{TRUE}).
#' 
#' @return A \code{SingleCellExperiment} class object with the imputed
#' expression values stored in the `"imputed"` assay.
#' 
#' @examples
#' library(BiocParallel)
#' library(splatter)
#' library(scater)
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' cores <- 2
#' BPPARAM = MulticoreParam(cores)
#' sce <- ccImpute(sce, BPPARAM=BPPARAM)
#' 
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames assay<-
#' @importFrom BiocParallel bplapply bpnworkers bpparam
#' @importFrom sparseMatrixStats rowVars
#' @importFrom Rcpp sourceCpp
#' @useDynLib ccImpute
#' @name ccImpute
#' @aliases ccImpute
#' @export
ccImpute.SingleCellExperiment <- 
    function(object, dist, nCeil = 2000, svdMaxRatio = 0.08,
            maxSets = 8, k, consMin=0.75, kmNStart, kmMax=1000,
            fastSolver = TRUE, BPPARAM=bpparam(), verbose = TRUE){
    startTime <- Sys.time() # Start reference point for elapsed time
    nCores <- bpnworkers(BPPARAM) # Number of cores used for openMP parallelism
    if (!"logcounts" %in% assayNames(object)) {
        stop("The `logcounts` slot in your input SingleCellExperiment
            object is required!")
    }
    
    # Step 1: Get the log normalized expression matrix
    isXSparse <- inherits(logcounts(object), "dgCMatrix")
    logX <- if (isXSparse) logcounts(object) else as.matrix(logcounts(object))
    
    n <- ncol(logX) #number of cells
    printer(verbose, sprintf('Running ccImpute on dataset (%d cells) with %d 
                            cores.\n', n, nCores))
    
    # Step 2: Transform to distance matrix form if dist is not provided
    if(missing(dist)){
        w <-if (isXSparse) sparseMatrixStats::rowVars(logX)
            else rowVars_fast(logX, nCores)
        dist <- getCorM(method = "spearman", logX, w=w, nCores=nCores)
        printer(verbose, "Distance matrix has been computed.", startTime)
    }
    
    # Step 3: Perform truncated SVD on the distance matrix
    v <- doSVD(x=dist, nCores=nCores)
    printer(verbose, "SVD completed.", startTime)
    
    # Step 4: Run k-means on each set
    consMtx <- runKM(logX, v, maxSets, k, consMin, kmNStart, kmMax, BPPARAM)
    printer(verbose, "Clustering completed.", startTime)
    
    # Step 5: Identify Dropouts
    dropIds <- findDropouts(logX, consMtx)
    printer(verbose, "Dropouts identified.", startTime)
    
    # Step 6: Compute Dropouts
    iLogX <- computeDropouts(consMtx, logX, dropIds, fastSolver, nCores)
    printer(verbose, "Dropouts imputed.", startTime)
    assay(object, "imputed") <- iLogX
    return(object)
}

#' ccImpute for SingleCellExperiment Objects
#'
#' Defines the generic function `ccImpute` and a specific method for
#' `SingleCellExperiment` objects.
#'
#' @rdname ccImpute
#'
#' @export
setGeneric("ccImpute", signature = "object", 
            function(object, dist, nCeil = 2000, svdMaxRatio = 0.08,
                    maxSets = 8, k, consMin=0.75, kmNStart, kmMax=1000,
                    fastSolver = TRUE, BPPARAM=bpparam(), verbose = TRUE)
                {standardGeneric("ccImpute")})

#' @rdname ccImpute
setMethod("ccImpute", signature(object = "SingleCellExperiment"), 
            ccImpute.SingleCellExperiment)

#' Calculate Column-wise Correlation Matrix
#'
#' Efficiently computes a column-wise correlation matrix for a 
#' given input matrix. Supports Pearson and Spearman correlations, with
#' optional weighting for features.
#'
#' @param method A character string specifying the correlation 
#' metric to use. Currently supported options are:
#'  - \code{"spearman"}: Spearman rank correlation
#'  - \code{"pearson"}: Pearson correlation
#' @param x A numeric matrix where each column represents a sample.
#' @param w (Optional) A numeric vector of weights for each feature (row) in 
#' \code{x}. If not provided, all features are equally weighted.
#' @param nCores The number of cores to use for parallel processing.
#' 
#' @return A correlation matrix of the same dimensions as the number of columns
#' in `x`. The values represent the pairwise correlations between samples
#' (columns) based on the chosen method and optional weights.
#'
#' @examples
#' library(scater)
#' library(splatter)
#' 
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' cores <- 2
#' logX <- as.matrix(logcounts(sce))
#' w <- rowVars_fast(logX, cores)
#' corMat <- getCorM("spearman", logcounts(sce), w, cores)
#' 
#' @export
getCorM <- function(method, x, w, nCores){
    isXSparse <- inherits(x, "dgCMatrix")
    if (method == "spearman"){
        if(missing(w)){
            if(isXSparse){
                ranked <- sparseColRanks_fast(x, nCores)
                return(cor_fast(ranked, nCores))
            }
            ranked <- colRanks_fast(x, nCores)
            return(cor_fast(ranked, nCores))
        }
        else{
            if(isXSparse){
                ranked <- sparseColRanks_fast(x, nCores)
                return(wCor_fast(ranked, w, nCores))
            }
            ranked <- colRanks_fast(x, nCores)
            return(wCor_fast(ranked, w, nCores))
        }    
    }
    else if (method == "pearson"){
        if(missing(w)){
            if(isXSparse){
                return(cor_fast(as.matrix(x), nCores))
            }
            x_copy <- x
            return(cor_fast(x_copy, nCores))
        }
        else{
            if(isXSparse){
                return(wCor_fast(as.matrix(x), w, nCores))
            }
            x_copy <- x
            return(wCor_fast(x_copy, w, nCores))
        }    
    }
    else{
        stop("Selected distance method not available.");
    }
}

#' Perform Truncated Singular Value Decomposition (SVD)
#'
#' Computes a truncated SVD on a matrix using the implicitly restarted 
#' Lanczos bidiagonalization algorithm (IRLBA).
#'
#' @inheritParams ccImpute
#' @param x A numeric matrix to perform SVD on.
#' @param nCores The number of cores to use for parallel processing. 
#' 
#' @return A matrix containing the right singular vectors of `x`.
#' 
#' @details
#' This function utilizes the `irlba` function from the `irlba` 
#' package to efficiently calculate the truncated SVD of the input matrix `x`. 
#' The returned matrix contains `nv` right singular vectors, which are often
#' used for dimensionality reduction and feature extraction in various
#' applications.
#'
#' @examples
#' library(scater)
#' library(splatter)
#' 
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' cores <- 2
#' logX <- as.matrix(logcounts(sce))
#' w <- rowVars_fast(logX, cores)
#' corMat <- getCorM("spearman", logcounts(sce), w, cores)
#' v <- doSVD(corMat, nCores=cores)
#'
#' @importFrom irlba irlba
#' @export
doSVD <- function(x, svdMaxRatio=.08, nCeil=2000, nCores){
    nv <- ceiling(svdMaxRatio*min(nCeil,ncol(x)))
    scale <- getScale(x, nCores)
    partial_svd <- (irlba::irlba(x, nv =nv,
                                    center = scale$means, scale = scale$sds))
    return(partial_svd$v)
}

#' Estimate the Number of Clusters (k) Using the Tracy-Widom Bound
#'
#' This function estimates the number of clusters (k) in a dataset using the
#' Tracy-Widom distribution as a bound for the eigenvalues of the scaled data
#' covariance matrix.
#'
#' @param x A numeric matrix or data frame where rows are observations and
#' columns are variables.
#'
#' @return The estimated number of clusters (k).
#' 
#' @examples
#' library(scater)
#' library(splatter)
#' 
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' logX <- as.matrix(logcounts(sce))
#' k <- estkTW(logX)
#' 
#' @export
estkTW <- function(x) {
    n <- nrow(x)
    p <- ncol(x)
    sqrt_term <- sqrt(n - 1) + sqrt(p)
    muTW <- sqrt_term^2
    sigmaTW <- sqrt_term * (1/sqrt(n - 1) + 1/sqrt(p))^(1/3)
    bd <- 3.273 * sigmaTW + muTW
    evals <- eigen(crossprod(scale(x)), symmetric = TRUE, 
                    only.values = TRUE)$values
    return(sum(evals > bd))
}

#' Perform Consensus K-Means Clustering
#'
#' Executes k-means clustering on multiple subsets of data defined by
#' singular value decomposition (SVD) components, and then aggregates the
#' results into a consensus matrix.
#'
#' @inheritParams ccImpute
#' @param logX A (sparse or dense) numeric matrix representing the transpose of
#' a log-normalized gene expression matrix. Rows correspond to cells, and 
#' columns correspond to genes.
#' @param v A matrix of right singular vectors obtained from SVD of a distance
#' matrix derived from `logX`.
#'
#' @return A consensus matrix summarizing the clustering results across
#' multiple sub-datasets.
#'
#' @examples
#' library(scater)
#' library(BiocParallel)
#' library(splatter)
#' 
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' cores <- 2
#' logX <- as.matrix(logcounts(sce))
#' w <- rowVars_fast(logX, cores)
#' corMat <- getCorM("spearman", logcounts(sce), w, cores)
#' v <- doSVD(corMat, nCores=cores)
#' BPPARAM = MulticoreParam(cores)
#' consMtx <- runKM(logX, v, BPPARAM=bpparam())
#' 
#' @importFrom stats kmeans
#' @importFrom BiocParallel bplapply
#' @export
runKM <- function(logX, v, maxSets = 8, k, consMin=0.75, kmNStart, kmMax=1000,
                    BPPARAM=bpparam()){
    # 1. Determine the number of cells (n) and initialize variables 
    lv <- ceiling(0.5*ncol(v))  # Lower bound of SVD vector
    rv <- ncol(v)               # Upper bound of SVD vector

    # 2. Calculate a sequence of SVD vector indices for defining subsets
    nDim <- seq(from=lv,to=rv, by=ceiling((rv-lv)/maxSets))
    if(length(nDim) > maxSets){
        nDim <- nDim[seq_len(maxSets)]
    }
    # 3. Estimate the optimal number of clusters (k) if not provided
    if(missing(k)){
        warning("For potentially better imputation results, please specify the
                number of clusters (k). Currently estimating k using the 
                Tracy-Widom bound.")
        k <- estkTW(logX)
    }

    # 4. Determine the number of k-means random starts (kmNStart)
    if(missing(kmNStart)){kmNStart <- ifelse(ncol(logX) >= 2000, 50, 1000)}

    # 5. Perform k-means clustering on each subset of data in parallel
    kmResults <- bplapply(nDim, BPPARAM = BPPARAM, FUN = 
                            function(i, v, k, kmNStart,kmMax){
                                    x <- v[, seq_len(i)]
                                    stats::kmeans(x, k, iter.max = kmMax, 
                                    nstart = kmNStart)$cluster
                            }, v = v, k = k, kmNStart=kmNStart, kmMax=kmMax)
    
    # 6. Create and return the consensus matrix from clustering results, the
    # matrix is filtered and column normalized.
    return(getConsMtx(matrix(unlist(kmResults), nrow = length(kmResults[[1]]))
            , consMin, bpnworkers(BPPARAM)))
}

#' Identify Dropout Events in Single-Cell Expression Data
#'
#' Determines which zero values within a transposed, log-normalized expression
#' matrix are likely dropout events. The identification is based on a weighted
#' cell voting scheme, where weights are derived from a processed consensus 
#' matrix.
#'
#' @param logX A (sparse or dense) numeric matrix representing the transpose of
#' a log-normalized gene expression matrix. Rows correspond to cells, and 
#' columns correspond to genes.
#' @param consMtx A numeric matrix representing the processed consensus matrix
#' obtained from clustering analysis.
#'
#' @return A two-column matrix (or data frame) where each row indicates the 
#' location (row index, column index) of a potential dropout event in the input
#' matrix `logX`.
#'
#' @examples 
#' library(scater)
#' library(BiocParallel)
#' library(splatter)
#' 
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' cores <- 2
#' logX <- as.matrix(logcounts(sce))
#' w <- rowVars_fast(logX, cores)
#' corMat <- getCorM("spearman", logcounts(sce), w, cores)
#' v <- doSVD(corMat, nCores=cores)
#' BPPARAM = MulticoreParam(cores)
#' consMtx <- runKM(logX, v, BPPARAM=bpparam())
#' dropIds <- findDropouts(logX, consMtx)
#' 
#' @importFrom Matrix which
#' @export
findDropouts <- function(logX, consMtx){
    zero_indices <- logX==0
    vote_m <- zero_indices*2-1
    votes <- (vote_m %*% consMtx)*zero_indices
    return(Matrix::which(votes<0, arr.ind = TRUE))
}

#' Impute Dropout Values in a Log-normalized Expression Count Matrix
#'
#' @description This function imputes dropout values (zeros) in a count matrix
#'   using either a fast numerical solver or a slower linear equations solver.
#'
#' @param consMtx  A numeric matrix representing the processed consensus matrix
#' obtained from clustering analysis.
#' @param logX A (sparse or dense) numeric matrix representing the transpose of
#' a log-normalized gene expression matrix. Rows correspond to cells, and 
#' columns correspond to genes.
#' @param dropIds A numeric vector containing the row/col indices of the 
#' dropouts to be imputed.
#' @param fastSolver A logical value indicating whether to use the fast solver
#'   (default) or the slow solver.
#' @param nCores An integer specifying the number of cores to use for parallel
#'   processing (if applicable).
#'
#' @return An imputed log-transformed count matrix (same dimensions as `logX`).
#'
#' @examples
#' library(scater)
#' library(BiocParallel)
#' library(splatter)
#' 
#' sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
#'         batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
#'         dropout.type = "experiment")
#' sce <- logNormCounts(sce)
#' cores <- 2
#' logX <- as.matrix(logcounts(sce))
#' w <- rowVars_fast(logX, cores)
#' corMat <- getCorM("spearman", logcounts(sce), w, cores)
#' v <- doSVD(corMat, nCores=cores)
#' BPPARAM = MulticoreParam(cores)
#' consMtx <- runKM(logX, v, BPPARAM=bpparam())
#' dropIds <- findDropouts(logX, consMtx)
#' impLogX <- computeDropouts(consMtx, logX, dropIds, nCores=cores)
#'
#' @export
computeDropouts <- function(consMtx, logX, dropIds, fastSolver=TRUE, nCores){
    isXSparse <- inherits(logX, "dgCMatrix")
    if(fastSolver){
        imp <- NULL;
        if (isXSparse)
            imp <- sparseSolver2(consMtx, logX, dropIds, nCores)
        else
            imp <- solver2(consMtx, logX, dropIds, nCores)

        iLogX <-logX
        iLogX[dropIds] <- imp
        rownames(iLogX) <- rownames(logX)
        colnames(iLogX) <- colnames(logX)
        return(iLogX)
    }
    else{
        if (isXSparse)
            stop("Slow solver not supported for sparse matrices.")
        warning("Slow solver selected: Linear Equations solver. This might 
                significantly increase processing time. Consider using the fast
                solver. For more details, refer to the documentation.")
        iLogX <- solver(consMtx, logX, dropIds, nCores)
        rownames(iLogX) <- rownames(logX)
        colnames(iLogX) <- colnames(logX)
        return(iLogX)
    }
}