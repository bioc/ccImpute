% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{getCorM}
\alias{getCorM}
\title{Calculate Column-wise Correlation Matrix}
\usage{
getCorM(method, x, w, nCores)
}
\arguments{
\item{method}{A character string specifying the correlation 
metric to use. Currently supported options are:
 - \code{"spearman"}: Spearman rank correlation
 - \code{"pearson"}: Pearson correlation}

\item{x}{A numeric matrix where each column represents a sample.}

\item{w}{(Optional) A numeric vector of weights for each feature (row) in 
\code{x}. If not provided, all features are equally weighted.}

\item{nCores}{The number of cores to use for parallel processing.}
}
\value{
A correlation matrix of the same dimensions as the number of columns
in `x`. The values represent the pairwise correlations between samples
(columns) based on the chosen method and optional weights.
}
\description{
Efficiently computes a column-wise correlation matrix for a 
given input matrix. Supports Pearson and Spearman correlations, with
optional weighting for features.
}
\examples{
library(scater)
library(splatter)

sce <- splatSimulate(group.prob = rep(1, 5)/5, sparsify = FALSE, 
        batchCells=100, nGenes=1000, method = "groups", verbose = FALSE, 
        dropout.type = "experiment")
sce <- logNormCounts(sce)
cores <- 2
logX <- as.matrix(logcounts(sce))
w <- rowVars_fast(logX, cores)
corMat <- getCorM("spearman", logcounts(sce), w, cores)

}
