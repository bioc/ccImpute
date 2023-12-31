% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{findNDim}
\alias{findNDim}
\title{Establish what subsets of loadings from PCA distance measure are used for
for measuring cluster instability}
\usage{
findNDim(n, distPCA, pcaMin, pcaMax)
}
\arguments{
\item{n}{number of samples}

\item{distPCA}{PCA reduced distance matrix}

\item{pcaMin}{This is used to establish the number of minimum PCA features
used for generating subsets. For small datasets up to \code{500} cells this
equals pcaMin*n minimum features, where \code{n} is number of cells. For
large datasets, this corresponds to the feature count that has proportion of
variance less than \code{pcaMin}. Both pcaMin and pcaMax must be specified
to be considered.}

\item{pcaMax}{This is used to establish the number of maximum PCA features
used for generating subsets. For small datasets up to \code{500} cells this
equals pcaMax*n maximum features, where \code{n} is number of cells. For
large datasets, this corresponds to the feature count that has proportion of
variance less than \code{pcaMax}. Both pcaMin and pcaMax must be specified
to be considered.}
}
\value{
list of numbers with each number corresponding to the number of
loadings to use for clustering.
}
\description{
Establish what subsets of loadings from PCA distance measure are used for
for measuring cluster instability
}
\keyword{internal}
