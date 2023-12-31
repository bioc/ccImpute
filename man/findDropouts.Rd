% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Methods.R
\name{findDropouts}
\alias{findDropouts}
\title{Establishes which zero values in \code{x} are dropout events based on
weighted cell voting with weights derived from processed consensus matrix
\code{consMtx}.}
\usage{
findDropouts(x, consMtx)
}
\arguments{
\item{x}{transpose of log normalized expression matrix}

\item{consMtx}{processed consensus matrix}
}
\value{
list of indices in x that are dropout events
}
\description{
Establishes which zero values in \code{x} are dropout events based on
weighted cell voting with weights derived from processed consensus matrix
\code{consMtx}.
}
\keyword{internal}
