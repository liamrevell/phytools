\name{to.matrix}
\alias{to.matrix}
\title{Convert a character vector to a binary matrix}
\usage{
to.matrix(x, seq)
}
\arguments{
	\item{x}{a vector of characters.}
	\item{seq}{the sequence for the columns in the output matrix.}
}
\description{
	Creates a binary matrix, normally from a factor.
}
\details{
	This function takes a vector of characters or a factor and computes a binary matrix.
	
	Primarily designed to be used internally by \code{\link{make.simmap}} and \code{\link{rerootingMethod}}.
}
\value{
	A binary matrix of dimensions \code{length(x)} by \code{length(seq)}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{utilities}
