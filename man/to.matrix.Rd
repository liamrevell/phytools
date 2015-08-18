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
	This function takes a vector of characters and computes a binary matrix. Primarily to be used internally by \code{\link{make.simmap}} and \code{\link{rerootingMethod}}.
}
\value{
	A binary matrix of dimensions \code{length(x)} by \code{length(seq)}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{utilities}
