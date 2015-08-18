\name{lambda.transform}
\alias{lambda.transform}
\title{Lambda transformation of matrix}
\usage{
lambda.transform(lambda, C)
}
\arguments{
	\item{lambda}{scalar, usually (but not necessarily) on the interval 0,1.}
	\item{C}{matrix probably returned by \code{\link{vcv.phylo}}.}
}
\description{
	Function multiplies the off-diagonals of a square matrix by \code{lambda}. Used internally in \code{\link{phyl.pca}} and other functions.
}
\details{
	Do not use unless you know what you're doing.
}
\value{
	A matrix.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{math}
\keyword{comparative method}
