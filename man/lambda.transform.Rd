\name{lambda.transform}
\alias{lambda.transform}
\title{\eqn{\lambda} transformation of matrix}
\usage{
lambda.transform(lambda, C)
}
\arguments{
	\item{lambda}{scalar, usually (but not necessarily) on the interval 0,1.}
	\item{C}{matrix probably returned by \code{\link{vcv.phylo}}.}
}
\description{
	Internal function for \code{\link{phyl.pca}} and others.
}
\details{
	Multiplies the off-diagonals of a square matrix by \code{lambda} and returns the result.
}
\value{
	Typically an among-species phylogenetic variance covariance matrix (e.g., \code{\link{vcv.phylo}}) in which the off-diagonal elements have been multiplied by \code{lambda}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{math}
\keyword{comparative method}
\keyword{utilities}
