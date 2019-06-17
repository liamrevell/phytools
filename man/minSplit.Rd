\name{minSplit}
\alias{minSplit}
\title{Finding the minimum (median) split in the posterior sample}
\usage{
minSplit(tree, split.list, method="sum", printD=FALSE)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{split.list}{either a matrix with two named columns, \code{"node"} and \code{"bp"}; a \code{$mcmc} matrix from \code{evol.rate.mcmc()}; or the entire raw output from \code{evol.rate.mcmc()}.}
	\item{method}{an optional string indicating the criterion to minimize: options are \code{"sum"} and \code{"sumsq"}.}
	\item{printD}{logical specifying whether to print distances to screen (\code{FALSE} by default).}
}
\description{
	This function takes a phylogenetic tree and a list of splits and identifies the split with the smallest summed or summed squared distances to all the other splits.	
}
\value{
	A list with the following components:
	\item{node}{node for the minimum split.}
	\item{bp}{location on the branch leading to \code{node} of the minimum split.}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J., D. L. Mahler, P. Peres-Neto, and B. D. Redelings (2012) A new method for identifying exceptional phenotypic diversification. \emph{Evolution}, \bold{66}, 135-146.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evol.rate.mcmc}}, \code{\link{posterior.evolrate}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
