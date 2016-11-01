\name{cospeciation}
\alias{cospeciation}
\alias{plot.cospeciation}
\alias{print.cospeciation}
\title{Conducts a statistical test of cospeciation between two trees}
\usage{
cospeciation(t1, t2, distance=c("RF","SPR"),
	method=c("simulation","permutation"), assoc=NULL,
	nsim=100, ...)
\method{plot}{cospeciation}(x, ...)
\method{print}{cospeciation}(x, ...)
}
\arguments{
	\item{t1}{object of class \code{"phylo"}.}
	\item{t2}{object of class \code{"phylo"}.}
	\item{distance}{distance method to compare trees.}
	\item{method}{method to use (simulation of pure-birth trees, or permutation of tip labels on a fixed tree) to obtain a null distribution of tree distances via \code{distance}.}
	\item{assoc}{matrix containing the tip labels in \code{t1} to match to the tip labels in \code{t2}. Note that not all labels in either tree need to be included; however, unlike \code{cophylo}, one label in \code{t1} cannot be matched with more than one label in \code{t2}, nor vice versa. If \code{NULL} then an exact match of tip labels will be sought.}
	\item{nsim}{number of simulations or permutations.}
	\code{x}{for \code{plot} and \code{print} methods, an object of class \code{"cospeciation"}.}
	\item{...}{optional arguments.}
}
\description{
	This function conducts a test for cospeciation based on tree distance, applying a distance metric selected by the user. Note that this method should be prone to be quite liberal as the null hypothesis is no similarity between trees!
}
\value{
	An object of class \code{"cospeciation"}, which includes the test-statistic, the null distribution, and a p-value for the test of the null hypothesis of no topological similarity between the two trees.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{cophylo}}
}
\keyword{phylogenetics}
\keyword{plotting}

