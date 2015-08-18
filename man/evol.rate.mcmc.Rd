\name{evol.rate.mcmc}
\alias{evol.rate.mcmc}
\title{Bayesian MCMC method for identifying exceptional phenotypic diversification in a phylogeny}
\usage{
evol.rate.mcmc(tree, x, ngen=10000, control=list())
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{ngen}{an optional integer value indicating the number of generations for the MCMC.}
	\item{control}{a list of control parameters containing the following elements: \code{sig1}: starting value for \eqn{\sigma(1)^2}; \code{sig2}: starting value for \eqn{\sigma(2)^2}; \code{a}: starting value for a; \code{sd1}: standard deviation for the normal proposal distribution for \eqn{\sigma(1)^2}; \code{sd2}: standard deviation for the normal proposal distribution for \eqn{\sigma(2)^2}; \code{kloc}: scaling parameter for tree move proposals - \eqn{1/\lambda} for the reflected exponential distribution; \code{sdlnr}: standard deviation on the log-normal prior on \eqn{\sigma(1)^2/\sigma(2)^2}; \code{rand.shift}: probability of proposing a random shift in the tree (improves mixing); \code{print}: print frequency for the MCMC; \code{sample}: sample frequency.}
}
\description{
	This function takes a phylogenetic tree and data for a single continuously valued character and uses a Bayesian MCMC approach to identify the phylogenetic location of a shift in the evolutionary rate through time.
}
\details{
	Default values of \code{control} are given in Revell et al. (2012).
}
\value{
	A list with the following components:
	\item{mcmc}{results from the MCMC run.}
	\item{tips}{list of stips in rate \eqn{\sigma(1)^2} for each sampled generation of MCMC (to polarize the rate shift).}
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J., D. L. Mahler, P. Peres-Neto, and B. D. Redelings. (2012) A new method for identifying exceptional phenotypic diversification. \emph{Evolution}, \bold{66}, 135-146.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.Bayes}}, \code{\link{brownie.lite}}, \code{\link{evol.vcv}}, \code{\link{minSplit}}, \code{\link{posterior.evolrate}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
