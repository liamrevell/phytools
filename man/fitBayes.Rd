\name{fitBayes}
\alias{fitBayes}
\title{Evolutionary model fitting with intraspecific variability using Bayesian MCMC}
\usage{
fitBayes(tree, x, ngen=10000, model="BM", method="reduced", control=list())
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of phenotypic values for individuals; \code{names(x)} should contain the species names (not individual IDs).}
	\item{ngen}{a integer indicating the number of generations for the MCMC.}
	\item{model}{an evolutionary model: either \code{"BM"} or \code{"lambda"}.}
	\item{method}{a method: either \code{"reduced"} or \code{"full"}.}
	\item{control}{a list of control parameters containing the following elements: \code{sig2}: starting value for \eqn{\sigma^2} (BM rate); \code{lambda}: starting value for the \eqn{\lambda} parameter; \code{a}: starting for the state at the root node; \code{xbar}: starting values for the states at the tips; \code{intV}: starting value for the intraspecific variance (reduced method); or \code{v}: starting value for the vector of intraspecific variances for all species (full method); \code{pr.mean}: means for the prior distributions in the following order - \code{sig2}, \code{lambda} (if applicable), \code{a}, \code{xbar}, \code{intV} or \code{v} (if applicable), note that the prior probability distribution is exponential for \code{sig2} and normal for \code{a} and \code{y}; \code{pr.var}: variances on the prior distributions, same order as \code{pr.mean}.}
}
\description{
	This function uses Bayesian MCMC to sample terminal states (species means) as well as evolutionary parameters.
}
\value{
	An object of class \code{"fitBayes"} that includes a matrix (\code{mcmc}) with a number of rows \code{ngen/control$sample+1} containing the posterior sample and likelihoods.  Matrix columns are labeled by species (for species means and variances), or by the corresponding evolutionary parameter.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. and R. G. Reynolds. (2012) A new Bayesian method for fitting evolutionary models to comparative data with intraspecific variation. \emph{Evolution}, 66, 2697-2707.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.Bayes}}, \code{\link{brownie.lite}}, \code{\link{evol.rate.mcmc}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
