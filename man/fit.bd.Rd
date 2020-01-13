\name{fit.bd}
\alias{fit.bd}
\alias{fit.yule}
\alias{lik.bd}
\alias{print.fit.bd}
\title{Fits birth-death (speciation/extinction) model to reconstructed phylogeny}
\usage{
fit.bd(tree, b=NULL, d=NULL, rho=1, ...)
fit.yule(tree, b=NULL, d=NULL, rho=1, ...)
lik.bd(theta, t, rho=1, N=NULL)
\method{print}{fit.bd}(x, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{b}{birth (speciation) rate. Presently doesn't do anything as the rate cannot be fixed.}
	\item{d}{death (extinction) rate. Presently doesn't do anything as the rate cannot be fixed.}
	\item{rho}{sampling fraction.}
	\item{theta}{vector of \code{b} and \code{d} for likelihood function.}
	\item{t}{branching times for calculation of the likelihood.}
	\item{N}{number of tips in the tree.}
	\item{x}{object of class \code{"fit.bd"} for \code{print} method.}
	\item{...}{optional arguments.}
}
\description{
	The function \code{fit.bd} fits a birth-death model to a phylogenetic tree with edge lengths and a (potentially) incomplete sampling fraction.
	
	The function \code{fit.yule} fits a pure-birth model with a (potentially) incomplete sampling fraction.
	
	The function \code{lik.bd} computes the likelihood of a set of birth & death rates given the set of branching times computed for a tree and a sampling fraction.
}
\value{
	\code{fit.bd} returns an object of class \code{"fit.bd"} which can be printed. This object is a list containing the fitted model parameters, likelihood, optimization conditions, a summary of the optimization, and a likelihood function.
	
	\code{fit.yule} returns an object of class \code{"fit.yule"}. This object is a list containing the fitted model parameter, likelihood, optimization conditions, a summary of the optimization, and a likelihood function.
}
\references{
	Nee, S., May, R. M. and Harvey, P. H. (1994) The reconstructed evolutionary process. \emph{Philosophical Transactions of the Royal Society of London B}, \bold{344}, 305-311.

	Stadler, T. (2012) How can we improve the accuracy of macroevolutionary rate estimates? \emph{Systematic Biology}, \bold{62}, 321-329.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{birthdeath}}
}
\examples{
\dontrun{
## fit BD & Yule models
bd<-fit.bd(tree,rho=0.8)
yule<-fit.yule(tree,rho=0.8)
## compare fitted models
AIC(yule,bd)
}
}
\keyword{comparative methods}
\keyword{diversification}
\keyword{phylogenetics}
