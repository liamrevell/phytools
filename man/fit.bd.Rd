\name{fit.bd}
\alias{fit.bd}
\alias{fit.yule}
\alias{lik.bd}
\alias{print.fit.bd}
\title{Fits birth-death (speciation/extinction) model to reconstructed phylogeny}
\usage{
fit.bd(tree, b=NULL, d=NULL, rho=1, ...)
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
	This function fits a birth-death model to a phylogenetic tree with edge lengths.
}
\details{
	This function duplicates \code{\link{birthdeath}} in ape and \code{make.bd} in the diversitree package.
}
\value{
	\code{fit.bd} generates an object of class \code{"fit.bd"} which can be printed.
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
\keyword{comparative methods}
\keyword{diversification}
\keyword{phylogenetics}
