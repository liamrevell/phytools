\name{fitMk}
\alias{fitMk}
\title{Fits Mk model}
\usage{
fitMk(tree, x, model="SYM", fixedQ=NULL, method="fitMk", ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{model}{model. See \code{make.simmap} or \code{ace} for details.}
	\item{fixedQ}{fixed value of transition matrix \code{Q}, if one is desired.}
	\item{...}{optional arguments, including \code{pi}, the prior distribution at the root node (defaults to \code{pi="equal"}).}
}
\description{
	This function fits a so-called M\emph{k} model for discrete character evolution. It is primarily designed to be used inside of \code{make.simmap}.
}
\value{
	An object of class \code{"fitMk"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{make.simmap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
