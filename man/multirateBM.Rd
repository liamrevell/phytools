\name{multirateBM}
\alias{multirateBM}
\title{Function to fit a multi-rate Brownian evolution model}
\usage{
multirateBM(tree, x, method=c("ML","REML"), 
	optim=c("Nelder-Mead","BFGS","CG"),
	maxit=2000, n.iter=1, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a named numerical vector.}
	\item{method}{method of optimization. Currently only \code{method="ML"}.}
	\item{optim}{optimization routine to be used by \code{\link{optim}}. If more than one is specified and \code{n.iter>1} then they will be alternated.}
	\item{maxit}{to be passed to \code{optim}.}
	\item{n.iter}{number of times to reiterate failed optimization.}
	\item{...}{optional arguments.}
}
\description{
	This function fits a flexible Brownian multi-rate model
}
\value{
	An object of class \code{"multirateBM"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{continuous character}
\keyword{maximum likelihood}
