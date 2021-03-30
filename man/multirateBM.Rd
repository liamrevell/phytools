\name{multirateBM}
\alias{multirateBM}
\title{Function to fit a multi-rate Brownian evolution model}
\usage{
multirateBM(tree, x, method=c("ML","REML"), 
	optim=c("Nelder-Mead","BFGS","CG"),
	maxit=NULL, n.iter=1, sig2BM=1, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a named numerical vector.}
	\item{method}{method of optimization. Currently only \code{method="ML"}.}
	\item{optim}{optimization routine to be used by \code{\link{optim}}. If more than one is specified and \code{n.iter>1} then they will be alternated.}
	\item{maxit}{to be passed to \code{optim}. If set to \code{maxit=NULL}, the default value of \code{maxit} will be used, depending on the optimization method.}
	\item{n.iter}{number of times to reiterate failed optimization.}
	\item{sig2BM}{fixed value of the instantaneous variance of the Brownian evolutionary process for the log-rates.}
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
