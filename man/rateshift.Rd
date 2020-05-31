\name{rateshift}
\alias{rateshift}
\alias{plot.rateshift}
\alias{likSurface.rateshift}
\title{Find the temporal position of one or more rate shifts}
\usage{
rateshift(tree, x, nrates=1, niter=10, method="ML", ...)
\method{plot}{rateshift}(x, ...)
likSurface.rateshift(tree, x, nrates=2, shift.range=NULL,
   density=20, plot=TRUE, ...)
}
\arguments{
 	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{vector of phenotypic trait values for species. \code{names(x)} should contain the species names and match \code{tree$tip.label}. For \code{plot} method, \code{x} is an object of class \code{"rateshift"}.}
	\item{nrates}{number of rates.}
	\item{niter}{number of iterations of optimization routine to ensure convergence.}
	\item{method}{optimization method. Can be \code{"ML"} (maximum likelihood) or \code{"REML"} (restricted maximum likelihood).}
	\item{...}{optional arguments. In the case of the \code{plot} method, these will be passed to \code{\link{plotSimmap}}. For \code{rateshift}, optional arguments include: \code{tol}, tolerance; \code{plot} & \code{print}, logical values indicating whether to plot or print the progress of the optimization (default to \code{FALSE}); \code{quiet}, logical argument indicating whether to suppress all notifications (defaults to \code{FALSE}); \code{minL}, numeric value; and \code{fixed.shift}, either a vector of fixed shift points, or a logical value.}
	\item{shift.range}{for \code{likSurface.rateshift}.}
	\item{density}{for \code{likSurface.rateshift}.}
	\item{plot}{logical argument for \code{likSurface.rateshift}. If \code{plot=FALSE} then the surface is returned.}
}
\description{
	Function finds the location of one or more rate shifts. \code{likSurface.rateshift} plots the likelihood surface (which is sometimes quite rugged).
}
\value{
	A fitted object of class \code{"rateshift"}, or, in the case of \code{likSurface.rateshift}, a likelihood surface for the shift points.
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
