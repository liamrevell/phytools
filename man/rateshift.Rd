\name{rateshift}
\alias{rateshift}
\alias{plot.rateshift}
\title{Find the temporal position of one or more rate shifts}
\usage{
rateshift(tree, x, nrates=1, niter=10, method="ML", ...)
\method{plot}{rateshift}(x, ...)
}
\arguments{
 	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{vector of phenotypic trait values for species. \code{names(x)} should contain the species names and match \code{tree$tip.label}.}
	\item{nrates}{number of rates.}
	\item{niter}{number of iterations of optimization routine to ensure convergence.}
	\item{method}{optimization method. Can be \code{"ML"} (maximum likelihood) or \code{"REML"} (restricted maximum likelihood).}
	\item{...}{optional arguments. In the case of the \code{plot} method, these will be passed to \code{\link{plotSimmap}}. For \code{rateshift}, optional arguments include: \code{tol}, tolerance; \code{plot} & \code{print}, logical values indicating whether to plot or print the progress of the optimization (default to \code{FALSE}); \code{quiet}, logical argument indicating whether to suppress all notifications (defaults to \code{FALSE}); \code{minL}, numeric value; and \code{fixed.shift}, either a vector of fixed shift points, or a logical value.}
}
\description{
	Function finds the location of one or more rate shifts.
}
\value{
	An object of class \code{"rateshift"}.
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
