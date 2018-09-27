\name{Dtest}
\alias{Dtest}
\title{Conducts correlational D-test from stochastic mapping}
\usage{
Dtest(t1, t2, nsim=100, ...)
}
\arguments{
	\item{t1}{set of stochastic map trees (i.e., object of class \code{"multiSimmap"} for character 1. Note that \code{t1} and \code{t2} should be of the same length.}
	\item{t2}{set of stochastic map trees (i.e., object of class \code{"multiSimmap"} for character 2. Note that \code{t1} and \code{t2} should be of the same length.}
	\item{nsim}{number of simulations to use in the test.}
	\item{...}{arguments to be passed internally to \code{make.simmap}. Note that (for now) these must be the same for both \code{t1} and \code{t2} (that is to say, we are not able to assume different trait evolution models for each tree).}
}
\description{
	This function conducts the 'D-test' of Huelsenbeck et al. (2003).
}
\details{
	Note that this function has been included without much testing, and so the user should be wary.
}
\value{
	An object of class \code{"Dtest"}.
}
\references{
	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{map.overlap}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{simulation}
\keyword{bayesian}
