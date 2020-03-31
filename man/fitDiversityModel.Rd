\name{fitDiversityModel}
\alias{fitDiversityModel}
\alias{logLik.fitDiversityModel}
\alias{print.fitDiversityModel}
\title{Fit diversity-dependent phenotypic evolution model}
\usage{
fitDiversityModel(tree, x, d=NULL, showTree=TRUE, tol=1e-6)
\method{logLik}{fitDiversityModel}(object, ...)
\method{print}{fitDiversityModel}(x, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector with tip values for a continuously distributed trait. For \code{print} method, an object of class \code{"fitDiversityModel"}.}
	\item{d}{a vector containing the inferred historical diversity at each node in the tree - if \code{d=NULL} (the default) function will treat the diversification as if it occurred in a single geographic area.}
	\item{showTree}{optional logical value indicating whether to plot the tree transformation implied by the model.}
	\item{tol}{some small value by which \code{d} is incremented during rescaling of \code{psi} for optimization. If R thinks your matrices are singular during optimization, try increasing \code{tol} slightly.}
	\item{object}{for \code{logLik} method, an object of class \code{"fitDiversityModel"}.}
	\item{...}{optional arguments for \code{logLik} and \code{print} methods. Note that for the \code{logLik} method the number of fitted parameters (\code{"df"}) is assumed to be 3 for the diversity dependent model (that is, if \code{psi} is estimated) and 2 for the diversity independent model, unless otherwise specified (using the argument \code{df}).}
}
\description{
	This function fits a diversity-dependent phenotypic evolution model (based on Mahler et al. 2010).
}
\value{
	An object of class \code{"fitDiversityModel"} consisting of the following components:
	\item{logL}{log-likelihood of the fitted model.}
	\item{sig0}{estimated starting value for the rate at the root of the tree, \eqn{\sigma_0^2}.}
	\item{psi}{the estimated rate of change in the rate associated with the addition of a lineage.}
	\item{vcv}{a matrix with the variances and covariance of the estimated parameters (from the Hessian).}
}
\references{
	Mahler, D. L, L. J. Revell, R. E. Glor, and J. B. Losos. 2010. Ecological opportunity and the rate of morphological evolution in the diversification of Greater Antillean anoles. \emph{Evolution}, \bold{64}, 2731-2745.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}, \code{\link{estDiversity}}, \code{\link{evol.rate.mcmc}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
