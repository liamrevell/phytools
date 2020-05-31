\name{evol.vcv}
\alias{evol.vcv}
\title{Likelihood test for variation in the evolutionary variance-covariance matrix}
\usage{
evol.vcv(tree, X, maxit=2000, vars=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}. If \code{tree} is an object of class \code{"phylo"} then a simple multivariate Brownian motion model will be fit to the data in \code{X}.}
	\item{X}{an \code{n} x \code{m} matrix of tip values for \code{m} continuously valued traits in \code{n} species - row names should be species names.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization. This quantity may need to be increased for difficult optimizations.}
	\item{vars}{an optional logical value indicating whether or not to estimate the variances of the parameter estimates from the Hessian matrix.}
	\item{...}{optional arguments.}
}
\description{
	This function takes an object of class \code{"simmap"} with a mapped binary or multistate trait and data for an arbitrary number of continuously valued character. It then fits the multiple evolutionary variance-covariance matrix (rate matrix) model of Revell & Collar (2009; \emph{Evolution}).
}
\details{
	This function performs optimization by maximizing the likelihood with respect to the Cholesky matrices using \code{\link{optim}}.	Optimization is by \code{method="Nelder-Mead"}. Using box constraints does not make sense here as they would be applied to the Cholesky matrix rather than the target parameters. Users may have to increase \code{maxit} for large trees and/or more than two traits.
}
\value{
	An object of class \code{"evol.vcv"} with the following components:
	\item{R.single}{vcv matrix for the single rate matrix model.}
	\item{vars.single}{optionally, a matrix containing the variances of the elements of \code{R.single}.}
	\item{logL1}{log-likelihood for single matrix model.}
	\item{k1}{number of parameters in the single marix model.}
	\item{R.multiple}{\code{m} x \code{m} x \code{p} array containing the \code{p} estimated vcv matrices for the \code{p} regimes painted on the tree.}
	\item{vars.multiple}{optionally, an array containing the variances of the parameter estimates in \code{R.multiple}.}
	\item{logL.multiple}{log-likelihood of the multi-matrix model.}
	\item{k2}{number of parameters estimated in this model.}
	\item{P.chisq}{P-value of the \eqn{\chi^2} test on the likelihood ratio.}
	\item{convergence}{logical value indicating whether or not the optimization has converged.}
}
\references{
	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evol.rate.mcmc}}, \code{\link{brownie.lite}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
