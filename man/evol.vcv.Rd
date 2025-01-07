\name{evol.vcv}
\alias{evol.vcv}
\title{Likelihood test for variation in the evolutionary variance-covariance matrix}
\usage{
evol.vcv(tree, X, maxit=2000, vars=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"}. If \code{tree} is an object of class \code{"phylo"} then a simple multivariate Brownian motion model will be fit to the data in \code{X}.}
	\item{X}{an \code{n} x \code{m} matrix of tip values for \code{m} continuously valued traits in \code{n} species - row names should be species names. If \code{X} is supplied as a data frame it will be coerced into a matrix without warning.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization. This quantity may need to be increased for difficult optimizations.}
	\item{vars}{an optional logical value indicating whether or not to estimate the variances of the parameter estimates from the Hessian matrix.}
	\item{...}{optional arguments. The most important optional argument at this time is \code{err_vcv} which should contain a list of matrices of sampling \emph{variances} and covariances for (and between) the means of each species. The sampling variance for the mean is just the square of the sampling error. Sampling covariances will tend to be zero (or close to zero) if error for different traits is uncorrelated, for instance, because different specimens were used to estimate the means for different traits, and non-zero otherwise.}
}
\description{
	Fits a multi-regime multivariate Brownian motion model following Revell & Collar (2009).
}
\details{
	The function takes an object of class \code{"simmap"} with a mapped binary or multi-state trait and data for an arbitrary number of continuously valued character. It then fits the multiple evolutionary variance-covariance matrix (rate matrix) model of Revell & Collar (2009; \emph{Evolution}).
	
	\code{evol.vcv} performs optimization by maximizing the likelihood with respect to the Cholesky matrices using \code{\link{optim}}.	Optimization is by \code{method="Nelder-Mead"}. Using box constraints does not make sense here as they would be applied to the Cholesky matrix rather than the target parameters. Users may have to increase \code{maxit} for large trees and/or more than two traits.
}
\value{
	An object of class \code{"evol.vcv"} with the following elements:
	\item{R.single}{vcv matrix for the single rate matrix model.}
	\item{vars.single}{optionally, a matrix containing the variances of the elements of \code{R.single}.}
	\item{logL1}{log-likelihood for single matrix model.}
	\item{k1}{number of parameters in the single matrix model.}
	\item{R.multiple}{\code{m} x \code{m} x \code{p} array containing the \code{p} estimated vcv matrices for the \code{p} regimes painted on the tree.}
	\item{vars.multiple}{optionally, an array containing the variances of the parameter estimates in \code{R.multiple}.}
	\item{logL.multiple}{log-likelihood of the multi-matrix model.}
	\item{k2}{number of parameters estimated in this model.}
	\item{P.chisq}{P-value of the \eqn{\chi^2} test on the likelihood ratio.}
	\item{convergence}{logical value indicating whether or not the optimization has converged.}
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Revell, L. J., and D. C. Collar (2009) Phylogenetic analysis of the evolutionary correlation using likelihood. \emph{Evolution}, \bold{63}, 1090-1100.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{evol.rate.mcmc}}, \code{\link{brownie.lite}}
}
\examples{
## load data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)

## fit multi-correlation model
sunfish.fit<-evol.vcv(sunfish.tree,sunfish.data[,2:3])
print(sunfish.fit)
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
