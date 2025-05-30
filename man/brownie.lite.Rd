\name{brownie.lite}
\alias{brownie.lite}
\title{Likelihood test for rate variation in a continuous trait}
\usage{
brownie.lite(tree, x, maxit=2000, test="chisq", nsim=100, se=NULL, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree either as an object of class \code{"phylo"} or \code{"simmap"}. (See \code{\link{read.simmap}}, \code{\link{make.simmap}}, or \code{\link{paintSubTree}} for more details about the latter object class.)}
	\item{x}{a vector of tip values for species. \code{names(x)} should be the species names.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization - may need to be increased for large trees.}
	\item{test}{an optional string indicating the method for hypothesis testing - options are \code{"chisq"} or \code{"simulation"}.}
	\item{nsim}{number of simulations (only used if \code{test="simulation"}).}
	\item{se}{a vector containing the standard errors for each estimated mean in \code{x}.}
	\item{...}{optional arguments.}
}
\description{
	Fits a multi-rate Brownian motion evolution model using maximum likelihood.
}
\details{
	This function takes an object of class \code{"phylo"} or class \code{"simmap"} with a mapped binary or multi-state trait (see \code{\link{read.simmap}}) and data for a single continuously valued character.  It then fits the Brownian rate variation ("noncensored") model of O'Meara et al. (2006; \emph{Evolution}). This is also the basic model implemented in Brian O'Meara's \emph{Brownie} software.
	
	Sampling error in the estimation of species means can also be accounted for by assigning the vector \code{se} with the species specific sampling errors for \code{x}.
}
\value{
	An object of class \code{"brownie.lite"} containing the following elements:
	\item{sig2.single}{is the rate, \eqn{\sigma^2}, for a single-rate model. This is usually the "null" model.}
	\item{a.single}{is the estimated state at the root node for the single rate model.}
	\item{var.single}{variance on the single rate estimator - obtained from the Hessian.}
	\item{logL1}{log-likelihood of the single-rate model.}
	\item{k1}{number of parameters in the single rate model (always 2).}
	\item{sig2.multiple}{is a length \emph{p} (for \emph{p} rates) vector of BM rates (\eqn{\sigma_1^2}, \eqn{\sigma_2^2}, and so on) from the multi-rate model.}
	\item{a.multiple}{is the estimated state at the root node for the multi-rate model.}
	\item{var.multiple}{\emph{p} x \emph{p} variance-covariance matrix for the \emph{p} rates - the square-roots of the diagonals should give the standard error for each rate.}
	\item{logL.multiple}{log-likelihood of the multi-rate model.}
	\item{k2}{number of parameters in the multi-rate model (\emph{p}+1).}
	\item{P.chisq}{P-value for a likelihood ratio test against the \eqn{\chi^2} distribution; or}
	\item{P.sim}{P-value for a likelihood ratio test against a simulated null distribution.}
	\item{convergence}{logical value indicating if the likelihood optimization converged.}
}
\references{
	O'Meara, B. C., C. Ane, M. J. Sanderson, and P. C. Wainwright (2006) Testing for different rates of continuous trait evolution using likelihood. \emph{Evolution}, \bold{60}, 922-933.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownieREML}}, \code{\link{evol.vcv}}, \code{\link{ratebytree}}
}
\examples{
## load data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)

## extract character of interest
buccal.length<-setNames(sunfish.data$buccal.length,
    rownames(sunfish.data))

## fit model
multiBM.fit<-brownie.lite(sunfish.tree,
    buccal.length)

## print fitted model
print(multiBM.fit)
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
