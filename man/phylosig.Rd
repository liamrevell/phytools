\name{phylosig}
\alias{phylosig}
\alias{plot.phylosig}
\title{Compute phylogenetic signal with two methods}
\usage{
phylosig(tree, x, method="K", test=FALSE, nsim=1000, se=NULL, start=NULL,
   control=list(), niter=10)
\method{plot}{phylosig}(x, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{x}{vector containing values for a single continuously distributed trait. In the case of the \code{plot} method, \code{x} is an object of class \code{"phylosig"}.}
	\item{method}{method to compute signal: can be \code{"K"} or \code{"lambda"}.}
	\item{test}{logical indicating whether or not to conduct a hypothesis test of \code{"K"} or \code{"lambda"}.}
	\item{nsim}{for \code{method="K"}, number of simulations in randomization test.}
	\item{se}{named vector containing the standard errors for each species.}
	\item{start}{vector of starting values for optimization of (respectively) \eqn{\sigma^2} and \eqn{\lambda}. Only used in \code{method="lambda"} and \code{se!=NULL}.}
	\item{control}{list of control parameters for multidimensional optimization, implemented in \code{\link{optim}}. Only used in \code{method="lambda"} and \code{se!=NULL}.}
	\item{niter}{number of iterations for likelihood optimization of \eqn{\lambda} (if \code{se!=NULL}), or the number of \emph{intervals} between 0 and the maximum possible value of \eqn{\lambda} for univariate optimization of \eqn{\lambda} (if \code{se==NULL}).}
	\item{...}{optional arguments for \code{plot} method.}
}
\description{
	Calculate phylogenetic signal using two different methods (Pagel, 1999; Blomberg et al., 2003).
}
\details{
	This function computes phylogenetic signal using two different methods. It can also conduct the hypothesis tests for significant phylogenetic signal, and estimate phylogenetic signal incorporating sampling error following Ives et al. (2007).

	\eqn{\lambda} optimization is performed using \code{\link{optimize}} with the range of \eqn{\lambda} set between 0 and the theoretical upper limit of \eqn{\lambda} (determined by the relative height of the most recent internal node on the tree).
	
	\code{plot.phylosig} creates either a plot of the null distribution of \emph{K} or a likelihood surface, depending on the value of \code{method}.
}
\value{
	The function returns an object of class \code{"phylosig"}. With default arguments (\code{method="K"}, \code{test=FALSE}, and \code{se=NULL}), this will be a single numeric value. Otherwise, if \code{(method="K")}, it will consist of a list with up to the following elements:
	\item{K}{value of the K-statistic.}
	\item{sig2}{rate of evolution, \eqn{\sigma^2}, for estimation with sampling error.}
	\item{logL}{log-likelihood, for estimation with sampling error.}
	\item{P}{optionally, the P-value from the randomization test.}
	If \code{(method="lambda")}, it will be a list with up to the following elements:
	\item{lambda}{fitted value of \eqn{\lambda}.}
	\item{sig2}{rate of evolution, \eqn{\sigma^2}, for estimation with sampling error.}
	\item{logL}{log-likelihood.}
	\item{logL0}{log-likelihood for \eqn{\lambda=0}.}
	\item{P}{P-value of the likelihood ratio test.}
	\item{convergence}{value for convergence, for estimation with sampling error only. (See \code{\link{optim}}).}
	\item{message}{message from \code{\link{optim}}, for estimation with sampling error only.}
}
\references{
	Blomberg, S. P., T. Garland Jr., and A. R. Ives (2003) Testing for phylogenetic signal in comparative data: Behavioral traits are more labile. \emph{Evolution}, \bold{57}, 717-745.

	Ives, A. R., P. E. Midford, and T. Garland Jr. (2007) Within-species variation and measurement error in phylogenetic comparative biology. \emph{Systematic Biology}, \bold{56}, 252-270.

	Pagel, M. (1999) Inferring the historical patterns of biological evolution. \emph{Nature}, \bold{401}, 877-884.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## extract characters of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
    rownames(mammal.data)))
ln.homeRange<-log(setNames(mammal.data$homeRange,
    rownames(mammal.data)))
## compute phylogenetic signal K
K.bodyMass<-phylosig(mammal.tree,ln.bodyMass,
    test=TRUE)
print(K.bodyMass)
plot(K.bodyMass)
K.homeRange<-phylosig(mammal.tree,ln.homeRange,
    test=TRUE)
print(K.homeRange)
plot(K.homeRange)
## compute phylogenetic signal lambda
lambda.bodyMass<-phylosig(mammal.tree,ln.bodyMass,
    method="lambda",test=TRUE)
print(lambda.bodyMass)
plot(lambda.bodyMass)
lambda.homeRange<-phylosig(mammal.tree,ln.homeRange,
    method="lambda",test=TRUE)
print(lambda.homeRange)
plot(lambda.homeRange)
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{simulation}
\keyword{continuous character}
