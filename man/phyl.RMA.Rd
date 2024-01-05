\name{phyl.RMA}
\alias{phyl.RMA}
\alias{coef.phyl.RMA}
\alias{plot.phyl.RMA}
\title{Phylogenetic reduced major axis (RMA) regression}
\usage{
phyl.RMA(x, y, tree, method="BM", lambda=NULL, fixed=FALSE, h0=1.0)
\method{coef}{phyl.RMA}(object, ...)
\method{plot}{phyl.RMA}(x, ...)
}
\arguments{
	\item{x}{vector with names. In the case of the S3 \code{plot} method \code{x} is an object of class \code{"phyl.RMA"}.}
	\item{y}{vector with names.}
	\item{tree}{a phylogenetic tree in \code{"phylo"} format.}
	\item{method}{method to obtain the correlation structure: can be \code{"BM"} or \code{"lambda"}.}
	\item{lambda}{value of \code{lambda} for fixed \eqn{\lambda}.}
	\item{fixed}{logical value indicating whether or not \eqn{\lambda} should be optimized using likelihood.}
	\item{h0}{null hypothesis for \eqn{\beta}. Defaults to 1.0. Note that a null hypothesis of 0.0 is not allowed.}
	\item{object}{for \code{coef} method, an object of class \code{"phyl.RMA"}.}
	\item{...}{optional arguments for S3 methods.}
}
\description{
	Phylogenetic reduced major axis (RMA) regression.
}
\details{
	Optionally jointly estimates \eqn{\lambda} if \code{method="lambda"}. Likelihood optimization of \eqn{\lambda} is performed using \code{\link{optimize}} on the interval (0,1).
	
	The statistical hypothesis testing is based on Clarke (1980; reviewed in McArdle 1988), which differs from some other implementations of non-phylogenetic major axis regression in R.
	
	Note that some statistician think there is \emph{never} a condition in which a reduced-major-axis regression is appropriate.
}
\value{
	An object of class \code{"phyl.RMA"} consisting of a list with the following elements:
	\item{RMA.beta}{a vector of RMA regression coefficients.}
	\item{V}{a VCV matrix for the traits.}
	\item{lambda}{fitted value of \eqn{\lambda} (\code{method="lambda"} only).}
	\item{logL}{log-likelihood (\code{method="lambda"} only).}
	\item{test}{a vector containing results for hypothesis tests on \eqn{\beta}.}
	\item{resid}{a vector of residuals for \code{y} given \code{x}.}
}
\references{
	Clarke, M. R. B. (1980) The reduced major axis of a bivariate sample. \emph{Biometrika}, \bold{67}, 441-446.
	
	McArdle, B. H. (1988) The structural relationship: Regression in biology. \emph{Can. J. Zool.}, \bold{66}, 2329-2339.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{phyl.cca}}, \code{\link{phyl.pca}}, \code{\link{phyl.resid}}
}
\examples{
## load data from Garland et al. (1992)
data(mammal.data)
data(mammal.tree)
## pull out & log transform variables
lnBodyMass<-setNames(log(mammal.data$bodyMass),
  rownames(mammal.data))
lnHomeRange<-setNames(log(mammal.data$homeRange),
  rownames(mammal.data))
## fit RMA regression & print results
fitted.rma<-phyl.RMA(lnBodyMass,lnHomeRange,
  mammal.tree)
print(fitted.rma)
## plot fitted RMA
par(las=1,bty="n")
plot(fitted.rma,las=1,bty="n")
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{statistics}
