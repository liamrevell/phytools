\name{multirateBM}
\alias{multirateBM}
\title{Function to fit a multi-rate Brownian evolution model}
\usage{
multirateBM(tree, x, method=c("ML","REML"), 
	optim=c("L-BFGS-B","Nelder-Mead","BFGS","CG"),
	maxit=NULL, n.iter=1, lambda=1, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a named numerical vector. Names should correspond to the species names of \code{tree}.}
	\item{method}{method of optimization. Currently only \code{method="ML"}.}
	\item{optim}{optimization routine to be used by \code{\link{optim}}. If more than one is specified and \code{n.iter>1} then they will be alternated. (This is recommended to improve optimization.)}
	\item{maxit}{to be passed to \code{optim}. If set to \code{maxit=NULL}, the default value of \code{maxit} will be used, depending on the optimization method.}
	\item{n.iter}{number of times to reiterate failed optimization.}
	\item{lambda}{lambda penalty term. High values of \code{lambda} correspond to high penalty for rate heterogeneity among edges. Low values of \code{lambda} correspond to low penalty.}
	\item{...}{optional arguments.}
}
\description{
	This function fits a flexible Brownian multi-rate model using penalized likelihood.
	
	The model that is being fit is one in which the rate of Brownian motion evolution itself evolves from edge to edge in the tree under a process of geometric Brownian evolution (i.e., Brownian motion evolution on a log scale).
	
	The penalty term, \code{lambda}, determines the cost of variation in the rate of evolution from branch to branch. If lambda is \emph{high}, then the rate of evolution will vary relative little between edges (and in the limiting case converge to the single-rate MLE estimate of the rate). By contrast, if the value of \code{lambda} is set to be low, then the rate of evolution can vary from edge to edge with relatively little penalty. 
	
	Decreasing the penalty term, however, is not with cost. As \code{lambda} is decreased towards zero, estimated rates will tend to become less and less accurate.
}
\value{
	An object of class \code{"multirateBM"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
	
	Revell, L. J. (Accepted) A variable-rate quantitative trait evolution model using penalized-likelihood. \emph{PeerJ}.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{brownie.lite}}
}
\examples{
\dontrun{
## load data
data(sunfish.tree)
data(sunfish.data)
## extract character of interest
gw<-setNames(sunfish.data$gape.width,
    rownames(sunfish.data))
## run penalized-likelihood optimization
## lambda=0.1 is arbitrary
fitBM<-multirateBM(sunfish.tree,gw,lambda=0.1)
## print and plot the results
print(fitBM)
plot(fitBM,ftype="i",fsize=0.8,lwd=4)
}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{continuous character}
\keyword{maximum likelihood}
