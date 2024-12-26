\name{anc.trend}
\alias{anc.trend}
\title{Ancestral character estimation with a trend}
\usage{
anc.trend(tree, x, maxit=2000)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{maxit}{an optional integer value indicating the maximum number of iterations for optimization.}
}
\description{
	This function estimates the evolutionary parameters and ancestral states for Brownian evolution with a directional trend.
}
\details{
	Note that this will generally only work and produce sensible results for a phylogeny with some non-contemporaneous tips (i.e., a tree with some fossil species).
	
	The function uses \code{\link{optim}} with \code{method=} \code{"L-BFGS-B"} internally; however, optimization is only constrained for the \code{sig2} which must be > 0.
}
\value{
	An object of class \code{"anc.trend"} with the following components:
	\item{ace}{a vector with the ancestral states.}
	\item{mu}{a trend parameter per unit time.}
	\item{sig2}{the variance of the BM process, \eqn{\sigma^2}.}
	\item{logL}{the log-likelihood.}
	\item{convergence}{the value of \code{$convergence} returned by \code{optim()} (0 is good).}
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{ace}}, \code{\link{anc.Bayes}}, \code{\link{anc.ML}}, \code{\link{optim}}
}
\examples{
## simulate tree & data using fastBM with a trend (m!=0)
tree<-rtree(n=26,tip.label=LETTERS)
x<-fastBM(tree,mu=4,internal=TRUE)
a<-x[as.character(1:tree$Nnode+Ntip(tree))]
x<-x[tree$tip.label]
## fit no trend model
fit.bm<-anc.ML(tree,x,model="BM")
print(fit.bm)
## fit trend model
fit.trend<-anc.trend(tree,x)
print(fit.trend)
## compare trend vs. no-trend models & estimates
AIC(fit.bm,fit.trend)
layout(matrix(c(3,4,1,2,5,6),3,2,byrow=TRUE),
    heights=c(1.5,3,1.5),widths=c(3,3))
xlim<-ylim<-range(c(a,fit.bm$ace,
    fit.trend$ace))
plot(a,fit.bm$ace,pch=19,
    col=make.transparent("blue",0.5),
    xlab="true ancestral states",
    ylab="ML estimates",
    main=paste("Comparison of true and estimated",
    "\nstates under a no-trend model"),font.main=3,
    cex.main=1.2,bty="l",cex=1.5,
    xlim=xlim,ylim=ylim)
lines(xlim,ylim,lty="dotted")
plot(a,fit.trend$ace,pch=19,
    col=make.transparent("blue",0.5),
    xlab="true ancestral states",
    ylab="ML estimates",
    main=paste("Comparison of true and estimated",
    "\nstates under a trend model"),font.main=3,
    cex.main=1.2,bty="l",cex=1.5,
    xlim=xlim,ylim=ylim)
lines(xlim,ylim,lty="dotted")
par(mfrow=c(1,1))
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
