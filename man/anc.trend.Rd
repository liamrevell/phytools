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
	This function estimates the evolutionary parameters and ancestral states for Brownian evolution with directional trend.
}
\details{
	Note that this will generally only work and produce sensible results for a phylogeny with some non-contemporary tips (i.e., a tree with some fossil species). The function uses \code{\link{optim}} with \code{method=} \code{"L-BFGS-B"}; however optimization is only constrained for the \code{sig2} which must be >0.
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
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{anc.Bayes}}, \code{\link{anc.ML}}, \code{\link{optim}}
}
\examples{
## simulate tree & data using fastBM with a trend (m!=0)
tree<-rtree(n=26,tip.label=LETTERS)
x<-fastBM(tree,mu=2) 
## fit no trend model
fit.bm<-anc.ML(tree,x,model="BM")
print(fit.bm)
 ## fit trend model
fit.trend<-anc.trend(tree,x)
print(fit.trend)
## compare trend vs. no-trend models & estimates
AIC(fit.bm,fit.trend)
xlim<-ylim<-range(c(fit.bm$ace,fit.trend$ace))
plot(fit.bm$ace,fit.trend$ace,pch=19,
    col=make.transparent("blue",0.5),
    xlab="no-trend estimates",
    ylab="trend estimates",
    main=paste("Comparison of ancestral state estimates",
    "\nunder no-trend vs. trend models"),font.main=3,
    cex.main=1.2,bty="l",cex=1.5,
    xlim=xlim,ylim=ylim)
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
