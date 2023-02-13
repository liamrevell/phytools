\name{fit.bd}
\alias{fit.bd}
\alias{fit.yule}
\alias{lik.bd}
\alias{print.fit.bd}
\title{Fits birth-death (speciation/extinction) model to reconstructed phylogeny}
\usage{
fit.bd(tree, b=NULL, d=NULL, rho=1, ...)
fit.yule(tree, b=NULL, d=NULL, rho=1, ...)
lik.bd(theta, t, rho=1, N=NULL)
\method{print}{fit.bd}(x, ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{b}{birth (speciation) rate. Presently doesn't do anything as the rate cannot be fixed.}
	\item{d}{death (extinction) rate. Presently doesn't do anything as the rate cannot be fixed.}
	\item{rho}{sampling fraction.}
	\item{theta}{vector of \code{b} and \code{d} for likelihood function.}
	\item{t}{branching times for calculation of the likelihood.}
	\item{N}{number of tips in the tree.}
	\item{x}{object of class \code{"fit.bd"} for \code{print} method.}
	\item{...}{optional arguments.}
}
\description{
	Fits a birth-death (\code{fit.bd}) or pure-birth (\code{fit.yule}) model to a reconstructed phylogenetic tree with branch lengths.
}
\details{
	The function \code{fit.bd} fits a birth-death model to a phylogenetic tree with edge lengths and a (potentially) incomplete sampling fraction.
	
	The function \code{fit.yule} fits a pure-birth model with a (potentially) incomplete sampling fraction.
	
	The function \code{lik.bd} computes the likelihood of a set of birth & death rates given the set of branching times computed for a tree and a sampling fraction.
}
\value{
	\code{fit.bd} returns an object of class \code{"fit.bd"} which can be printed. This object is a list containing the fitted model parameters, likelihood, optimization conditions, a summary of the optimization, and a likelihood function.
	
	\code{fit.yule} returns an object of class \code{"fit.yule"}. This object is a list containing the fitted model parameter, likelihood, optimization conditions, a summary of the optimization, and a likelihood function.
}
\references{
	Nee, S., May, R. M. and Harvey, P. H. (1994) The reconstructed evolutionary process. \emph{Philosophical Transactions of the Royal Society of London B}, \bold{344}, 305-311.

	Stadler, T. (2012) How can we improve the accuracy of macroevolutionary rate estimates? \emph{Systematic Biology}, \bold{62}, 321-329.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{birthdeath}}
}
\examples{
data(salamanders)
## compute sampling fraction based on 55 species of Plethodon
sampling.f<-Ntip(salamanders)/55
## fit birth-death model
bd.fit<-fit.bd(salamanders,rho=sampling.f)
print(bd.fit)
## fit Yule model
yule.fit<-fit.yule(salamanders,rho=sampling.f)
print(yule.fit)
## compare b-d and yule models
anova(yule.fit,bd.fit)
## create a likelihood surface for b-d model
ngrid<-100
b<-seq(0.01,0.06,length.out=ngrid)
d<-seq(0.005,0.03,length.out=ngrid)
logL<-sapply(d,function(d,b) sapply(b,function(b,d) 
    bd.fit$lik(c(b,d)),d=d),b=b)
contour(x=b,y=d,logL,nlevels=100,
    xlab=expression(lambda),
    ylab=expression(mu),bty="l")
title(main="Likelihood surface for plethodontid diversification",
    font.main=3)
points(bd.fit$b,bd.fit$d,cex=1.5,pch=4,
    col="blue",lwd=2)
legend("bottomright","ML solution",pch=4,col="blue",
    bg="white",pt.cex=1.5,pt.lwd=2)
}
\keyword{comparative methods}
\keyword{diversification}
\keyword{phylogenetics}
\keyword{diversification}
\keyword{maximum likelihood}
