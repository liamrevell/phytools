\name{anc.Bayes}
\alias{anc.Bayes}
\alias{plot.anc.Bayes}
\alias{density.anc.Bayes}
\title{Bayesian ancestral character estimation}
\usage{
anc.Bayes(tree, x, ngen=10000, control=list(), ...)
\method{plot}{anc.Bayes}(x, ...)
\method{density}{anc.Bayes}(x, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names. In the case of the \code{plot} and \code{density} methods, an object of class \code{"anc.Bayes"}.}
	\item{ngen}{a integer indicating the number of generations for the MCMC.}
	\item{control}{a list of control parameters containing the following elements: \code{sig2}: starting value for \eqn{\sigma^2} (BM rate); \code{a}: starting for the state at the root node; \code{y}: starting values for the states at all internal nodes excluding the root (should be labeled with node numbers); \code{pr.mean}: means for the prior distributions in the following order - \code{sig2}, \code{a}, \code{y}, note that the prior probability distribution is exponential for \code{sig2} and normal for \code{a} and \code{y}; \code{pr.var}: variances on the prior distributions, same order as \code{pr.mean} (but the variance is not used for \code{sig2}); \code{prop}: variances on the normal proposal distributions in the same order as \code{pr.mean}; \code{sample}: sample frequency from the MCMC.}
	\item{...}{optional arguments, including to be passed to \code{plot} and \code{density} methods.}
}
\description{
	This function uses Bayesian MCMC to sample from the posterior distribution for the states at internal nodes in the tree.
}
\value{
	\code{anc.Bayes} returns an object of class \code{"anc.Bayes"} including at least two components: 
	\item{mcmc}{a data frame with rows \code{ngen/sample+1} containing the posterior sample and likelihoods.  Matrix columns are labeled either \code{sig2} or by the node number of the internal node.}
	\item{tree}{our input phylogeny.}
}
\details{
	\code{plot.anc.Bayes} generates a likelihood profile plot of the MCMC by default, but can also create a profile plot for any of the sampled variables by using the optional argument \code{what}. For instance, \code{what=40} (or, equivalently, \code{what="40"}) will create a profile plot of the MCMC for node \code{40}. Additional arguments are passed to \code{\link{plot}}.
	
	\code{density.anc.Bayes} computes a posterior density from the MCMC sample. Like \code{plot.anc.Bayes} takes the optional argument \code{what}, but unlike \code{plot.anc.Bayes} computes the posterior density for the root node by default. The object computed by this function is of class \code{"density"} and can be visualized using \code{\link{plot.density}}. Burn-in (in generations) can be set using the optional argument \code{burnin}, otherwise it will be assumed to be 20% of the sample.

	The \code{print} and \code{summary} methods for this object class also return (invisibly) a vector of estimated ancestral states based on a user-supplied burn-in (or 20\% of the number of generations of MCMC, if no burn-in is provided). Burn-in can be specified with the optional argument \code{burnin}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{anc.ML}}, \code{\link{anc.trend}}, \code{\link{evol.rate.mcmc}}, \code{\link{fastAnc}}
}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## extract character of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
rownames(mammal.data)))
## run MCMC
mcmc<-anc.Bayes(mammal.tree,ln.bodyMass,
    ngen=100000)
print(mcmc,printlen=20) ## estimates
par(mfrow=c(2,1))
plot(mcmc,bty="l",main="Likelihood-profile from MCMC",
    font.main=3) ## likelihood-profile
plot(density(mcmc,what=Ntip(mammal.tree)+1,
    burnin=20000),bty="l",
    main="Posterior density for root state of log(body mass)",
    font.main=3)
par(mfrow=c(1,1)) ## reset par to default
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
\keyword{continuous character}