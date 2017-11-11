\name{anc.Bayes}
\alias{anc.Bayes}
\title{Bayesian ancestral character estimation}
\usage{
anc.Bayes(tree, x, ngen=10000, control=list())
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{ngen}{a integer indicating the number of generations for the MCMC.}
	\item{control}{a list of control parameters containing the following elements: \code{sig2}: starting value for \eqn{\sigma^2} (BM rate); \code{a}: starting for the state at the root node; \code{y}: starting values for the states at all internal nodes excluding the root (should be labeled with node numbers); \code{pr.mean}: means for the prior distributions in the following order - \code{sig2}, \code{a}, \code{y}, note that the prior probability distribution is exponential for \code{sig2} and normal for \code{a} and \code{y}; \code{pr.var}: variances on the prior distributions, same order as \code{pr.mean} (but the variance is not used for \code{sig2}); \code{prop}: variances on the normal proposal distributions in the same order as \code{pr.mean}; \code{sample}: sample frequency from the MCMC.}
}
\description{
	This function uses Bayesian MCMC to sample from the posterior distribution for the states at internal nodes in the tree.
}
\value{
	An object of class \code{"anc.Bayes"} that includes at least two components: 
	\item{mcmc}{a data frame with rows \code{ngen/sample+1} containing the posterior sample and likelihoods.  Matrix columns are labeled either \code{sig2} or by the node number of the internal node.}
	\item{tree}{our input phylogeny.}
}
\details{
	The \code{print} method also returns (invisibly) a vector of estimated ancestral states based on a user-supplied burn-in (or 20% of the number of generations of MCMC, if no burn-in is provided).
	
	The object class \code{plot} method by default plots a likelihood profile from the MCMC.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{anc.ML}}, \code{\link{anc.trend}}, \code{\link{evol.rate.mcmc}}, \code{\link{fastAnc}}
}
\examples{
tree<-pbtree(n=50)
x<-fastBM(tree,sig2=2) # simulate using fastBM
obj<-anc.Bayes(tree,x,ngen=10000) # sample ancestral states
print(obj,printlen=20) ## estimates
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{bayesian}
