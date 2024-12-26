\name{fastAnc}
\alias{fastAnc}
\title{(Reasonably) fast estimation of ML ancestral states}
\usage{
fastAnc(tree, x, vars=FALSE, CI=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{vars}{a logical value indicating whether or not to compute variances on the ancestral state estimates. Variances are based on Equation (6) of Rohlf (2001).}
	\item{CI}{a logical value indicating whether or not to compute 95\% confidence intervals on state estimates.}
	\item{...}{optional arguments. Presently this consists of \code{anc.states}, a named vector containing ancestral states to fix. Names should correspond to node numbers in the input tree.}
}
\description{
	Estimates ancestral states for a continuous character under maximum likelihood.
}
\details{
	This function performs (reasonably) fast estimation of the ML ancestral states for a continuous trait by taking advantage of the fact that the state computed for the root node of the tree during Felsenstein's (1985) contrasts algorithm is also the MLE of the root node. Thus, the function re-roots the tree at all internal nodes and computes the contrasts state at the root each time.
	
	The function can also (optionally) compute variances or 95\% confidence intervals on the estimates.
}
\value{
	An object of class \code{"fastAnc"} consisting of either: a named vector containing the states at internal nodes - names are node numbers; or a list containing ancestral state estimates (\code{ace}), variances on the estimates (\code{var}), and/or 95\% confidence intervals (\code{CI95}).
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{ace}}, \code{\link{anc.Bayes}}, \code{\link{anc.ML}}, \code{\link{ancr}} \code{\link[ape]{pic}}
}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## extract character of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
    rownames(mammal.data)))
## estimate ancestral body sizes
fit.BM<-fastAnc(mammal.tree,ln.bodyMass,CI=TRUE)
print(fit.BM,printlen=10)
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
