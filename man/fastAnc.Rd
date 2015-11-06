\name{fastAnc}
\alias{fastAnc}
\title{Fast estimation of ML ancestral states}
\usage{
fastAnc(tree, x, vars=FALSE, CI=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a vector of tip values for species; \code{names(x)} should be the species names.}
	\item{vars}{a logical value indicating whether or not to compute variances on the ancestral state estimates. Variances are based on Equation (6) of Rohlf (2001).}
	\item{CI}{a logical value indicating whether or not to compute 95-percent confidence intervals on state estimates.}
	\item{...}{optional arguments. Presently this consists of \code{anc.states}, a named vector containing ancestral states to fix. Names should correspond to node numbers in the input tree.}
}
\description{
	This function performs fast estimation of the ML ancestral states for a continuous trait by taking advantage of the fact that the state computed for the root node of the tree during Felsenstein's (1985) contrasts algorithm is also the MLE of the root node. Thus, the function reroots the tree at all internal nodes and computes the contrasts state at the root each time. The function can also (optionally) compute variances or 95-percent confidence intervals on the estimates.
}
\value{
	A named vector containing the states at internal nodes - names are node numbers; or a list containing ancestral state estimates (\code{ace}), variances on the estimates (\code{var}), and/or 95-percent confidence intervals (\code{CI95}).
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{ace}}, \code{\link{anc.Bayes}}, \code{\link{anc.ML}}, \code{\link{pic}}
}
\examples{
tree<-pbtree(n=50)
x<-fastBM(tree) # simulate using fastBM
fastAnc(tree,x) # estimate states
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
