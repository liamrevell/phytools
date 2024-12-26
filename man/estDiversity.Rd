\name{estDiversity}
\alias{estDiversity}
\title{Estimate diversity at each node of the tree}
\usage{
estDiversity(tree, x, method=c("asr","simulation"), model="ER", ...)
}
\arguments{
	\item{tree}{is a phylogenetic tree in \code{"phylo"} format.}
	\item{x}{a vector containing the biogeographic area for each of the tip taxa.}
	\item{method}{method for reconstructing ancestral biogeography.}
	\item{model}{model for ancestral character estimation. In theory, any model from \code{\link[ape]{ace}}; however only symmetric models permitted for \code{method="asr"}.}
	\item{...}{optional arguments. So far, this includes only \code{nsim}, the number of stochastic mappings to conduct using \code{\link{make.simmap}} for \code{method="simulation"}.}
}
\description{
	Estimates the lineage density at each node in the tree based on a biogeographic model (similar to Mahler et al. 2010).
}
\details{
	Two different methods are implemented in the current version.
	
	For \code{method="asr"} the state at the current node, and at each position along each co-extant internal edge, is computed as the marginal (empirical Bayesian) ancestral state reconstruction using the re-rooting method of Yang (2006). The lineage density is then computed as the sum of the marginal reconstructions (posterior probabilities) times the summed marginal ancestral reconstructions across co-extant edges. 
	
	For \code{method="simulation"}, stochastic character mapping is used to generate optional argument \code{nsim} stochastic maps of ancestral biogeography. Then the lineage density at each node is computed as the number of co-existing lineages with the same biogeography as the focal node, averaged across stochastic maps.
	
	The importance of this distinction may depend on the degree to which reconstructions at internal nodes are independent, which relates to the distinction between marginal and joint reconstruction (e.g., see Yang 2006).
}
\value{
	A vector containing the estimated lineage density at each node
}
\references{
	Mahler, D. L, L. J. Revell, R. E. Glor, and J. B. Losos. (2010) Ecological opportunity and the rate of morphological evolution in the diversification of Greater Antillean anoles. \emph{Evolution}, \bold{64}, 2731-2745.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Yang, Z. (2006) \emph{Computational Molecular Evolution}. Oxford University Press.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{fitDiversityModel}}
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{maximum likelihood}
\keyword{continuous character}
\keyword{diversification}
