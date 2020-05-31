\name{phylo.heatmap}
\alias{phylo.heatmap}
\title{Creates a phylogenetic heat map}
\usage{
phylo.heatmap(tree, X, fsize=1, colors=NULL, standardize=FALSE, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a matrix containing data for multiple continuous characters in which \code{rownames} correspond to the tip labels of the tree.}
	\item{fsize}{an integer or vector of length 3 containing the font size for the tip labels, the trait labels, and the legend text. (If a single integer is supplied, then the value will be recycled.)}
	\item{colors}{a vector of colors to be passed to \code{\link{image}}. Can be a function call (e.g., \code{heat.colors(n=} \code{200)[200:1]}).}
	\item{standardize}{a logical value indicating whether or not to standardize each column of \code{X} to have the same variance & mean prior to analysis.}
	\item{...}{optional arguments. So far these include: \code{legend}, a logical value indicating whether or not to plot a figure legend (defaults to \code{legend=TRUE}); \code{labels}, a logical value indicating whether or not to plot trait labels (defaults to \code{labels=TRUE}); \code{split}, a numeric vector indicating the fraction of the horizontal dimension to use for the tree & heatmap, respectively (defaults to \code{split=c(0.5,0.5)}); \code{xlim}, \code{ylim}, & \code{mar}, defined as in \code{\link{par}}; and \code{ftype}, \code{lwd}, and \code{pts} as defined in \code{\link{plotSimmap}}.}
}
\description{
	Functions creates a multivariate phylogenetic \code{\link{heatmap}}.
}
\value{
	Function creates a plot.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
# simulate tree
tree <- pbtree(n=20,scale=1)
# simulate continuous character
X <- fastBM(tree, nsim=5)
phylo.heatmap(tree, X, grid=TRUE)
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{continuous character}
