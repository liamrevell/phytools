\name{plotTree}
\alias{plotTree}
\title{Plots rooted phylogenetic tree}
\usage{
plotTree(tree, ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format; or multiple trees as an object of class \code{"multiPhylo"}.}
	\item{...}{optional arguments.}
}
\details{
	Plots a rooted phylogram or cladogram.
}
\description{
	Essentially a wrapper for \code{\link{plotSimmap}}. Arguments in \code{...} are passed to \code{\link{plotSimmap}}, with the exception of optional argument \code{color} which is used to determine the plotted color of the branch lengths of the tree.
}
\value{
	This function plots a rooted phylogram or cladogram.
}
\references{
	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
## load tree
data(vertebrate.tree)

## create plot
plotTree(vertebrate.tree,fsize=1.2,ftype="i")

## reset margins to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
