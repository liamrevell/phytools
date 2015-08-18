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
\description{
	This function plots a rooted phylogram. Arguments in \code{...} are passed to \code{\link{plotSimmap}}, with the exception of optional argument \code{color} which is used to determine the plotted color of the branch lengths of the tree.
}
\value{
	This function plots a rooted phylogram.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
tree<-pbtree(n=25)
plotTree(tree,color="blue",ftype="i")
}
\keyword{phylogenetics}
\keyword{plotting}
