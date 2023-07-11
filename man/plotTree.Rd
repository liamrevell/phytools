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
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plot.phylo}}, \code{\link{plotSimmap}}
}
\examples{
data(vertebrate.tree)
plotTree(vertebrate.tree,fsize=1.2,ftype="i")
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
}
\keyword{phylogenetics}
\keyword{plotting}
