\name{collapseTree}
\alias{collapseTree}
\title{Interactive tree visualizer}
\usage{
collapseTree(tree, ...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
 	\item{...}{optional arguments. These \emph{mostly} match the arguments of \code{\link{plotSimmap}}, but also include the argument \code{drop.extinct=TRUE} which will (if the input tree is ultrametric) drop any 'extinct' lineages from the tree that is returned by the function.}
}
\description{
	Function creates an interactive visualization of collapsing & expanding clades on the tree.
}
\details{
	Function first plots a fan style tree, and then the user collapses node on the tree by clicking on them. Collapsed nodes are collapsed to the common ancestor of the clade. Nodes that have been collapsed can also be expanded by clicking. Right-click (or select Finish button in \pkg{RStudio}) to end.
	
	Note, for the animation to work probably when run in the \pkg{RStudio} IDE the "zoom" level should be set to 100\%.
}
\value{
	Returns the final plotted tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotTree}}, \code{\link{plotSimmap}} 
}
\examples{
\dontrun{
## load phylogeny
data(liolaemid.tree)

## run collapseTree
pruned<-collapseTree(liolaemid.tree,fsize=0.5)

## reset par to default
par(mar=c(5.1,4.1,4.1,2.1))}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{utilities}
