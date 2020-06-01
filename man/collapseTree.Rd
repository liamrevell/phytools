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
	Function first plots a fan style tree, and then the user collapses node on the tree by clicking on them. Collapsed nodes are collapsed to the common ancestor of the clade. Nodes that have been collapsed can also be expanded by clicking. Right-click to end.
}
\value{
	Returns the final plotted tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotTree}}, \code{\link{plotSimmap}} 
}
\examples{
\dontrun{
data(anoletree)
pruned<-collapseTree(anoletree)}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{utilities}
