\name{plotTree.wBars}
\alias{plotTree.wBars}
\title{Plot a tree with bars at the tips}
\usage{
plotTree.wBars(tree, x, scale=1, width=NULL, type="phylogram", method="plotTree", 
	tip.labels=FALSE, col="grey", border=NULL, ...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{x}{a named vector of trait values (normally > 0).}
 	\item{scale}{scaling factor for the tip bars (relative to the total tree height).}
	\item{width}{width of the tip bars.}
	\item{type}{plot type. Can be \code{"phylogram"} or \code{"fan"}.}
	\item{method}{plotting method to use. Can be \code{"plotTree"} (for \code{\link{plotTree}}) or \code{"plotSimmap"} (for \code{\link{plotSimmap}}).}
	\item{tip.labels}{argument indicating whether or not tip labels should be plotted. Defaults to \code{tip.labels=FALSE}.}
	\item{col}{colors of the plotted bars. Can be a single value or a vector with length equal to the number of tips in the tree.}
	\item{border}{single value specifying the color of the border for the plotted bars. Defaults to \code{border=NULL}, which means that black borders will be plotted.}
	\item{...}{optional arguments to be passed to code{\link{plotTree}} or \code{\link{plotSimmap}}.}
}
\description{
	Function plots a phylogeny in phylogram or fan style with bars at the tips representing the values for a phenotypic trait.
}
\value{
	Plots a tree.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}, \code{\link{plotTree}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
