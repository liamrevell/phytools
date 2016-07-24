\name{plotTree.wBars}
\alias{plotTree.wBars}
\alias{plotTree.barplot}
\title{Plot a tree with bars at the tips}
\usage{
plotTree.wBars(tree, x, scale=1, width=NULL, type="phylogram", method="plotTree", 
	tip.labels=FALSE, col="grey", border=NULL, ...)
plotTree.barplot<-function(tree,x,args.plotTree=list(),args.barplot=list())
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
	\item{args.plotTree}{in \code{plotTree.barplot}, arguments to be passed to \code{\link{plotTree}}.}
	\item{args.barplot}{in \code{plotTree.barplot}, arguments to be passed to \code{\link{barplot}}.}
	\item{...}{optional arguments to be passed to \code{\link{plotTree}} or \code{\link{plotSimmap}}.}
}
\description{
	\code{plotTree.wbars} plots a phylogeny in phylogram or fan style with bars at the tips representing the values for a phenotypic trait.
	
	\code{plotTree.barplot} creates a split plot in which a right-facing phylogram is on the left, and a bar plot is shown on the right.
}
\value{
	Plots a tree with an associated bar plot for a continuously valued character at the tips.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{barplot}}, \code{\link{dotTree}}, \code{\link{plotSimmap}}, \code{\link{plotTree}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
