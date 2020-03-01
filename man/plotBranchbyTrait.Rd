\name{plotBranchbyTrait}
\alias{plotBranchbyTrait}
\title{Plot branch colors by a quantitative trait or value}
\usage{
plotBranchbyTrait(tree, x, mode=c("edges","tips","nodes"), palette="rainbow", 
   legend=TRUE, xlims=NULL, ...)
}
\arguments{
	\item{tree}{an object of class "phylo".}
	\item{x}{either a vector of states for the edges, tips, or nodes of the tree (for \code{mode="edges"}, \code{"tips"}, and \code{"nodes"}, respectively).}
	\item{mode}{string indicating plotting mode. \code{mode="edges"}, the default, requires that the mapping state of each edge in the tree should be provided. \code{mode="tips"} takes the tip values and estimates the state at each internal node. The mapped character value along each branch is the average of the nodes subtending that branch. \code{mode="nodes"} similar to \code{"tips"}, except that the node values are provided instead of estimated.}
	\item{palette}{color palette to translate character values to color. Options are presently \code{"rainbow"} (the default), \code{"heat.colors"}, and \code{"gray"}. \code{palette} can also be a function produced by \code{\link{colorRampPalette}}.}
	\item{legend}{can be a logical value (\code{TRUE} or \code{FALSE}) or a numeric value greater than 0. In the latter case the numeric value gives the length of the plotted legend, which also acts as a scale bar for the branch lengths of the tree.}
	\item{xlims}{range for the translation map between trait values and the color map. Should be inclusive of all the values in \code{x}.}
	\item{...}{other optional arguments to be passed to \code{\link{plot.phylo}} - pretty much all arguments are available. In addition, there \code{plotBranchbyTrait} has the following additional optional arguments: \code{tol} a small tolerance value to be added to the range of \code{x}; \code{prompt} for \code{legend=TRUE}, a logical value indicating whether to prompt for the position of the legend (or not) - the default is to put the legend in the lower left hand size of the plot; \code{title} for \code{legend=TRUE}, the title of the legend; and \code{digits} for \code{legend=TRUE}, the number of digits in the quantitative scale of the legend. Finally, \code{cex} can be supplied as either a single numeric value, or as a vector of two different values. If the latter is true than the second element of \code{cex} will be passed internally to the function \code{\link{add.color.bar}}.}
}
\description{
	Function plots a tree with branches colored by the value for a quantitative trait or probability, by various methods. Unlike most other tree plotting functions in phytools, this function calls \code{\link{plot.phylo}} (not \code{plotSimmap}) internally.
}
\details{
	Note that if \code{prompt=TRUE}, the function will prompt for the position of the legend.
}
\value{
	Plots a phylogeny.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
