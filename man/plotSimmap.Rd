\name{plotSimmap}
\alias{plotSimmap}
\alias{plot.simmap}
\alias{plot.multiSimmap}
\title{Plot stochastic character mapped tree}
\usage{
plotSimmap(tree, colors=NULL, fsize=1.0, ftype="reg", lwd=2, pts=FALSE, 
   node.numbers=FALSE, mar=NULL, add=FALSE, offset=NULL,
   direction="rightwards", type="phylogram", setEnv=TRUE, part=1.0, 
   xlim=NULL, ylim=NULL, nodes="intermediate", tips=NULL, maxY=NULL, 
   hold=TRUE, split.vertical=FALSE,
   lend=2, asp=NA, outline=FALSE, plot=TRUE)
\method{plot}{simmap}(x, ...)
\method{plot}{multiSimmap}(x, ...)
}
\arguments{
 	\item{tree}{a modified object of class "phylo" or "multiPhylo" containing a stochastic mapping or set of mappings (e.g., see \code{\link{read.simmap}} & \code{\link{make.simmap}}).}
 	\item{colors}{a vector with names translating the mapped states to colors - see \code{Examples}.}
 	\item{fsize}{relative font size for tip labels.}
	\item{ftype}{font type - options are \code{"reg"}, \code{"i"} (italics), \code{"b"} (bold), or \code{"bi"} (bold-italics).}
	\item{lwd}{line width for plotting.}
	\item{pts}{logical value indicating whether or not to plot filled circles at each vertex of the tree, as well as at transition points between mapped states.  Default is \code{FALSE}.}
	\item{node.numbers}{a logical value indicating whether or not node numbers should be plotted.}
	\item{mar}{vector containing the margins for the plot to be passed to \code{\link{par}}. If not specified, the default margins are [0.1,0.1,0.1,0.1].}
	\item{add}{a logical value indicating whether or not to add the plotted tree to the current plot (\code{TRUE}) or create a new plot (\code{FALSE}, the default).}
	\item{offset}{offset for the tip labels in character widths.}
	\item{direction}{plotting direction. Options are \code{"rightwards"} (the default), \code{"leftwards"}, \code{"upwards"} or \code{"downwards"}. For \code{method="fan"} \code{direction} is ignored.}
	\item{type}{plot type. Can be \code{"phylogram"}, \code{"fan"}, or \code{"cladogram"}. Only a subset of options are presently available for \code{type="fan"}.}
	\item{setEnv}{logical value indicating whether or not to set the environment \code{.PlotPhyloEnv}. Setting this to \code{TRUE} (the default) will allow compatibility with ape labeling functions such as \code{\link{nodelabels}}.}
	\item{part}{value between 0 and 1 for \code{type="fan"} indicating what fraction of the full circular tree to use as plotting area. For instance, \code{part=0.5} will plot a half fan phylogeny. It also affects the axis scaling used.}
	\item{xlim}{x-limits for the plot.}
	\item{ylim}{y-limits for the plot.}
	\item{nodes}{node placement following Felsenstein (2004; pp. 574-576). Can be \code{"intermediate"}, \code{"centered"}, \code{"weighted"}, or \code{"inner"}. So far only works for \code{type="phylogram"}.}
	\item{tips}{labeled vector containing the vertical position of tips. Normally this will be \code{1:N} for \code{N} tips in the tree.}
	\item{maxY}{maximum value of y to use before rotating a tree into fan configuration. This will only make a difference if different from \code{Ntip(tree)}.}
	\item{hold}{logical argument indicating whether or not to hold the output to the graphical device before plotting. Defaults to \code{hold=TRUE}.}
	\item{split.vertical}{split the color of the vertically plotted edges by the state of the daughter edges. Only applies if the edge state changes exactly at a node.}
	\item{lend}{line end style. See \code{\link{par}}.}
	\item{asp}{aspect ratio. See \code{\link{plot.window}}.}
	\item{outline}{logical value indicating whether or not to draw a black outline around the plotted edges of the tree.}
	\item{plot}{logical value indicating whether or not to actually plot the tree. (See equivalent argument in \code{\link{plot.phylo}}.)}
	\item{x}{for S3 plotting method, object of class \code{"simmap"} or \code{"multiSimmap"}.}
	\item{...}{for S3 plotting method, other arguments to be passed to \code{\link{plotSimmap}}.}
}
\description{
	Function plots one or multiple stochastic character mapped trees.
}
\details{
	The underscore character \code{"_"} is automatically swapped for a space in tip labels, as in \code{\link{plot.phylo}}.
}
\value{
	Plots a tree.
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Felsenstein, J. (2004) \emph{Inferring Phylogenies}. Sinauer.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback. (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{densityMap}}, \code{\link{make.simmap}}, \code{\link{read.simmap}} 
}
\examples{
# simulate a mapped tree
Q<-matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
rownames(Q)<-colnames(Q)<-letters[1:3]
tree<-sim.history(pbtree(n=100,scale=1),Q)
cols<-setNames(c("blue","red","green"),letters[1:3])
# plot the mapping
plot(tree,cols,ftype="i",fsize=0.7)
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{discrete character}
