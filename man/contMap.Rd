\name{contMap}
\alias{contMap}
\alias{plot.contMap}
\alias{errorbar.contMap}
\title{Map continuous trait evolution on the tree}
\usage{
contMap(tree, x, res=100, fsize=NULL, ftype=NULL, lwd=4, legend=NULL,
   lims=NULL, outline=TRUE, sig=3, type="phylogram", direction="rightwards", 
   plot=TRUE, ...)
\method{plot}{contMap}(x, ...)
errorbar.contMap(obj, ...)
}
\arguments{
 	\item{tree}{object of class \code{"phylo"}.}
	\item{x}{vector of phenotypic trait values for species. \code{names(x)} should contain the species names and match \code{tree$tip.label}. Or, for \code{plot.contMap}, an object of class \code{"contMap"}.}
	\item{res}{resolution for gradient plotting. Larger numbers indicate a finer (smoother) gradient.}
	\item{fsize}{relative font size - can be a vector with the second element giving the font size for the legend.}
	\item{ftype}{font type - see options in \code{\link{plotSimmap}}. As with \code{fsize}, this can be a vector with the second element giving font type for the legend.}
	\item{lwd}{line width for branches. Can be a single integer number or a vector. In the latter case, the second number will be taken to be the desired legend width.}
	\item{legend}{if \code{FALSE} no legend is plotted; if a numeric value, it gives the length of the legend in units of branch length. Default is 0.5 times the total tree length.}
	\item{lims}{range for the color map. By default, this will be \code{c(min(x),max(x))}, and should always include this range.}
	\item{outline}{logical value indicating whether or not to outline the branches of the tree in black.}
	\item{sig}{the number of decimal places to show on the legend limits.}
	\item{type}{type of plot desired. Options are \code{"phylogram"} for a rightward square phylogram; and \code{"fan"} for a circular phylogram.}
	\item{direction}{plotting direction for \code{type="phylogram"}.}
	\item{plot}{logical value indicating whether or not to plot the tree. If \code{plot=FALSE} then an object of class \code{"contMap"} will be returned without plotting.}
	\item{obj}{object of class \code{"contMap"}.}
	\item{...}{optional arguments for \code{plot.contMap} which include all the arguments of \code{contMap} except for \code{tree}, \code{x}, \code{res}, and \code{lims}. Also \code{method}, \code{"fastAnc"}, \code{"anc.ML"}, or \code{"user"} (for user-supplied states) specifying which function to use for ancestral state estimation; \code{hold} specifies whether or not to hold output to graphical device before plotting (defaults to \code{hold=TRUE}); and \code{anc.states} a vector containing some or multiple ancestral user-supplied ancestral states at nodes. Some other plotting arguments, such as \code{xlim} and \code{ylim}, may also work. Optional arguments for \code{errorbar.contMap} include \code{x}, a vector containing the original trait values mapped onto the tree (otherwise these will be obtained from \code{obj}), \code{scale.by.ci}, a logical argument (defaulting to \code{TRUE}) that determines whether or not the length of the error bars will be scaled by the CI width, and \code{lwd}, which detemines the line width of the plotted error bars.}
}
\description{
	Function plots a tree with a mapped continuous character. The mapping is accomplished by estimating states at internal nodes using ML with \code{\link{fastAnc}}, and the interpolating the states along each edge using equation [2] of Felsenstein (1985).
	
	\code{errorbar.contMap} adds error bars to an existing plot.
}
\value{
	Plots a tree. An object of class \code{"contMap"} is returned invisibly.
	
	\code{errorbar.contMap} adds colorful error bars to a plotted tree.
}
\references{
	Felsenstein, J. 1985. Phylogenies and the comparative method. \emph{American Naturalist}, \bold{125}, 1-15. 

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Revell, L. J. 2013. Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.ML}}, \code{\link{densityMap}}, \code{\link{fastAnc}}, \code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
