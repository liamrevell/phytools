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
	\item{x}{a numerical vector of phenotypic trait values for species. \code{names(x)} should contain the species names and match \code{tree$tip.label}. Or, for \code{plot.contMap}, an object of class \code{"contMap"}.}
	\item{res}{resolution for gradient plotting. Larger numbers (to a point) indicate a finer (smoother) gradient.}
	\item{fsize}{relative font size - can be a vector of length 2 in which the first element gives the font size for the tip labels & the second element giving the font size for the legend.}
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
	\item{...}{optional arguments for \code{plot.contMap} which include all the arguments of \code{contMap} except for \code{tree}, \code{x}, \code{res}, and \code{lims}. Also \code{method}, \code{"fastAnc"}, \code{"anc.ML"}, or \code{"user"} (for user-supplied states) specifying which function to use for ancestral state estimation; \code{hold} specifies whether or not to hold output to graphical device before plotting (defaults to \code{hold=TRUE}); and \code{anc.states} a vector containing some or multiple ancestral user-supplied ancestral states at nodes. Some other plotting arguments, such as \code{xlim} and \code{ylim}, may also work. Optional arguments for \code{errorbar.contMap} include \code{x}, a vector containing the original trait values mapped onto the tree (otherwise these will be obtained from \code{obj}), \code{scale.by.ci}, a logical argument (defaulting to \code{TRUE}) that determines whether or not the length of the error bars will be scaled by the CI width, and \code{lwd}, which determines the line width of the plotted error bars.}
}
\description{
	Projects the observed and reconstructed values of a continuous trait onto the edges of a tree using a color gradient.
}
\details{
	Function plots a tree with a mapped continuous character.
	
	The mapping is accomplished by estimating states at internal nodes using ML with \code{\link{fastAnc}}, and then interpolating the states along each edge using equation [2] of Felsenstein (1985). This makes these interpolated values equal to the maximum likelihood estimates under a Brownian evolutionary process.
	
	The default color palette is \emph{not} recommended as it is not colorblind friendly and does not render well into gray scale; however, this can be updated using the helper function \code{\link{setMap}}.
	
	\code{errorbar.contMap} adds error bars to an existing plot.
}
\value{
	Plots a tree. An object of class \code{"contMap"} is returned invisibly.
	
	\code{errorbar.contMap} adds colorful error bars to a plotted tree.
}
\references{
	Felsenstein, J. (1985) Phylogenies and the comparative method. \emph{American Naturalist}, \bold{125}, 1-15. 

	Revell, L. J. (2013) Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
	
	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
	
	Revell, L. J., K. Schliep, E. Valderrama, and J. E. Richardson (2018) Graphs in phylogenetic comparative analysis: Anscombe's quartet revisited. \emph{Methods in Ecology and Evolution}, \bold{9}, 2145-2154.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{anc.ML}}, \code{\link{densityMap}}, \code{\link{fastAnc}}, \code{\link{plotSimmap}}
}
\examples{
## load data from Garland et al. (1992)
data(mammal.tree)
data(mammal.data)
## extract character of interest
ln.bodyMass<-log(setNames(mammal.data$bodyMass,
    rownames(mammal.data)))
## create "contMap" object
mammal.contMap<-contMap(mammal.tree,
    ln.bodyMass,plot=FALSE,res=200)
## change color scheme
mammal.contMap<-setMap(mammal.contMap,
    c("white","#FFFFB2","#FECC5C","#FD8D3C",
    "#E31A1C"))
plot(mammal.contMap,fsize=c(0.7,0.8),
    leg.txt="log(body mass)")
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{continuous character}
\keyword{maximum likelihood}
