\name{densityMap}
\alias{densityMap}
\alias{plot.densityMap}
\title{Plot posterior density of stochastic mapping on a tree}
\usage{
densityMap(trees, res=100, fsize=NULL, ftype=NULL, lwd=3, check=FALSE,
   legend=NULL, outline=FALSE, type="phylogram", direction="rightwards",
   plot=TRUE, ...)
\method{plot}{densityMap}(x, ...)
}
\arguments{
 	\item{trees}{set of phylogenetic trees in a modified \code{"multiPhylo"} object. Values for a two-state discrete character are mapped on the tree. See \code{\link{make.simmap}} and \code{\link{read.simmap}} for details.}
	\item{res}{resolution for gradient plotting. Larger numbers indicate a finer (smoother) gradient.}
	\item{fsize}{relative font size - can be a vector with the second element giving the font size for the legend.}
	\item{ftype}{font type - see options in \code{\link{plotSimmap}}. As with \code{fsize}, can be a vector with the second element giving font type for the legend.}
	\item{lwd}{line width for branches. If a vector of two elements is supplied, the second element will be taken to be the desired width of the legend bar.}
	\item{check}{check to make sure that the topology and branch lengths of all phylogenies in \code{trees} are equal.}
	\item{legend}{if \code{FALSE} no legend is plotted; if a numeric value, it gives the length of the legend in units of branch length. Default is 0.5 times the total tree length.}
	\item{outline}{logical value indicating whether or not to outline the branches of the tree in black.}
	\item{type}{type of plot desired. Options are \code{"phylogram"} for a rightward square phylogram; and \code{"fan"} for a circular phylogram.}
	\item{plot}{logical value indicating whether or not to plot the tree. If \code{plot=FALSE} then an object of class \code{"densityMap"} will be returned without plotting.}
	\item{direction}{plotting direction for \code{type="phylogram"}.}
	\item{x}{for \code{plot.densityMap}, an object of class \code{"densityMap"}.}
	\item{...}{optional arguments for \code{plot.densityMap}. These include all the arguments of \code{densityMap} except \code{trees} and \code{res}. Additional optional arguments include \code{mar} (margins), \code{offset} (tip label offset: in units of the edge length or character widths, as in \code{\link{plotSimmap}}), and \code{hold} (whether or not to use \code{dev.hold} to hold output to graphical device before plotting; defaults to \code{hold=TRUE}). Also, the argument \code{states} can be used to 'order' the states on the probability axis (that is, which state should correspond to a posterior probability of 0 or 1). Some other plotting arguments, such as \code{xlim} and \code{ylim}, may also work.}
}
\description{
	Visualize posterior probability density from stochastic mapping using a color gradient on the tree.
}
\details{
	Function plots a tree with the posterior density for a mapped character from stochastic character mapping on the tree. Since the mapped value is the probability of being in state "1", only binary [0,1] characters are allowed.
}
\value{
	Plots a tree and returns an object of class \code{"densityMap"} invisibly.
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Huelsenbeck, J. P., R. Neilsen, and J. P. Bollback (2003) Stochastic mapping of morphological characters. \emph{Systematic Biology}, \bold{52}, 131-138. 

	Revell, L. J. (2013) Two new graphical methods for mapping trait evolution on phylogenies. \emph{Methods in Ecology and Evolution}, \bold{4}, 754-759.
	
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{plotSimmap}}, \code{\link{read.simmap}}
}
\examples{
\dontrun{
## load tree and data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)
## extract discrete character (feeding mode)
fmode<-setNames(sunfish.data$feeding.mode,
     rownames(sunfish.data))
## do stochastic mapping
smap.trees<-make.simmap(sunfish.tree,fmode,model="ER",
     nsim=100)
## compute "densityMap" object
sunfish.dmap<-densityMap(smap.trees,plot=FALSE,
	res=50) ## res should be higher
## plot density map
plot(sunfish.dmap,lwd=5,outline=TRUE)
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{discrete character}
