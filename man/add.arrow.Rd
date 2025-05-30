\name{add.arrow}
\alias{add.arrow}
\title{Add an arrow pointing to a tip or node on the tree}
\usage{
add.arrow(tree=NULL, tip, ...)
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}, \code{"contMap"}, or \code{"densityMap"}. If not supplied, the function will obtain the last plotted phylogeny from the environmental variable \code{last_plot.phylo}.}
	\item{tip}{label of tip or tip or node number; or vector of such values. If \code{tree=NULL} then the tip label(s) or node number(s) must be supplied.}
	\item{...}{optional arguments to control the shape and size of the arrow including: its length (\code{arrl}) in the units of the plot; the length of the arrowhead (\code{hedl}); the total angle between the wings in the arrowhead (\code{angle}); the line width for the plotted lines (\code{lwd}); the offset from the tip or end of tip label, in character widths (\code{offset}); and the color (\code{col}).}
}
\description{
	Adds an arrow or a set of arrows to a plotted tree.
}
\details{
	Trees can be plotted using \pkg{phytools} function \code{plotTree}, \code{plotSimmap}, \code{contMap}, \code{densityMap}, and \pkg{ape} S3 method \code{plot.phylo}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{nodelabels}}
}
\examples{
## load tree
data(anoletree)

## plot tree
plotTree(anoletree,type="fan",fsize=0.7,ftype="i")

## add lwd=5 arrows in black
add.arrow(anoletree,tip=c("cuvieri","krugi",
    "pulchellus","poncensis","stratulus",
    "evermanni","cooki","cristatellus",
    "gundlachi","occultus"),lwd=5,arrl=1)

## overlay thinner arrows to create outline effect
add.arrow(anoletree,tip="cuvieri",col="green",
    lwd=3,arrl=1)
add.arrow(anoletree,tip=c("krugi","pulchellus",
    "poncensis"),col="#E4D96F",lwd=3,arrl=1)
add.arrow(anoletree,tip=c("stratulus","evermanni"),
    col="darkgreen",lwd=3,arrl=1)
add.arrow(anoletree,tip=c("cooki","cristatellus",
    "gundlachi"),col="brown",lwd=3,arrl=1)
add.arrow(anoletree,tip="occultus",col="darkgrey",
    lwd=3,arrl=1)
	
## add legend
legend(x="topleft",c("crown-giant","grass-bush","trunk-crown","trunk-ground",
    "twig"),pch=22,pt.bg=c("green","#E4D96F","darkgreen",
    "brown","darkgrey"),cex=0.9,
    pt.cex=2,title="PR ecomorphs",bty="n")
	
## reset plotting parameters
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{plotting}
