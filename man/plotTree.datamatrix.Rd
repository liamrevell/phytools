\name{plotTree.datamatrix}
\alias{plotTree.datamatrix}
\alias{plotFanTree.wTraits}
\title{Plot a tree with a discrete (or continuous) character data matrix at the tips}
\usage{
plotTree.datamatrix(tree, X, ...)
plotFanTree.wTraits(tree, X, type=c("arc","fan"), ...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a data frame with columns as factors for \code{plotTree.datamatrix}, or a matrix of data frame containing continuous and/or discrete characters for \code{plotFanTree.wTraits}.}
	\item{type}{plot type for \code{plotFanTree.wTraits}.}
	\item{...}{optional arguments.}
}
\description{
	\code{plotTree.datamatrix} plots a phylogeny next to a matrix of discrete characters. \code{plotFanTree.wTraits} plots an arc or fan style tree with discrete or continuous data at the tips.
}
\value{
	\code{plotTree.datamatrix} invisibly returns a list containing the font size, a list of the colors used for each column of the plotted data matrix, and the x-coordinate of the rightmost edge of the matrix.
	
	\code{plotFanTree.wTraits} invisibly returns the color palette used for plotting.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{dotTree}}, \code{\link{phylo.heatmap}}
}
\examples{
data("liolaemid.data")
data("liolaemid.tree")
liolaemid.data[liolaemid.tree$tip.label,]->liolaemid.data
colors<-list(
  c("blue","white","red"),
  terrain.colors(n=10),
  setNames(c("#F0EAD6","#DF536B"),c("O","V")))
cols<-plotFanTree.wTraits(liolaemid.tree,
  liolaemid.data[,3:1],lwd=12,colors=colors,ftype="off",
  spacer=0)
legend(x=0,y=0.7*max(nodeHeights(liolaemid.tree)),
  names(colors[[3]]),lwd=8,col=colors[[3]],
  title="parity mode",bty="n",xjust=0.5,yjust=0.5)
add.color.bar(1.5*max(nodeHeights(liolaemid.tree)),cols[[2]],
  title="maximum altitude (m)",
  lims=range(liolaemid.data[,2]),digits=2,prompt=FALSE,
  x=-0.75*max(nodeHeights(liolaemid.tree)),
  y=0.2*max(nodeHeights(liolaemid.tree)),subtitle="",
  lwd=8,outline=FALSE)
add.color.bar(1.5*max(nodeHeights(liolaemid.tree)),cols[[1]],
  title="environmental temp.",
  lims=range(liolaemid.data[,3]),digits=2,prompt=FALSE,
  x=-0.75*max(nodeHeights(liolaemid.tree)),
  y=-0.15*max(nodeHeights(liolaemid.tree)),subtitle="",
  lwd=8,outline=FALSE)
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margins to default
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{discrete character}
