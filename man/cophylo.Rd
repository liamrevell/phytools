\name{cophylo}
\alias{cophylo}
\alias{plot.cophylo}
\alias{cotangleplot}
\title{Creates a co-phylogenetic plot}
\usage{
cophylo(tr1, tr2, assoc=NULL, rotate=TRUE, ...)
\method{plot}{cophylo}(x, ...)
cotangleplot(tr1, tr2, type=c("cladogram","phylogram"),
	use.edge.length=TRUE, tangle=c("both","tree1","tree2"), ...)
}
\arguments{
	\item{tr1}{object of class \code{"phylo"}.}
	\item{tr2}{object of class \code{"phylo"}.}
	\item{assoc}{matrix containing the tip labels in \code{tr1} to match to the tip labels in \code{tr2}. Note that not all labels in either tree need to be included; and, furthermore, one label in \code{tr1} can be matched with more than one label in \code{tr2}, or vice versa.}
	\item{rotate}{logical argument indicating whether nodes on both trees should be rotated to attempt to match in vertical position.}
	\item{x}{in the case of \code{plot.cophylo}, an object of class \code{"cophylo"} to be plotted.}
	\item{type}{for \code{cotangleplot}, the tree plotting style.}
	\item{use.edge.length}{for \code{cotangleplot}, a logical value indicating whether or not to plot trees with edge lengths.}
	\item{tangle}{for \code{cotangleplot}, whether to tangle the left tree, the right tree, or both.}
	\item{...}{optional arguments to be passed to \code{\link{tipRotate}}, or, in the case of \code{plot.cophylo}, to the internally used tree plotting function, \code{phylogram}. \code{phylogram} takes similar arguments to \code{\link{plotSimmap}}, such as \code{fsize}, \code{ftype}, \code{lwd}, and \code{pts}, though not all options from \code{plotSimmap} and \code{plotTree} are available. If \code{fsize} is supplied as a vector, different size fonts for the left & right facing trees may be used. In addition, the optional argument \code{scale.bar}, which should be a vector containing the lengths of the scale bars desired for the right & left trees, will add scale bars to the plot when supplied to \code{plot.cophylo}. If either tree contains polytomies, the \code{cophylo} argument \code{rotate.multi} should be set to \code{TRUE}. If curved linking lines are desired, the \code{plot.cophylo} argument \code{link.type} should be set to \code{"curved"}. Other arguments for the \code{plot} method include \code{link.col}, \code{link.lty}, and \code{link.lwd}, which can be supplied as a scalar or a vector in which the order of the elements corresponds to the order of the associations in \code{assoc}. Finally, \code{edge.col}, a list consisting of two vectors (\code{left} and \code{right}) can be used to specify the edge colors of the two left & right plotted trees. Note that the edge order is the same as in the \emph{rotated} trees, assuming that a rotation has been performed on \code{x}.}
}
\description{
	This function creates an object of class \code{"cophylo"} or, in the case of \code{plot.cophylo}, plots that object. The function can (optionally) first attempt to rotate the nodes of both trees to optimize vertical matching of tips.
	
	\code{cotangleplot} creates a co-phylogenetic plot in which the edges of the matched trees are crossing.
}
\details{
	If no matrix of associations, \code{assoc}, is provided, then \code{cophylo} will look for exact matches of tip labels between trees.
}
\value{
	An object of class \code{"cophylo"} which includes the following components or a pair of plotted facing phylogenies with links between tips as specified in \code{assoc}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{cophyloplot}}, \code{\link{plotSimmap}}
}
\examples{
## load data from Lopez-Vaamonde et al. (2001)
data(wasp.trees)
data(wasp.data)
## create co-phylogenetic object
wasp.cophylo<-cophylo(wasp.trees[[1]],wasp.trees[[2]],
    assoc=wasp.data)
## plot co-phylogenies
plot(wasp.cophylo,link.type="curved",link.lwd=4,
     link.lty="solid",link.col=make.transparent("red",
     0.25))
par(mar=c(5.1,4.1,4.1,2.1))
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{diversification}
\keyword{co-phylogenetics}
