\name{phylomorphospace}
\alias{phylomorphospace}
\alias{project.phylomorphospace}
\title{Creates phylomorphospace plot}
\usage{
phylomorphospace(tree, X, A=NULL, label=c("radial","horizontal","off"),
    control=list(), ...)
project.phylomorphospace(tree, X, nsteps=200, sleep=0,
    direction=c("to","from","both"), ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format, or a modified \code{"phylo"} object with a mapped discrete character.}
	\item{X}{an \code{n} x 2 matrix of tip values for two characters in \code{n} species.}
	\item{A}{an optional \code{m} x 2 matrix (for \code{m} nodes) of values for two traits at internal nodes in the tree - if not supplied, these values will be estimated using \code{\link{fastAnc}}.}
	\item{label}{string indicating whether to plot the tip labels in the same direction as the terminal edge (\code{label="radial"}), horizontally \code{label="horizontal"}, or not at all \code{"off"}. \code{label=TRUE} and \code{label=FALSE} are also acceptable, for compatibility with phytools <= 0.3-03.}
	\item{control}{a list containing the following optional control parameters: \code{col.edge}: a vector of edge colors; and \code{col.node}: a vector of node colors.}
	\item{nsteps}{for \code{project.phylomorphospace} the number of frames in the animation between the phylogeny & the phylomorphospace or \emph{vice versa}.}
	\item{sleep}{for \code{project.phylomorphospace} the time between frames.}
	\item{direction}{for \code{project.phylomorphospace} whether to morph \code{"to"} a phylomorphospace, \code{"from"} a phylomorphospace, or there & back again (\code{"both"}).}
	\item{...}{optional arguments for plotting, including \code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, \code{lwd}, \code{colors}, \code{fsize}, and \code{node.by.map}. \code{colors} is only used when there is a mapped discrete character on the tree, in which case \code{control$col.edge} is ignored. \code{fsize} is relative to the default, which is \code{textxy(...,cx=0.75)}. \code{node.by.map} is a logical value (defaults to \code{FALSE} which tells the function whether or not to plot the node colors using the colors of the mapped discrete character. Setting this option to \code{TRUE} will cause \code{control$col.node} to be ignored. \code{node.size} is a vector containing the point size relative to the default (see \code{\link{par}} for plotted internal nodes and tips, respectively. Defaults to \code{node.size=c(1,1.3)}. If only one number is provided it will be recycled. \code{axes} is a logical value indicating whether or not axes should be plotted (see \code{\link{plot.default}}. Finally, \code{add} indicates whether to add the phylomorphospace to the current plot.}
}
\description{
	Project a phylogeny into morphospace.
}
\details{
	This function creates a phylomorphospace plot (a projection of the tree into morphospace) for two characters following Sidlauskas (2008; \emph{Evolution}). It will also plot a discrete character mapped on tree.
	
	\code{project.phylomorphospace} animates the phylomorphospace projection.
}
\value{
	This function creates a phylomorphospace plot
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-90.

	Revell, L. J. (2014) Graphical methods for visualizing comparative data on phylogenies. Chapter 4 in \emph{Modern phylogenetic comparative methods and their application in evolutionary biology: Concepts and practice} (L. Z. Garamszegi ed.), pp. 77-103.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.

	Sidlauskas, B. (2008) Continuous and arrested morphological diversification in sister clades of characiform fishes: A phylomorphospace approach. \emph{Evolution}, \bold{62}, 3135-3156.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## load tree & data from Revell & Collar (2009)
data(sunfish.tree)
data(sunfish.data)

## set colors for mapped discrete character
cols<-setNames(c("blue","red"),
    levels(sunfish.data$feeding.mode))

## graph phylomorphospace
phylomorphospace(sunfish.tree,sunfish.data[,3:2],
    colors=cols,bty="l",ftype="off",node.by.map=TRUE,
    node.size=c(0,1.2),xlab="relative buccal length",
    ylab="relative gape width")
title(main="Phylomorphospace of buccal morphology in Centrarchidae",
    font.main=3)
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{plotting}
\keyword{continuous character}
