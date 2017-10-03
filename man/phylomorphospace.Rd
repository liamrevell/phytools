\name{phylomorphospace}
\alias{phylomorphospace}
\title{Creates phylomorphospace plot}
\usage{
phylomorphospace(tree, X, A=NULL, label=c("radial","horizontal","off"),
   control=list(), ...)
}
\arguments{
	\item{tree}{a phylogenetic tree in \code{"phylo"} format, or a modified \code{"phylo"} object with a mapped discrete character.}
	\item{X}{an \code{n} x 2 matrix of tip values for two characters in \code{n} species.}
	\item{A}{an optional \code{m} x 2 matrix (for \code{m} nodes) of values for two taits at internal nodes in the tree - if not supplied, these values will be estimated using \code{\link{fastAnc}}.}
	\item{label}{string indicating whether to plot the tip labels in the same direction as the terminal edge (\code{label="radial"}), horizontally \code{label="horizontal"}, or not at all \code{"off"}. \code{label=TRUE} and \code{label=FALSE} are also acceptable, for compatibility with phytools <= 0.3-03.}
	\item{control}{a list containing the following optional control parameters: \code{col.edge}: a vector of edge colors; and \code{col.node}: a vector of node colors.}
	\item{...}{optional arguments for plotting, including \code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, \code{lwd}, \code{colors}, \code{fsize}, and \code{node.by.map}. \code{colors} is only used when there is a mapped discrete character on the tree, in which case \code{control$col.edge} is ignored. \code{fsize} is relative to the default, which is \code{textxy(...,cx=0.75)}. \code{node.by.map} is a logical value (defaults to \code{FALSE} which tells the function whether or not to plot the node colors using the colors of the mapped discrete character. Setting this option to \code{TRUE} will cause \code{control$col.node} to be ignored. \code{node.size} is a vector containing the point size relative to the default (see \code{\link{par}} for plotted internal nodes and tips, respectively. Defaults to \code{node.size=c(1,1.3)}. If only one number is provided it will be recycled. \code{axes} is a logical value indicating whether or not axes should be plotted (see \code{\link{plot.default}}. Finally, \code{add} indicates whether to add the phylmorphospace to the current plot.}
}
\description{
	This function creates a phylomorphospace plot (a projection of the tree into morphospace) for two characters following Sidlauskas (2008; \emph{Evolution}). It will also plot a discrete character mapped on tree.
}
\value{
	This function creates a phylomorphospace plot
}
\references{
	Paradis, E., J. Claude, and K. Strimmer (2004) APE: Analyses of phylogenetics and evolution in R language. \emph{Bioinformatics}, \bold{20}, 289-90.

	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.

	Sidlauskas, B. (2008) Continuous and arrested morphological diversification in sister clades of characiform fishes: A phylomorphospace approach. \emph{Evolution}, \bold{62}, 3135-3156.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
tree<-pbtree(n=25)
X<-fastBM(tree,nsim=2)
phylomorphospace(tree,X,xlab="trait 1",ylab="trait 2")
}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{plotting}
