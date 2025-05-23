\name{add.random}
\alias{add.random}
\title{Add tips at random to the tree}
\usage{
add.random(tree, n=NULL, tips=NULL, edge.length=NULL, order=c("random","input"))
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{n}{a number of tips to add to the tree. If \code{NULL}, will use \code{length(tips)}. If \code{tips} is specified, then \code{n} is also ignored, regardless of its value.}
	\item{tips}{a set of tip names for the added tips. If \code{NULL}, names will be supplied}
	\item{edge.length}{terminal edge length for the added tips. If \code{NULL}, and \code{is.ultrametric(tree)} \code{==TRUE}, then edge lengths will be assigned to keep the tree ultrametric. Note that if edge lengths are assigned and \code{n>1}, then the assigned terminal edge lengths are not guaranteed as subsequent random tip addition could occur along the new terminal edge.}
	\item{order}{addition order for the new tips.}
}
\description{
	This function adds new tips at random to a tree with branch lengths.
}
\details{
	If no edge lengths are provided, and the tree is ultrametric, then edge lengths are assigned to keep the tree ultrametric.
	
	The probability that at new tip is added along any branch is directly proportional to the length of the branch.

	Note that sometimes the resultant tree plotted with \code{\link[ape]{plot.phylo}} or \code{\link{plotSimmap}} may display with branches crossing. If so, the tree can be 'untangled' using \code{\link{untangle}}.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{allFurcTrees}}, \code{\link{add.everywhere}}
}
\keyword{phylogenetics}
\keyword{utilities}
