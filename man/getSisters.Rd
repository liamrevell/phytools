\name{getSisters}
\alias{getSisters}
\title{Get the sister node number, label, or set of nodes for a node or tip}
\usage{
getSisters(tree, node, mode=c("number","label"))
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{node}{a node number, tip number, node label, or tip label.}
	\item{mode}{an optional string indicating whether to return the node or tip number(s) or the node or tip label(s), if applicable.}
}
\description{
	Computes the sister taxon or node.
}
\details{
	This function takes a tree and node or tip number of label and returns the number or label of the sister or sisters to that node or tip.
}
\value{
	If \code{mode="number"} this function returns an integer or vector containing the node number of numbers of the sister node or tip. If \code{mode="label"} then this function returns a list containing up to two vectors: one for the node numbers of labels of sister nodes (if applicable); and the other containing the tip labels of the sister tips.
}
\seealso{
	\code{\link{getDescendants}}, \code{\link[phangorn]{Siblings}}
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
