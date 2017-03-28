\name{getDescendants}
\alias{getDescendants}
\title{Get descendant node numbers}
\usage{
getDescendants(tree, node, curr=NULL)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{node}{an integer specifying a node number in the tree.}
	\item{curr}{the set of previously stored node numbers - used in recursive function calls.}
}
\description{
	This function returns the set of node & tip numbers descended from \code{node}.
}
\value{
	The set of node and tip numbers for the nodes and tips descended from \code{node} in a vector.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{Descendants}}, \code{\link{paintSubTree}}
}
\keyword{phylogenetics}
\keyword{utilities}
