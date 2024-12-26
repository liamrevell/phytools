\name{getDescendants}
\alias{getDescendants}
\alias{getParent}
\title{Get descendant node numbers}
\usage{
getDescendants(tree, node, curr=NULL)
getParent(tree, node)
}
\arguments{
	\item{tree}{a phylogenetic tree as an object of class \code{"phylo"}.}
	\item{node}{an integer specifying a node number in the tree.}
	\item{curr}{the set of previously stored node numbers - used in recursive function calls.}
}
\description{
	Returns the descendants or parent of a specified node.
}
\details{
	\code{getDescendants} returns the set of node & tip numbers descended from \code{node}.
	
	\code{getParent} returns the \emph{single} parent node of a specified node number (or \code{NULL} if \code{node} is already the root).
}
\value{
	The set of node and tip numbers for the nodes and tips descended from \code{node} in a vector, or for \code{getParent} the single node preceding \code{node} in the tree.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[phangorn]{Descendants}}, \code{\link{paintSubTree}}
}
\keyword{phylogenetics}
\keyword{utilities}
