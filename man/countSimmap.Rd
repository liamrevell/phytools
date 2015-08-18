\name{countSimmap}
\alias{countSimmap}
\title{Counts the number of character changes on a SIMMAP style tree or set of trees}
\usage{
countSimmap(tree, states=NULL, message=TRUE)
}
\arguments{
	\item{tree}{a single tree or a set of trees with a mapped discrete character (e.g, see \code{\link{make.simmap}} or \code{\link{read.simmap}}.}
	\item{states}{optional argument with the states for the mapped character. If not provided, these will be computed from the tree. This is useful if averaging across many trees, some of which may lack certain states.}
	\item{message}{optional logical argument indicating whether or not to return an informative message about the function output.}
}
\description{
	This function takes a tree or a set of trees with a mapped discrete character (SIMMAP style, e.g., see \code{\link{make.simmap}} or \code{\link{read.simmap}}), and computes the total number of character changes as well as the number of character changes between all states.
}
\value{
	A list with up to three elements: \code{N} is an integer value giving the total number of character changes on the tree; \code{Tr} gives the number of of transitions between row and column states (or a matrix containing both \code{N} and the transitions between states, in rows, for an object of class \code{"multiPhylo"}); and (optionally) \code{message} contains an explanatory message about the function  output.
}
\examples{
tree<-pbtree(n=100,scale=1)
Q<-matrix(c(-2,1,1,1,-2,1,1,1,-2),3,3)
colnames(Q)<-rownames(Q)<-c("A","B","C")
mtree<-sim.history(tree,Q)
countSimmap(mtree,states=rownames(Q))
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
