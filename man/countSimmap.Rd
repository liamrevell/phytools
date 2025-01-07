\name{countSimmap}
\alias{countSimmap}
\title{Counts the number of character changes on a object of class \code{"simmap"} or \code{"multiSimmap"}}
\usage{
countSimmap(tree, states=NULL, message=TRUE)
}
\arguments{
	\item{tree}{an object of class \code{"simmap"} or \code{"multiSimmap"}.}
	\item{states}{optional argument with the states for the mapped character. If not provided, these will be computed from the tree. This is useful if averaging across many trees, some of which may lack certain states.}
	\item{message}{optional logical argument indicating whether or not to return an informative message about the function output.}
}
\description{
	Count the number of changes of different types on a stochastically mapped trees or set of trees (objects of class \code{"simmap"} or \code{"multiSimmap"}).
}
\details{
	This function takes a tree or a set of trees with a mapped discrete character (that is, an object of class \code{"simmap"} or \code{"multiSimmap"}), and computes the total number of character changes as well as the number of character changes between all states.
}
\value{
	A list with up to three elements: \code{N} is an integer value giving the total number of character changes on the tree; \code{Tr} gives the number of of transitions between row and column states (or a matrix containing both \code{N} and the transitions between states, in rows, for an object of class \code{"multiPhylo"}); and (optionally) \code{message} contains an explanatory message about the function  output.
}
\examples{
## load anole tree with ecomorphs mapped
data(anoletree)
anoletree

## countSimmap
countSimmap(anoletree)
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{utilities}
\keyword{discrete character}
