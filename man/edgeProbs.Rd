\name{edgeProbs}
\alias{edgeProbs}
\title{Compute the relative frequencies of state changes along edges}
\usage{
edgeProbs(trees)
}
\arguments{
	\item{trees}{an object of class \code{"multiSimmap"} containing a sample of trees that are identical in topology & branch lengths with different stochastically mapped character histories.}
}
\description{
	This function computes the relative frequencies of character state changes along edges from a sample of stochastically mapped character histories. This function assumes that all trees in the sample differ only in their mapped histories & not at all in topology or branch lengths. Note that it only asks whether the starting and ending states of the edge differ in a particular way, and thus ignores multiple-hits along a single edge.
}
\value{
	The object that is returned is a matrix with the state changes & the relative frequency of each state change. Rows are in the order of the matrix \code{edge} for any of the mapped trees.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}
}
\keyword{phylogenetics}
\keyword{plotting}
