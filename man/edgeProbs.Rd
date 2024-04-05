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
	Computes the relative frequencies of character state changes along edges from a sample of stochastically mapped character histories. 
}
\details{	
	The function assumes that all trees in the sample differ only in their mapped histories & not at all in topology or branch lengths.
	
	Note that \code{edgeProbs} only asks whether the starting and ending states of the edge \emph{differ} in a particular way, and thus ignores multiple-hits along a single edge.
}
\value{
	The object that is returned is a matrix with the state changes & the relative frequency of each state change. Rows are in the order of the matrix \code{edge} for any of the mapped trees.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{plotSimmap}}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{discrete character}
