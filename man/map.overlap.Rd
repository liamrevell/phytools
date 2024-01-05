\name{map.overlap}
\alias{map.overlap}
\alias{Map.Overlap}
\title{Proportional overlap between two mapped character histories on a tree}
\usage{
map.overlap(tree1, tree2, tol=1e-6, ...)
Map.Overlap(tree1, tree2, tol=1e-06, standardize=TRUE, ...)
}
\arguments{
	\item{tree1}{an object of class \code{"simmap"}.}
	\item{tree2}{an object of class \code{"simmap"}.}
	\item{tol}{an optional tolerance value.}
	\item{standardize}{for \code{Map.Overlap}, a logical value indicating whether or not to standardize overlap by dividing by the summed branch length of the tree.}
	\item{...}{optional arguments, such as \code{check.equal}, a logical value indicating whether or not to check if \code{tree1} and \code{tree2} match in underlying topology and branch lengths (they should). This value is \code{TRUE} by default, but can be set to \code{FALSE} if \code{tree1} and \code{tree2} are known to be equal to speed up calculation.}
}
\description{
	Calculates the similarity of two different stochastically mapped character histories.
}
\details{	
	\code{map.overlap} computes a single quantity giving the overall similarity of the maps, consequently this measure only makes sense of some or all of the states are shared between the two mapped tress. 
	
	In \code{Map.Overlap} what is computed instead is a matrix in which the rows correspond to the states observed in \code{tree1} and the columns give the states for \code{tree2}, with the numerical values of the matrix showing the total overlap between each pair of mapped states on the two trees.
}
\value{
	A numerical value on the interval (0, 1), for \code{map.overlap}; or a matrix whose elements should sum to 1.0 (\code{Map.Overlap}).
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link{read.simmap}}
}
\keyword{ancestral states}
\keyword{phylogenetics}
\keyword{comparative method}
\keyword{discrete character}
