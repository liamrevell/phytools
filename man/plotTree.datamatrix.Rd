\name{plotTree.datamatrix}
\alias{plotTree.datamatrix}
\title{Plot a tree with a discrete character data matrix at the tips}
\usage{
plotTree.datamatrix(tree, X, ...)
}
\arguments{
 	\item{tree}{an object of class \code{"phylo"}.}
	\item{X}{a data frame with columns as factors.}
	\item{...}{optional arguments.}
}
\description{
	Plots a phylogeny next to a matrix of discrete character data.
}
\value{
	Invisibly returns a list containing the font size, a list of the colors used for each column of the plotted data matrix, and the x-coordinate of the rightmost edge of the matrix.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{dotTree}}, \code{\link{phylo.heatmap}}
}
\keyword{phylogenetics}
\keyword{plotting}
\keyword{comparative method}
\keyword{discrete character}
