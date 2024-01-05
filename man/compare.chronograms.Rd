\name{compare.chronograms}
\alias{compare.chronograms}
\title{Compares two chronograms with precisely matching nodes in a visual manner}
\usage{
compare.chronograms(t1, t2, ...)
}
\arguments{
	\item{t1}{object of class \code{"phylo"}.}
	\item{t2}{object of class \code{"phylo"} that matches \code{t1} precisely in topology & node rotations, but differs in edge lengths.}
	\item{...}{optional arguments.}
}
\description{
	Creates a visual comparison of two chronograms that differ in edge lengths but not topology.
}
\details{
	This function plots two trees, with semi-transparent colors by default, & uses arrows to highlight differences in depth of corresponding nodes between the trees.
}
\value{
	Function creates a plot.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
