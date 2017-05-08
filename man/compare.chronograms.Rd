\name{compare.chronograms}
\alias{compare.chronograms}
\title{Compares to chronograms with precisely matching nodes in a visual manner}
\usage{
compare.chronograms(t1, t2, ...)
make.transparent(color, alpha)
}
\arguments{
	\item{t1}{object of class \code{"phylo"}.}
	\item{t2}{object of class \code{"phylo"} that matches \code{t1} precisely in topology & node rotations, but differs in edge lengths.}
	\item{...}{optional arguments.}
}
\description{
	This function plots two trees, with semi-transparent colors by default, & uses arrows to highlight differences in depth of corresponding nodes between the trees.
}
\value{
	Function creates a plot.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
