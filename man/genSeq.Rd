\name{genSeq}
\alias{genSeq}
\title{Simulate a DNA alignment on the tree under a model}
\usage{
genSeq(tree, l=1000, Q=NULL, rate=1, format="DNAbin", ...)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{l}{length of desired sequences.}
	\item{Q}{transition matrix for the simulation. Row and column names (\code{"a"}, \code{"c"}, \code{"g"}, and \code{"t"}, although not necessarily in that order) should be provided. If \code{NULL}, a single rate is assumed.}
	\item{rate}{multiplier for \code{Q}, or a vector for \eqn{\Gamma} rate heterogeneity.}
	\item{format}{format of the output object. Can be \code{"DNAbin"}, \code{"phyDat"}, or \code{"matrix"}.}
	\item{...}{optional arguments.}
}
\description{
	Simulates DNA sequence on \code{tree} under the specified model.
}
\details{
	Uses \code{\link{sim.Mk}} internally.
}
\value{
	An object of class \code{"DNAbin"} or \code{"phyDat"}, or a matrix of nucleotides.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\examples{
## load tree
data(mammal.tree)

## rescale edge lengths
mammal.tree$edge.length<-mammal.tree$edge.length/
    max(nodeHeights(mammal.tree))*0.2

## simulate gamma rate heterogeneity
gg<-rgamma(n=100,shape=0.25,rate=0.25)

## simulate DNA sequences with rate heterogeneity
dna<-genSeq(mammal.tree,l=100,rate=gg)
}
\keyword{phylogenetics}
\keyword{simulation}
\keyword{discrete character}
