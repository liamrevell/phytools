\name{writeNexus}
\alias{writeNexus}
\title{Write a tree to file in Nexus format}
\usage{
writeNexus(tree, file="")
}
\arguments{
 	\item{tree}{object of class \code{"phylo"} or \code{"multiPhylo"}.}
	\item{file}{file name for output.}
}
\description{
	This function writes one or multiple phylogenetic trees to file in NEXUS format. Redundant with \code{ape::\link{write.nexus}}.
}
\value{
	Trees written to file.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{write.simmap}}, \code{\link{write.nexus}}
}
\keyword{phylogenetics}
\keyword{input/output}
