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
	Writes one or multiple phylogenetic trees to file in NEXUS format. 
	
	Somewhat redundant with \code{\link[ape]{write.nexus}}.
}
\value{
	Trees written to file.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
  \code{\link{write.simmap}}, \code{\link[ape]{write.nexus}}
}
\keyword{phylogenetics}
\keyword{input/output}
