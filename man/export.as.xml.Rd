\name{export.as.xml}
\alias{export.as.xml}
\title{Export trees & data in XML format}
\usage{
export.as.xml(file, trees, X)
}
\arguments{
	\item{file}{filename for export.}
	\item{trees}{a phylogenetic tree or trees in \code{"phylo"} or \code{"multiPhylo"} format.}
	\item{X}{a matrix of class \code{"DNAbin"} or a matrix with discretely valued non-DNA character data.}
}
\description{
	Exports trees & character data in XML format.
}
\details{
	Can be used to create input file for the program SIMMAP v1.5 (Bollback 2006).
}
\value{
	A file.
}
\references{
	Bollback, J. P. (2006) Stochastic character mapping of discrete traits on phylogenies. \emph{BMC Bioinformatics}, \bold{7}, 88.

	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{make.simmap}}, \code{\link[ape]{read.nexus}}, \code{\link{read.simmap}}, \code{\link{write.simmap}}
}
\keyword{phylogenetics}
\keyword{utilities}
