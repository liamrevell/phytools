\name{add.species.to.genus}
\alias{add.species.to.genus}
\title{Add species to genus on a phylogeny}
\usage{
add.species.to.genus(tree, species, genus=NULL, where=c("root","random"))
}
\arguments{
	\item{tree}{an object of class \code{"phylo"}.}
	\item{species}{string contain species name in the format \code{"Genus_species"} or \code{"Genus species"}.}
	\item{genus}{optional argument containing the genus to which \code{species} is to be attached. If \code{NULL} then \code{genus} will be extracted from \code{species}.}
	\item{where}{location to attach \code{species} to the tree. \code{where="root"} will cause the species to be attached to the MRCA of all members of \code{genus}. \code{where="random"} will cause \code{species} to be attached at random to the subtree descended from the MRCA of all members of \code{genus}.}
}
\description{
	This function adds an additional species to a genus on a phylogeny.
}
\details{
	If \code{genus} contains only one species and \code{where="root"}, then \code{species} will be attached midway along the branch leading to the one species. If \code{where="random"} then \code{species} will be added at a random position along the edge. If \code{genus} cannot be found in the tree, then the original tree is returned and a warning printed. If the tree is not ultrametric, then the resultant tree may not contain branch lengths and a warning will be printed. If \code{genus} is non-monophyletic then \code{species} will be attached to the most inclusive group containing members of \code{genus} and a warning will be printed.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{add.random}}, \code{\link{bind.tip}}
}
\keyword{phylogenetics}
\keyword{utilities}
