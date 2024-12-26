\name{starTree}
\alias{starTree}
\title{Create star phylogeny}
\usage{
starTree(species, branch.lengths=NULL)
}
\arguments{
	\item{species}{a list of species.}
	\item{branch.lengths}{an optional list of branch lengths in the same order as \code{species}.}
}
\description{
	Creates a star phylogeny.
}
\details{
	Creates a star phylogeny with (optionally) user specified branch lengths.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link[ape]{stree}}
}
\keyword{phylogenetics}
\keyword{utilities}
