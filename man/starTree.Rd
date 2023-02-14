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
	Revell, L. J. (2012) phytools: An R package for phylogenetic comparative biology (and other things). \emph{Methods Ecol. Evol.}, \bold{3}, 217-223.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\seealso{
	\code{\link{stree}}
}
\keyword{phylogenetics}
\keyword{utilities}
