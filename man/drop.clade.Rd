\name{drop.clade}
\alias{drop.clade}
\title{Drop a clade from a tree}
\usage{
drop.clade(tree, tip)
}
\arguments{
	\item{tree}{object of class \code{"phylo"}.}
	\item{tip}{set of tips in a clade.}
}
\description{
	Internal function for \code{\link{posterior.evolrate}}.
}
\details{
	Function drops the clade containing the species in \code{tip}.

	Probably should not use unless you know what you're doing.
}
\value{
	An object of class \code{"phylo"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{utilities}
